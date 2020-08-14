from utils.contact_io import load_npz_contacts as lnc
from utils.links import link_features
from utils.misc import split_by_chr, process_regions

import numpy as np
import pandas as pd
import os
import itertools
from scipy.sparse import coo_matrix as coo
from multiprocessing import Pool
from functools import partial

CHROMS = [str(i+1) for i in np.arange(19)] + ['X']

def single_chrom_ABC(pregs,
                     eregs,
                     llinks,
                     contacts,
                     chrom
                    ):
    
    prom_bins = np.mean(pregs[chrom][:,:2],axis = 1).astype('int')
    enh_bins = np.mean(eregs[chrom][:,:2],axis = 1).astype('int')
    
    prom_exp = pregs[chrom][:,-1]
    enh_act = eregs[chrom][:,-1]
    
    prom_bins = (prom_bins/5000).astype('int')
    enh_bins = (enh_bins/5000).astype('int')
    
    #Make sure we have any contact data for our promoter and enhancer bins 
    prom_bins[prom_bins>=contacts[chrom].shape[0]] = 0
    enh_bins[enh_bins>=contacts[chrom].shape[0]] = 0
    
    clipped_contacts = contacts[chrom].todense()[prom_bins,:]
    clipped_contacts = clipped_contacts[:,enh_bins]
    
    print("Clipped contact data for chromosome {}".format(chrom))
    
    ABC_scores = np.multiply(llinks[chrom].todense(), np.repeat(enh_act[None,:], prom_bins.shape[0], axis = 0))
    ABC_scores = np.multiply(ABC_scores, clipped_contacts)
    ABC_scores[prom_exp == 0,:] = 0
    
    ABC_score_sums = np.array(np.sum(ABC_scores, axis = 1))[:,0]
    
    ABC_scores[ABC_score_sums>0,:] = np.divide(ABC_scores[ABC_score_sums>0,:],
                                               np.repeat(ABC_score_sums[ABC_score_sums>0,None],
                                                         ABC_scores.shape[1],
                                                         axis = 1)
                                              ) 
    
    ABC_scores = coo(ABC_scores)
    
    print("Got scores for chromosome {}".format(chrom))
    return ABC_scores.data, np.append(ABC_scores.row[:,None], ABC_scores.col[:,None], axis = 1), chrom

def multi_chrom_ABC(pregs,
                    eregs,
                    llinks,
                    contacts,
                    threshold
                   ):
    
    outlinks = {}
    tot = 0
    fn = partial(single_chrom_ABC,pregs,eregs, llinks, contacts)
    p = Pool()
    t_outs = p.imap(fn, (chrom for chrom in CHROMS))
    for t_out in t_outs:
        scores = t_out[0]
        links = t_out[1]
        chrom = t_out[2]
        tot += np.sum(scores>threshold)
        outlinks[chrom + "_links"] = links[scores>threshold,:]
        outlinks[chrom + "_scores"] = scores[scores>threshold]
        
    p.close()
    p.terminate()
    
    print("TOTAL # of promoter-enhancer LINKS with positive ABC score: {}".format(tot))

    return outlinks

if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='Read in a file or set of files, and return the result.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p","--promoters",
                        help="bed file containing the first set of regions. Should have columns: chrom, start, end",
                        default = "data/raw/promoters/promoter_info.csv",
                        type=str)
    parser.add_argument("-e","--enhancers",
                        help="bed file containing the second set of regions. Should have columns: chrom, start, end",
                        default = "data/raw/enhancers/enhancer_info.csv",
                        type=str)
    parser.add_argument("-s","--strength",
                        help="bed file detailing the strength of k27ac signal at each enhancer",
                        default = "data/raw/enhancers/naive_k27ac.csv",
                        type=str)
    parser.add_argument("-l","--linearlinks",
                        help= ".npz file generated using linear_links.py detailing the associations between enhancers and promoters in terms of linear distance along the chromosome",
                        default = "data/processed/links/linear_promoter_enhancer_links.npz"
                       )
    parser.add_argument("-c","--contacts",
                        help="contact file in .npz format generated using nuc_tools bin_csv function and VC-SQRT normed using norm_contacts.py",
                        default = "data/processed/contacts/normed_SLX-7671_haploid_5kb.npz")
    parser.add_argument("-o", "--outpath",
                        help="path of the output directory for the numpy archive detailing p-e enhancer links and ABC scores",
                        default = "data/processed/links/",
                        type=str)
    parser.add_argument("-r", "--rnaexpression",
                        help="optional path to a bed detailing the expression from the associated genes of each promoter. If provided then ABC links are only calculated for those promoters with positive rna-expression",
                        default = None
                       )
    parser.add_argument("-t", "--threshold",
                        help="threshold of ABC scores to store. Only links with ABC scores > thresh will be stored",
                        default= 0.5,
                        type = float)
    
    args = parser.parse_args()
    
    print("Processing promoters...")
    pregs = pd.read_csv(args.promoters,
                     sep = "\t",
                     na_filter = False,
                     dtype = {"chr":"str",
                              "start":"int64",
                              "end":"int64",
                              "id": "str"
                             }
                    )
    print("Processing enhancers...")
    eregs = pd.read_csv(args.enhancers,
                     sep = "\t",
                     na_filter = False,
                     dtype = {"chr":"str",
                              "start":"int64",
                              "end":"int64",
                              "id": "str"
                             }
                    )
    print("Processing enhancer strengths...")
    estrength = pd.read_csv(args.strength)
    try:
        eregs = eregs.join(estrength)
    except:
        print("Failed to join enhancer strengths dataframe to the enhancer region dataframe. Proceeding with identical enhancer strengths")
        eregs['strength'] = np.ones(eregs.values.shape[0]).astype(float)
        
    if args.rnaexpression is not None:
        rna = pd.read_csv(args.rnaexpression)
        try:
            pregs = pregs.join(rna)
        except:
            print("Failed to join RNA expression dataframe to the promoter region dataframe. Proceeding with all promoters...")
            pregs["expression"] = np.ones(pregs.values.shape[0]).astype(float)
    else:
        pregs["expression"] = np.ones(pregs.values.shape[0]).astype(float)
     
    print("Splitting enhancers and promoters by chromosome")
    
    
    pregs = split_by_chr(pregs)
    eregs = split_by_chr(eregs)
    
    print("Parsing E-P links...")
    llinks = np.load(args.linearlinks, allow_pickle = True)
    llinks = {key: coo((np.ones(llinks[key].shape[0]),
                        (llinks[key][:,0], llinks[key][:,1])
                       ),
                       shape = (pregs[key].shape[0],eregs[key].shape[0])) for key in llinks}
    
    
    print("Loading normed contacts...")
    contacts = np.load(args.contacts, allow_pickle = True)
    contacts = {key: contacts[key][()] for key in contacts}
    
    outlinks = multi_chrom_ABC(pregs, eregs, llinks, contacts, args.threshold)
    
    np.savez(os.path.join(args.outpath,
                          "ABC_EP_links_{}".format(os.path.split(args.contacts)[-1][:-4])
                         ),
             **outlinks
            )
                              
                              
                              
                              
