from utils.misc import split_by_chr, process_regions, safe_divide, quantile_normalise
from utils.straw_contacts import infer_gamma, theoretical_contact_strengths, clip_contacts, straw_extract_chr_matrix
from utils.datatracks import DataTrack_rvp as dtrvp

import numpy as np
import pandas as pd
import os
import itertools
from scipy.sparse import coo_matrix as coo
from multiprocessing import Pool
from functools import partial

CHROMS = [str(i+1) for i in np.arange(19)] + ['X']
    
def single_chrom_ABC(pregs,
                     prna,
                     eregs,
                     ek27ac,
                     llinks,
                     contact_path,
                     binsize,
                     normalisation,
                     gamma,
                     chrom
                    ):
    
    print("loading contact data, chromosome {}".format(chrom))
    contacts = straw_extract_chr_matrix(contact_path, 
                                        chrom,
                                        binsize = binsize,
                                        normalisation = normalisation)
    
    
    prom_bins = np.mean(pregs[chrom],axis = 1).astype('int')
    enh_bins = np.mean(eregs[chrom],axis = 1).astype('int')
    
    enh_act = ek27ac[chrom]
    
    prom_bins = (prom_bins/binsize).astype('int')
    enh_bins = (enh_bins/binsize).astype('int')
    
    if contact_path is not None:
        clipped_contacts, logs = clip_contacts(prom_bins, enh_bins, contacts, gamma,chrom)
        print("Clipped contact data for chromosome {}".format(chrom))
    else:
        clipped_contacts, logs = theoretical_contact_strengths(prom_bins, enh_bins, gamma,chrom)
        print("Using theoretical contact strengths for chromosome {}".format(chrom))
    
    ABC_scores = np.multiply(llinks[chrom].todense(), np.repeat(enh_act[None,:], prom_bins.shape[0], axis = 0))
    ABC_scores = np.multiply(ABC_scores, clipped_contacts)
    ABC_scores[~prna[chrom],:] = 0

    
    ABC_score_sums = np.array(np.sum(ABC_scores, axis = 1))[:,0]
    
    ABC_scores[ABC_score_sums>0,:] = np.divide(ABC_scores[ABC_score_sums>0,:],
                                               np.repeat(ABC_score_sums[ABC_score_sums>0,None],
                                                         ABC_scores.shape[1],
                                                         axis = 1)
                                              ) 
    
    ABC_scores = coo(ABC_scores)
    
    print("Got scores for chromosome {}".format(chrom))
    return ABC_scores.data, np.append(ABC_scores.row[:,None], ABC_scores.col[:,None], axis = 1), chrom, logs

def multi_chrom_ABC(pregs,
                    prna,
                    eregs,
                    ek27ac,
                    llinks,
                    contact_path,
                    binsize,
                    normalisation,
                    gamma,
                    threshold
                   ):
    
    outlinks = {}
    logs = []
    tot = 0
    fn = partial(single_chrom_ABC,pregs, prna, eregs,
                 ek27ac,llinks, contact_path, binsize, normalisation, gamma)
    p = Pool()
    t_outs = p.imap(fn, (chrom for chrom in CHROMS))
    for t_out in t_outs:
        scores = t_out[0]
        links = t_out[1]
        chrom = t_out[2]
        logs.append(t_out[3])
        tot += np.sum(scores>threshold)
        outlinks[chrom + "_links"] = links[scores>threshold,:]
        outlinks[chrom + "_scores"] = scores[scores>threshold]
        
    p.close()
    p.terminate()
    
    print("######################################")
    print("######################################")
    for item in logs:
        print(item)
        
    print("TOTAL # of promoter-enhancer LINKS with positive ABC score: {}".format(tot))

    return outlinks

if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='Calculate Promoter-Regulatory Element links using the ABC model presented in https://www.nature.com/articles/s41588-019-0538-0', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p","--promoters",
                        help="bed file containing the first set of regions. Should have columns: chrom, start, end, id",
                        default = "data/processed/promoter_info.csv",
                        type=str)
    parser.add_argument("-re","--regelements",
                        help="bed file containing the second set of regions. Should have columns: chrom, start, end, regtype",
                        default = "data/processed/candidate_elements.csv",
                        type=str)
    parser.add_argument("-k","--k27ac",
                        help="H3K27ac ChIPseq NarrowPeak file",
                        default = "data/raw/tracks/0h_H3K27ac_R1_ERR2588258.Q30.srt.nodup.bam_peaks.narrowPeak",
                        type=str)
    parser.add_argument("-l","--linearlinks",
                        help= ".npz file generated using linear_links.py detailing the associations between enhancers and promoters in terms of linear distance along the chromosome",
                        default = "data/processed/links/linear_promoter_regulatory_element_links.npz"
                       )
    parser.add_argument("-c","--contacts",
                        help="contact file in .hic format generated using juicer from the Lieberman-Aiden Lab",
                        default = None)
    parser.add_argument("-o", "--outpath",
                        help="path of the output directory for the numpy archive detailing p-e enhancer links and ABC scores",
                        default = "data/processed/links/",
                        type=str)
    parser.add_argument("-b", "--binsize",
                        help="Binsize to bin our contact data at in basepairs",
                        type=str,
                        default = 5000)
    parser.add_argument("-g", "--gamma",
                        help="power law relationship to use when inferring contact strengths for sparse portions of the contact data",
                        default = None,
                        type = float)
    parser.add_argument("-n", "--normalisation",
                        help="What matrix normalisation to use for link prediction. Can choose from 'KR', 'NONE',VC','SQRT_VC'",
                        type = str,
                        default = 'KR')
    parser.add_argument("-r", "--rnaexpression",
                        help="optional path to a bed detailing the expression from the associated genes of each promoter. If provided then ABC links are only calculated for those promoters with positive rna-expression",
                        default = None
                       )
    parser.add_argument("-t", "--threshold",
                        help="threshold of ABC scores to store. Only links with ABC scores > thresh will be stored",
                        default= 0.02,
                        type = float)
    parser.add_argument("-on","--outname",
                        help="If no contacts are provided then this provides an identifying tag for the experiments which the links are related to",
                        default = None,
                        type = str
                       )
    
    args = parser.parse_args()
    
    print("Processing k27ac data...")
    k27ac = dtrvp('k27ac').from_bed(args.k27ac,chrom_col=0,
                                    region_cols=(1,2),value_col=6)
    
    print("Processing promoters...")
    pregs = process_regions(args.promoters)
    pregs = split_by_chr(pregs)
    pregs = {chrom: pregs[chrom][:,:2] for chrom in pregs}
    
    if args.rnaexpression is not None:
        try:
            rna = dtrvp("rna").from_bed(args.rnaexpression,
                                    chrom_col = 0,
                                    region_cols = (1,2),
                                    value_col = 6,
                                    header = None)
            expressed_proms = {chrom: rna.stats(chrom, pregs[chrom])>0 for chrom in pregs}
        except:
            print("Couldn't create RNA datatrack from provided file, using k27ac datatrack to estimate expressed genes")
            promk27ac = {chrom: k27ac.stats(chrom, pregs[chrom], stats_type = 'sum') for chrom in pregs}
            thresh = np.percentile(np.concatenate([promk27ac[chrom] for chrom in promk27ac]),40)
            
            expressed_proms = {chrom: promk27ac[chrom] > thresh for chrom in promk27ac}
            
    else:
        print("No RNA data provided, using k27ac datatrack to estimate expressed genes")
        promk27ac = {chrom: k27ac.stats(chrom, pregs[chrom], stats_type = 'sum') for chrom in pregs}
        thresh = np.percentile(np.concatenate([promk27ac[chrom] for chrom in promk27ac]),40)
            
        expressed_proms = {chrom: promk27ac[chrom] > thresh for chrom in promk27ac}
    
    print("Processing regulatory elements...")
    eregs = pd.read_csv(args.regelements,sep="\t",
                        na_filter = False,
                        dtype = {"chrom":"str",
                                 "start":"int64",
                                 "end":"int64",
                                 "strength":"float",
                                 "id": "str"
                                })
                    
    eregs = split_by_chr(eregs[["chrom","start","end","strength"]])
    
    print("Calculating k27ac strength at enhancer elements")
    estrength = {chrom: k27ac.stats(chrom, eregs[chrom][:,:2]) for chrom in eregs}  
    estrength = {chrom: quantile_normalise(estrength[chrom],eregs[chrom][:,2]) for chrom in estrength}
    #remove the reference strengths from the enhancer regions dictionary
    eregs = {chrom: eregs[chrom][:,:2] for chrom in eregs}
    
    
    print("Parsing E-P links...")
    llinks = np.load(args.linearlinks, allow_pickle = True)
    llinks = {key: coo((np.ones(llinks[key].shape[0]),
                        (llinks[key][:,0], llinks[key][:,1])
                       ),
                       shape = (pregs[key].shape[0],eregs[key].shape[0])) for key in llinks}
    
    
    if args.gamma is None:
        if args.contacts is not None:
            print("No value for gamma provided. Inferring power law relationship from contacts...")
            args.gamma = infer_gamma(args.contacts,
                                     binsize = 5e3,
                                     thresh = 1e6,
                                     normalisation = args.normalisation)
            print("Inferred gamma: {}".format(args.gamma))
        else:
            print("No Gamma or contacts provided. Proceeding with Gamma = -1")
            args.gamma = -1
    
    outlinks = multi_chrom_ABC(pregs,
                               expressed_proms,
                               eregs,
                               estrength,
                               llinks,
                               args.contacts,
                               args.binsize,
                               args.normalisation,
                               args.gamma,
                               args.threshold)
    
    
    if args.outname is not None:
        print("Using provided name for contact file")
        np.savez(os.path.join(args.outpath,
                              "ABC_EP_links_{}_{}".format(os.path.split(args.k27ac)[-1].split(".")[0],
                                                          args.outname)
                             ),
                 **outlinks
                )
    elif args.contacts is not None:
        np.savez(os.path.join(args.outpath,
                              "straw_ABC_EP_links_{}_{}".format(os.path.split(args.k27ac)[-1].split(".")[0],
                                                          os.path.split(args.contacts)[-1].split(".")[0])
                             ),
                 **outlinks
                )
    else:
        np.savez(os.path.join(args.outpath,
                              "ABC_EP_links_{}".format(os.path.split(args.k27ac)[-1].split(".")[0])
                             ),
                 **outlinks
                )
        
                              
                              
                              
                              
