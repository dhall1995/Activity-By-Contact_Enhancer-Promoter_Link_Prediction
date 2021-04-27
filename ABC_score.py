from utils.misc import split_by_chr, process_regions,safe_divide,quantile_normalise,make_link_df
from utils.straw_contacts import infer_gamma,theoretical_contact_strengths,clip_contacts,straw_extract_chr_matrix
from utils.datatracks import DataTrack_rvp as dtrvp
from utils.datatracks import DataTrack_bigwig as dtbw

import numpy as np
import pandas as pd
import os
import itertools
from scipy.sparse import coo_matrix as coo
from multiprocessing import Pool
from functools import partial

CHROMS = [str(i+1) for i in np.arange(19)] + ['X','Y']
    
def single_chrom_ABC(pregs,
                     prna,
                     eregs,
                     ek27ac,
                     llinks,
                     contact_path,
                     binsize,
                     normalisation,
                     gamma,
                     distance_thresh,
                     chrom
                    ):
    
    print("loading contact data, chromosome {}".format(chrom))
    contacts = straw_extract_chr_matrix(contact_path, 
                                        chrom,
                                        binsize = binsize,
                                        normalisation = normalisation)
    
    
    prom_bins = np.mean(pregs[chrom],axis = 1).astype('int')
    enh_bins = np.mean(eregs[chrom],axis = 1).astype('int')
    
    ep_dists = abs(prom_bins[:,None] - enh_bins[None,:])
    distance_mask = ep_dists < distance_thresh
    prom_bins = (prom_bins/binsize).astype('int')
    enh_bins = (enh_bins/binsize).astype('int')
    
    enh_act = ek27ac[chrom]
    
    
    if contact_path is not None:
        clipped_contacts, logs = clip_contacts(prom_bins, enh_bins, contacts, gamma,chrom)
        print("Clipped contact data for chromosome {}".format(chrom))
    else:
        clipped_contacts, logs = theoretical_contact_strengths(prom_bins, enh_bins, gamma,chrom)
        print("Using theoretical contact strengths for chromosome {}".format(chrom))
    
    ABC_scores = np.multiply(llinks[chrom].todense(), np.repeat(enh_act[None,:], prom_bins.shape[0], axis = 0))
    ABC_scores = np.multiply(ABC_scores, clipped_contacts)
    ABC_scores[~prna[chrom],:] = 0
    ABC_scores[distance_mask] = 0

    
    ABC_score_sums = np.array(np.sum(ABC_scores, axis = 1))[:,0]
    if np.sum(ABC_score_sums)>0:
        ABC_score_sums[ABC_score_sums==0] = np.percentile(ABC_score_sums[ABC_score_sums>0],1)
    
        ABC_scores[ABC_score_sums>0,:] = np.divide(ABC_scores[ABC_score_sums>0,:],
                                               np.repeat(ABC_score_sums[ABC_score_sums>0,None],
                                                         ABC_scores.shape[1],
                                                         axis = 1)
                                              ) 
    
        ABC_scores = coo(ABC_scores)
    
        print("Got scores for chromosome {}".format(chrom))
        return ABC_scores.data, np.append(ABC_scores.row[:,None], ABC_scores.col[:,None], axis = 1), chrom, logs
    else:
        return np.zeros(1), np.append(np.zeros(1)[:,None], np.zeros(1)[:,None], axis = 1), chrom, logs

def multi_chrom_ABC(pregs,
                    prna,
                    eregs,
                    ek27ac,
                    llinks,
                    contact_path,
                    binsize,
                    normalisation,
                    gamma,
                    distance_thresh,
                    threshold
                   ):
    
    outlinks = {}
    logs = []
    tot = 0
    fn = partial(single_chrom_ABC,pregs, prna, eregs,
                 ek27ac,llinks, contact_path, binsize, normalisation, gamma,distance_thresh)
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
    parser.add_argument("-p",
                        "--promoters",
                        help="bed file containing the first set of regions. Should have columns: chrom, start, end, id",
                        default = "data/processed/promoter_info.csv",
                        type=str)
    parser.add_argument("-re",
                        "--regelements",
                        help="bed file containing the second set of regions. Should have columns: chrom, start, end, regtype",
                        default = "data/processed/candidate_elements.csv",
                        type=str)
    parser.add_argument("-k",
                        "--k27ac",
                        help="H3K27ac ChIPseq file",
                        default = "data/raw/tracks/k27ac/0h_H3K27ac_R1_ERR2588258.Q30.srt.nodup.bam_peaks.narrowPeak",
                        type=str)
    parser.add_argument("-cf",
                        "--chip_format",
                        help="ChIPseq signals for H3K27ac can be provided in either bed format or bigwig format",
                        type = str,
                        default = 'bigwig')
    parser.add_argument("-l",
                        "--linearlinks",
                        help= ".npz file generated using linear_links.py detailing the associations between enhancers and promoters in terms of linear distance along the chromosome",
                        default = "data/processed/links/linear_promoter_regulatory_element_links.npz"
                       )
    parser.add_argument("-c",
                        "--contacts",
                        help="contact file in .hic format generated using juicer from the Lieberman-Aiden Lab",
                        default = None)
    parser.add_argument("-o", 
                        "--outpath",
                        help="path of the output directory for the numpy archive detailing p-e enhancer links and ABC scores",
                        default = "data/processed/links/",
                        type=str)
    parser.add_argument("-b", "--binsize",
                        help="Binsize to bin our contact data at in basepairs",
                        type=str,
                        default = 5000)
    parser.add_argument("-g", 
                        "--gamma",
                        help="power law relationship to use when inferring contact strengths for sparse portions of the contact data",
                        default = None,
                        type = float)
    parser.add_argument("-n", 
                        "--normalisation",
                        help="What matrix normalisation to use for link prediction. Can choose from 'KR', 'NONE',VC','VC_SQRT'",
                        type = str,
                        default = 'KR')
    parser.add_argument("-rc", 
                        "--rnaexpression_column",
                        help="Optional column number detailing rna expression for a given promoter. Provided as a column (>4) of the promoter_info.csv",
                        default = None,
                        type=int 
                       )
    parser.add_argument("-t", 
                        "--threshold",
                        help="threshold of ABC scores to store. Only links with ABC scores > thresh will be stored",
                        default= 0.02,
                        type = float)
    parser.add_argument("-d", 
                        "--distance_threshold",
                        help="Minimum distance for a link to be considered a link",
                        default= 5e3,
                        type = float)
    parser.add_argument("-on",
                        "--outname",
                        help="If no contacts are provided then this provides an identifying tag for the experiments which the links are related to",
                        default = None,
                        type = str
                       )
    parser.add_argument("-tsv", 
                        "--tsv",
                        help="Boolean to determine whether to output the results as a tsv. By default, this tsv is saved in the same directory as the raw links",
                        default = 1,
                        type = int)
    
    args = parser.parse_args()
    
    parser.tsv = bool(args.tsv)
    
    print("Processing k27ac data...")
    if args.chip_format == 'narrowPeak':
        k27ac = dtrvp('k27ac').from_bed(args.k27ac,
                                    chrom_col=0,
                                    region_cols=(1,2),
                                    value_col=6)
    elif args.chip_format == 'bigwig':
        k27ac = dtbw('k27ac').from_bw(args.k27ac)
    else:
        print("Chip format must either be 'narrowPeak' or 'bigwig'")
        raise
    
    print("Processing promoters...")
    pregs = process_regions(args.promoters)
    pregs = split_by_chr(pregs)
    pregs = {chrom: pregs[chrom][:,:2] for chrom in pregs}
    
    pids = process_regions(args.promoters)
    pids = split_by_chr(pids)
    pids = {chrom: pids[chrom][:,2] for chrom in pregs}
    
    if args.rnaexpression_column is not None:
        try:
            print("Filtering promoters by expression")
            rna = process_regions(args.promoters)
            rna = split_by_chr(rna,accepted_cols = [args.rnaexpression_column])
            rna = {chrom: rna[chrom].astype('float')>0 for chrom in rna}
            
            n1 = np.sum([rna[chrom].shape[0] for chrom in rna])
            n2 = np.sum([np.sum(rna[chrom]>0) for chrom in rna])
            print("Filtered from {} to {} promoters".format(n1, n2))
            
            expressed_proms = {chrom: rna[chrom] for chrom in rna}
            
            
            print("For genes with non-zero k27ac signal, filtering out promoters with no k27ac signal...")
            promk27ac = {chrom: k27ac.stats(chrom, 
                                            pregs[chrom], 
                                            stats_type = 'sum')>0 for chrom in pregs}
            
            unique_pids = {chrom: np.unique(pids[chrom][promk27ac[chrom]]) for chrom in pids}
            
            active_proms = {chrom: np.ones(rna[chrom].shape).astype('bool') for chrom in rna}
            for chrom in unique_pids:
                for pid in unique_pids[chrom]:
                    idxs = (pids[chrom] == pid)&(~promk27ac[chrom])
                    active_proms[chrom][idxs] = 0
                    if np.sum(idxs)>0:
                        n = np.sum(idxs)
                        print("Filtering {} promoters from {}".format(n, pid))
                        
            active_proms = {chrom: active_proms[chrom].astype('bool') for chrom in active_proms}
            n3 = np.sum([np.sum((rna[chrom]>0)&(active_proms[chrom]>0)) for chrom in rna])
            print("Filtered from {} to {} promoters".format(n2, n3))
            expressed_proms = {chrom: expressed_proms[chrom][:,0]&active_proms[chrom][:,0] for chrom in expressed_proms}
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
    estrength = {chrom: k27ac.stats(chrom, 
                                    eregs[chrom][:,:2], 
                                    stats_type='sum') for chrom in eregs}  
    estrength = {chrom: quantile_normalise(estrength[chrom],
                                           eregs[chrom][:,2]) for chrom in estrength}
    #remove the reference strengths from the enhancer regions dictionary
    eregs = {chrom: eregs[chrom][:,:2] for chrom in eregs}
    
    for chrom in CHROMS:
        if chrom not in eregs:
            eregs[chrom] = np.zeros((1,2)).astype('int32')
        if chrom not in pregs:
            pregs[chrom] = np.zeros((1,2)).astype('int32')
        if chrom not in estrength:
            estrength[chrom] = np.zeros((1)).astype('int32')
        if chrom not in expressed_proms:
            expressed_proms[chrom] = np.zeros((1)).astype('int32')
    
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
                               args.distance_threshold,
                               args.threshold
                              )
    
    
    if args.outname is not None:
        print("Using provided name for contact file")
        outfilename = os.path.join(args.outpath,
                              "ABC_EP_links_{}_{}".format(os.path.split(args.k27ac)[-1].split(".")[0],
                                                          args.outname)
                             )
        np.savez(outfilename,
                 **outlinks
                )
    elif args.contacts is not None:
        outfilename = os.path.join(args.outpath,
                              "ABC_EP_links_{}_{}".format(os.path.split(args.k27ac)[-1].split(".")[0],
                                                          os.path.split(args.contacts)[-1].split(".")[0])
                             )
        np.savez(outfilename,
                 **outlinks
                )
    else:
        outfilename = os.path.join(args.outpath,
                              "ABC_EP_links_{}".format(os.path.split(args.k27ac)[-1].split(".")[0])
                             )
        np.savez(outfilename,
                 **outlinks
                )
        
                              
    if args.tsv:
        print("Writing output to tsv")
        df = make_link_df(args.regelements, args.promoters, outfilename+".npz", args.binsize)
        df.to_csv(outfilename+".tsv", index = False, sep = "\t")
                              
                              
