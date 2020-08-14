from utils.misc import parse_bed, buffer_vec, load_obj, save_obj
from utils.links import link_features

from multiprocessing import Pool
from functools import partial

import numpy as np
import pandas as pd
import itertools
import os
import glob
from scipy.sparse import coo_matrix as coo


def single_chrom_links(regs1,
                       regs2,
                       thresh,
                       chrom):
    
    vec1 = buffer_vec(np.mean(regs1[chrom], axis = 1), buffer = thresh)
    vec2 = buffer_vec(np.mean(regs2[chrom], axis = 1), buffer = thresh)
    
    all_links = link_features(vec1, vec2)
    
    print("Done Chromosome {}".format(chrom))
    return all_links, chrom
    


def multi_chrom_links(regs1,
                      regs2,
                      thresh,
                      featname1,
                      featname2):
    
    outlinks = {}
    tot = 0
    fn = partial(single_chrom_links,regs1,regs2, thresh)
    p = Pool()
    t_outs = p.imap(fn, (chrom for chrom in regs1))
    for t_out in t_outs:
        links = t_out[0]
        chrom = t_out[1]
        tot += links.shape[0]
        outlinks[chrom] = links
    p.close()
    p.terminate()
    
    print("TOTAL # of linear {}-{} LINKS: {}".format(featname1, featname2, tot))

    return outlinks

if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='Read in a file or set of files, and return the result.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r1","--regions1",
                        default = "data/raw/promoters/promoter_info.csv",
                        help="bed file containing the first set of regions. Should have columns: chr, start, end",
                        type=str)
    parser.add_argument("-r2","--regions2",
                        default = "data/raw/enhancers/enhancer_info.csv",
                        help="bed file containing the second set of regions. Should have columns: chr, start, end",
                        type=str)
    parser.add_argument("-o", "--outpath",
                        help="path of the output directory for the numpy archive detailing the links between regions",
                        default = "data/processed/links/",
                        type=str)
    parser.add_argument("-r1n", "--region1names",
                        default = "promoter",
                        help="name of the first regions",
                        type = str)
    parser.add_argument("-r2n", "--region2names",
                        default = "enhancer",
                        help="name of the second regions",
                        type = str)
    parser.add_argument("-t", "--distance_thresh",
                        default = 5000,
                        help = "distance threshold to consider links, in kb. defaults to 5000 (5Mb)",
                        type = int)
    
    args = parser.parse_args()
    
    print("Parsing regions...")
    regs1 = parse_bed(args.regions1)
    regs2 = parse_bed(args.regions2)
        
    
    
    print("####################################################")
    print("Linking {}s to {}s via linear separation".format(args.region1names, args.region2names))
    ep_links = multi_chrom_links(regs1,
                                 regs2,
                                 args.distance_thresh*1000,
                                 args.region1names,
                                 args.region2names
                                )
    
    print("Done!")
    print("####################################################")
        
        
    np.savez(os.path.join(args.outpath,"linear_{}_{}_links".format(args.region1names,
                                                               args.region2names)), **ep_links)
                        
