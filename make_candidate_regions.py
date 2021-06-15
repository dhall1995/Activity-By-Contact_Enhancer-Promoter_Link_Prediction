from utils.misc import split_by_chr, to_bool, buffer_vec, name_chr
from utils.cython.dtrack_utils import non_overlapping, pairRegionsIntersection
from utils.links import link_features
from utils.datatracks import DataTrack_rvp as dtrvp

import numpy as np
import pandas as pd
from multiprocessing import Pool
from functools import partial
CHROMS = [str(i+1) for i in np.arange(19)] + ['X','Y']


def _single_chrom_candidate_elements_from_tracks(tracks, promoters, gene_ids,chrom):
    print("Calculating promoter regions, chromosome {}..".format(chrom))
    if chrom in promoters:
        promoters[chrom] = promoters[chrom].astype('int32')
    else:
        promoters[chrom] = np.zeros((1,2)).astype('int32')
        
        
    print("Generating unclassified candidate elements for chromosome {} from k27ac peak files".format(chrom))
    regs = np.concatenate([track.regions[chrom] for track in tracks],axis = 0)
    #Temporarily extend k27ac peaks by 250bp in either direction
    regs[:,0] -= 250
    regs[:,1] += 250
    #merge peaks which were closer than 500bp
    regs = non_overlapping(regs.astype('int32'))
    #trim back the edges of the peaks again
    regs[:,0] += 250
    regs[:,1] -= 250
    #get peak lengths
    diffs = np.diff(regs, axis = 1)
    #discard peaks <500bp in length
    regs = np.append(regs[diffs[:,0] > 500,:], promoters[chrom], axis = 0)
    regs = non_overlapping(regs.astype('int32'))
    
    print("Calculating reference k27ac strengths for each candidate element on chromosome {}...".format(chrom))
    strengths = np.mean(np.stack([track.stats(chrom, regs) for track in tracks]), axis = 0)
    print("Unclassified candidate elements generated for chromosome {}. Summary:".format(chrom))
    diagnostics = ["\tChromosome {}".format(chrom),
                   "\t# regions: {}".format(regs.shape[0]),
                   "\tmean strength: {}".format(np.mean(strengths)),
                   "#####################################"]
    print("\n".join(diagnostics))
                                             
    candidate_elements = {}
    candidate_elements['regions'] = {}
    candidate_elements['regions']['unclassified'] = regs
    candidate_elements['strength'] = {}
    candidate_elements['strength']['unclassified'] = strengths
    
    print("Separating promoter, intergenic candidate elements, chromosome {}...".format(chrom))
    candidate_elements['promoter_overlaps'] = pairRegionsIntersection(candidate_elements['regions']['unclassified'],
                                                                      promoters[chrom],
                                                                      allow_partial = True)
    candidate_elements['promoter_overlaps'] = to_bool(candidate_elements['promoter_overlaps'],
                                                      candidate_elements['regions']['unclassified'].shape[0])
    
    conds = {'intergenic':~(candidate_elements['promoter_overlaps']),
             'promoter': candidate_elements['promoter_overlaps']
            }
    
    for cond in ['intergenic','promoter']:
        candidate_elements['regions'][cond] = candidate_elements['regions']['unclassified'][conds[cond],:]
        candidate_elements['strength'][cond] = candidate_elements['strength']['unclassified'][conds[cond]]
    
    conds = ['intergenic','promoter']
    
    print("Linking promoter candidate elements with their respective genes, chromosome {}...".format(chrom))
    
    print("Promoter links chromosome {}...".format(chrom))
    candidate_elements['element_promoter_links'] = link_features(candidate_elements['regions']['promoter'],
                                                                 promoters[chrom])
    
    print("Assigning IDs chromosome {}...".format(chrom))
    candidate_elements['ids'] = {cond: [] for cond in conds}
    candidate_elements['ids']['intergenic'] = ['intergenic' for reg in candidate_elements['regions']['intergenic']]
    
    for idx, reg in enumerate(candidate_elements['regions']['promoter']):
        links = candidate_elements['element_promoter_links']
        links = links[links[:,0]==idx, 1]
        candidate_elements['ids']['promoter'].append("promoter_{}".format(gene_ids[chrom][links[0],0]))
            
    
    candidate_elements['regions']['classified'] = np.concatenate([candidate_elements['regions'][cond] for cond in conds],
                                                                 axis = 0)
    candidate_elements['strength']['classified'] = np.concatenate([candidate_elements['strength'][cond] for cond in conds])
    candidate_elements['ids']['classified'] = np.concatenate([candidate_elements['ids'][cond] for cond in conds])
    
    out = {'regions': candidate_elements['regions']['classified'],
           'strength': candidate_elements['strength']['classified'],
           'ids': candidate_elements['ids']['classified']
          }
    
    return out, chrom

def candidate_elements_from_tracks(tracks,
                                   promoters,
                                   gene_ids,
                                   chroms = CHROMS):
    
    out = {'regions': {chrom: None for chrom in chroms},
           'strength': {chrom: None for chrom in chroms},
           'ids': {chrom: None for chrom in chroms}}

    fn = partial(_single_chrom_candidate_elements_from_tracks,tracks, promoters, gene_ids)
    p = Pool()
    t_outs = p.imap(fn, (chrom for chrom in chroms))
    for t_out in t_outs:
        chrdict = t_out[0]
        chrom = t_out[1]
        out['regions'][chrom] = chrdict['regions']
        out['strength'][chrom] = chrdict['strength']
        out['ids'][chrom] = chrdict['ids']
        
    p.close()
    p.terminate()

    return out
    
if __name__ == "__main__":
    import sys
    import argparse
    import glob
    import os
    
    parser = argparse.ArgumentParser(description='Generate candidate regulatory element regions using k27ac ChIP data and TSS sites', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p",
                        "--promoters",
                        help="bed file containing promoter info. Can be generated using make_promoter_info.py. Should have columns: chrom,start,end,id",
                        default = "data/processed/WT_promoter_info.csv",
                        type=str)
    parser.add_argument("-k",
                        "--k27ac_files", 
                        nargs="+",
                        help="Paths to files containing k27ac ChIP data in narrowPeak format. All narrowPeak files in the folder will be treated as files with which to construct the reference enhancer strengths. These strengths will be used in downstream analysis to quantile normalise k27ac ChIP from individual timepoints. This assumes non-quantiative ChIP data"
                       )
    parser.add_argument('-o',
                        '--outpath',
                        default="data/processed/", 
                        help="Path to save candidate element information in")
    parser.add_argument("-chr","--chromosomes",
                        nargs = "+",
                        help="Chromsomes to retrieve",
                        default = CHROMS,
                        type = str)
    
    args = parser.parse_args()
    args.chromosomes = [name_chr(chrom) for chrom in args.chromosomes]

    # Parse paths
    k27ac_files = args.k27ac_files
        
    proms = pd.read_csv(args.promoters, sep = "\t",
                              na_filter = False,
                              dtype = {"chrom":"str",
                                       "start":"int64",
                                       "end":"int64",
                                       "id":"str"
                                      }
                             )

    proms = proms[[item in CHROMS for item in proms['chrom'].values]]
     
    promoters = split_by_chr(proms,accepted_cols = np.array([1,2]))
    gene_ids = split_by_chr(proms,accepted_cols = np.array([3]))
    
    print("Got promoters and IDs")
    print("Making k27ac tracks...")
    k27ac_tracks = [dtrvp('k27ac').from_bed(path,chrom_col=0,
                                            region_cols=(1,2),value_col=6) for path in k27ac_files]
    
    print("Made k27ac tracks.")
    candidate_elements = candidate_elements_from_tracks(k27ac_tracks, promoters, gene_ids)
    candidate_chroms = {chrom: np.array([chrom for reg in candidate_elements['regions'][chrom]]) for chrom in CHROMS}
    
    
    print("Writing candidate element data to file...")
    df_data = {chrom: np.concatenate([candidate_chroms[chrom][:,None],
                                  candidate_elements['regions'][chrom],
                                  candidate_elements['strength'][chrom][:,None],
                                  candidate_elements['ids'][chrom][:,None]
                                 ], axis = 1) for chrom in CHROMS}
    df_data = np.concatenate([df_data[chrom] for chrom in CHROMS], axis = 0)

    df = pd.DataFrame(data=df_data, columns = ["chrom","start","end","strength","id"])
    
    #df['chrom'] = ['chr'+ item for item in df['chrom']]
    df.to_csv(os.path.join(args.outpath,"candidate_elements.csv"),
              sep="\t",
              index = False)
    print("Done. Written to: {}".format(os.path.join(args.outpath,
                                                     f"candidate_elements.csv")))
