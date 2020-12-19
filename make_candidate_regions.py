from utils.misc import split_by_chr, to_bool, buffer_vec
from utils.cython.dtrack_utils import non_overlapping, pairRegionsIntersection
from utils.links import link_features
from utils.datatracks import DataTrack_rvp as dtrvp

import numpy as np
import pandas as pd
from multiprocessing import Pool
from functools import partial
CHROMS = [str(i+1) for i in np.arange(19)] + ['X','Y']


def _single_chrom_candidate_elements_from_tracks(tracks, promoters, transcripts, gene_ids,chrom):
    print("Calculating promoter regions, chromosome {}..".format(chrom))
    promoters[chrom] = promoters[chrom].astype('int32')
    transcripts[chrom] = transcripts[chrom].astype('int32')
        
        
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
    #discard peaks <300bp in length
    regs = np.append(regs[diffs[:,0] > 300,:], promoters[chrom], axis = 0)
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
    
    print("Separating promoter, genic and intergenic candidate elements, chromosome {}...".format(chrom))
    candidate_elements['promoter_overlaps'] = pairRegionsIntersection(candidate_elements['regions']['unclassified'],
                                                                      promoters[chrom],
                                                                      allow_partial = True)
    candidate_elements['promoter_overlaps'] = to_bool(candidate_elements['promoter_overlaps'],
                                                      candidate_elements['regions']['unclassified'].shape[0])
    candidate_elements['transcript_overlaps'] = pairRegionsIntersection(candidate_elements['regions']['unclassified'],
                                                                        transcripts[chrom],
                                                                        allow_partial = True)
    candidate_elements['transcript_overlaps'] = to_bool(candidate_elements['transcript_overlaps'],
                                                              candidate_elements['regions']['unclassified'].shape[0])
    
    conds = {'intergenic':~(candidate_elements['promoter_overlaps']|candidate_elements['transcript_overlaps']),
             'promoter': candidate_elements['promoter_overlaps'],
             'genic': candidate_elements['transcript_overlaps']&~candidate_elements['promoter_overlaps']
            }
    
    for cond in ['intergenic','genic','promoter']:
        candidate_elements['regions'][cond] = candidate_elements['regions']['unclassified'][conds[cond],:]
        candidate_elements['strength'][cond] = candidate_elements['strength']['unclassified'][conds[cond]]
    
    conds = ['intergenic','genic','promoter']
    
    print("Linking promoter and genic candidate elements with their respective genes, chromosome {}...".format(chrom))
    
    print("Genic links chromosome {}...".format(chrom))
    candidate_elements['element_transcript_links'] = link_features(candidate_elements['regions']['genic'],
                                                                   transcripts[chrom])
    print("Promoter links chromosome {}...".format(chrom))
    candidate_elements['element_promoter_links'] = link_features(candidate_elements['regions']['promoter'],
                                                                 promoters[chrom])
    
    print("Assigning IDs chromosome {}...".format(chrom))
    candidate_elements['ids'] = {cond: [] for cond in conds}
    candidate_elements['ids']['intergenic'] = ['intergenic' for reg in candidate_elements['regions']['intergenic']]
    
    for idx, reg in enumerate(candidate_elements['regions']['genic']):
        links = candidate_elements['element_transcript_links']
        links = links[links[:,0] == idx,1]
        candidate_elements['ids']['genic'].append("genic_{}".format(gene_ids[chrom][links[0],0]))
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
                                   transcripts,
                                   gene_ids,
                                   chroms = CHROMS):
    
    out = {'regions': {chrom: None for chrom in chroms},
           'strength': {chrom: None for chrom in chroms},
           'ids': {chrom: None for chrom in chroms}}

    fn = partial(_single_chrom_candidate_elements_from_tracks,tracks, promoters, transcripts, gene_ids)
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
    parser.add_argument("-t","--tss_sites",
                        help="bed file containing TSS sites. Should have columns: chrom, transcript_start, transcript_end, strand, id",
                        default = "data/raw/ensembl_TSS_sites_GRCm38_p6.tsv",
                        type=str)
    parser.add_argument("-k","--k27ac_path", nargs="+",
                        help="Path to folder containing k27ac ChIP data in narrowPeak format. All narrowPeak files in the folder will be treated as files with which to construct the reference enhancer strengths. These strengths will be used in downstream analysis to quantile normalise k27ac ChIP from individual timepoints. This assumes non-quantiative ChIP data"
                       )
    parser.add_argument('-e', '--extension', default='narrowPeak', help='File extension to filter by.')
    parser.add_argument('-o','--outpath',default="data/processed/", help="Path to save candidate element information in")
    
    args = parser.parse_args()

    # Parse paths
    full_paths = [os.path.join(os.getcwd(), path) for path in args.k27ac_path]
    k27ac_files = set()
    for path in full_paths:
        if os.path.isfile(path):
            k27ac_files.add(path)
        else:
            k27ac_files |= set(glob.glob(path + '/*' + args.extension))
            
    k27ac_files = list(k27ac_files)
        
    TSS = pd.read_csv(args.tss_sites, sep = "\t",
                              na_filter = False,
                              dtype = {"Chromosome/scaffold name":"str",
                                       "Transcription start site (TSS)":"int64",
                                       "Transcript start":"int64",
                                       "Transcript end":"int64",
                                       "Gene stable ID": "str"
                                      }
                             )
    TSS = TSS[[item in CHROMS for item in TSS['Chromosome/scaffold name'].values]]
     
    TSSregs = buffer_vec(TSS["Transcription start site (TSS)"].values, buffer = 500)
    TSS['Promoter start'] = TSSregs[:,0]
    TSS['Promoter end'] = TSSregs[:,1]
    
    promoters = split_by_chr(TSS,accepted_cols = np.array([5,6]))
    transcripts = split_by_chr(TSS, accepted_cols = np.array([2,3]))
    gene_ids = split_by_chr(TSS, accepted_cols = np.array([4]))
    
    print("Got transcripts, strands and IDs")
    print("Making k27ac tracks...")
    k27ac_tracks = [dtrvp('k27ac').from_bed(path,chrom_col=0,
                                            region_cols=(1,2),value_col=6) for path in k27ac_files]
    
    print("Made k27ac tracks.")
    candidate_elements = candidate_elements_from_tracks(k27ac_tracks, promoters, transcripts, gene_ids)
    candidate_chroms = {chrom: np.array([chrom for reg in candidate_elements['regions'][chrom]]) for chrom in CHROMS}
    
    
    print("Writing candidate element data to file...")
    df_data = {chrom: np.concatenate([candidate_chroms[chrom][:,None],
                                  candidate_elements['regions'][chrom],
                                  candidate_elements['strength'][chrom][:,None],
                                  candidate_elements['ids'][chrom][:,None]
                                 ], axis = 1) for chrom in CHROMS}
    df_data = np.concatenate([df_data[chrom] for chrom in CHROMS], axis = 0)

    df = pd.DataFrame(data=df_data, columns = ["chrom","start","end","strength","id"])
    
    df.to_csv(os.path.join(args.outpath,"candidate_elements.csv"),
              sep="\t",
              index = False)
    print("Done. Written to: {}".format(os.path.join(args.outpath,"candidate_elements.csv")))
