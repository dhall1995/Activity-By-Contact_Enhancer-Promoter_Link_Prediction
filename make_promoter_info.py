import numpy as np
import pandas as pd
import sys
import argparse
import glob
import os

from utils.misc import split_by_chr
CHROMS = [str(i+1) for i in np.arange(19)] + ['X']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate candidate regulatory element regions using k27ac ChIP data and TSS sites', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t","--tss_sites",
                        help="bed file containing TSS sites. Should have columns: chrom, transcript_start, transcript_end, strand, id",
                        default = "data/raw/gene_info/ensembl_TSS_sites_GRCm38_p2.csv",
                        type=str)
    parser.add_argument('-o','--outpath',default="data/processed/", help="Path to save candidate element information in")
    
    args = parser.parse_args()
    
    print("Calculating promoter regions")
    tss_regions = pd.read_csv(args.tss_sites, sep = "\t",
                              na_filter = False,
                              dtype = {"chrom":"str",
                                       "transcript_start":"int64",
                                       "transcript_end":"int64",
                                       "strand":"int64",
                                       "id": "str"
                                      }
                             )
    transcripts = split_by_chr(tss_regions, accepted_cols = np.array([1,2]))
    strands = split_by_chr(tss_regions, accepted_cols = np.array([3]))
    gene_ids = split_by_chr(tss_regions, accepted_cols = np.array([4]))
    
    promoters = {chrom: np.zeros(transcripts[chrom].shape) for chrom in transcripts}
    for chrom in promoters:
        for idx,prom in enumerate(promoters[chrom]):
            if strands[chrom][idx] == 1:
                promoters[chrom][idx,:] = [transcripts[chrom][idx,0]-250, transcripts[chrom][idx,0]+250]
            else:
                promoters[chrom][idx,:] = [transcripts[chrom][idx,1]-250, transcripts[chrom][idx,1]+250]
        
        promoters[chrom] = promoters[chrom].astype('int32')
        transcripts[chrom] = transcripts[chrom].astype('int32')
        
    
    promoter_chroms = {chrom: np.array([chrom for reg in promoters[chrom]]) for chrom in CHROMS}
    
    print("Writing promoter data to file...")
    df_data = {chrom: np.concatenate([promoter_chroms[chrom][:,None],
                                  promoters[chrom],
                                  gene_ids[chrom]
                                 ], axis = 1) for chrom in CHROMS}
    df_data = np.concatenate([df_data[chrom] for chrom in CHROMS], axis = 0)

    df = pd.DataFrame(data=df_data, columns = ["chrom","start","end","id"])
    
    df.to_csv(os.path.join(args.outpath,"promoter_info.csv"),
              sep="\t",
              index = False)
    print("Done. Written to: {}".format(os.path.join(args.outpath,"promoter_info.csv")))