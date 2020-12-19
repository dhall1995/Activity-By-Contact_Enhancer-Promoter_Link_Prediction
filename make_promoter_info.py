import numpy as np
import pandas as pd
import sys
import argparse
import glob
import os

from utils.misc import split_by_chr, buffer_vec
CHROMS = [str(i+1) for i in np.arange(19)] + ['X','Y']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate candidate regulatory element regions using k27ac ChIP data and TSS sites', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t","--tss_sites",
                        help="bed file containing TSS sites. Should have columns: 'Chromosome/scaffold name', 'Transcription start site (TSS)', 'Transcript start', 'Transcript end', 'Gene stable ID'",
                        default = "data/raw/ensembl_TSS_sites_GRCm38_p6.tsv",
                        type=str)
    parser.add_argument('-o','--outpath',default="data/processed/", help="Path to save candidate element information in")
    
    args = parser.parse_args()
    
    print("Calculating promoter regions")
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
    
    df = TSS[["Chromosome/scaffold name","Promoter start","Promoter end","Gene stable ID"]]
    
    df.columns = ['chrom','start','end','id']
    
    df.to_csv(os.path.join(args.outpath,"promoter_info.csv"),
              sep="\t",
              index = False)
    print("Done. Written to: {}".format(os.path.join(args.outpath,"promoter_info.csv")))
