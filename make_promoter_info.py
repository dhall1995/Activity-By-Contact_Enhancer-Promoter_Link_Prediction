import numpy as np
import pandas as pd
import sys
import argparse
import glob
import os

from utils.misc import split_by_chr, buffer_vec, name_chr
from utils.datatracks import DataTrack_rvp as dtrvp

CHROMS = [str(i+1) for i in np.arange(19)] + ['X','Y']


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate candidate regulatory element regions using k27ac ChIP data and TSS sites', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t",
                        "--transcripts",
                        help="bed file containing TSS sites. Should have columns: 'chromosome', 'start', 'end', 'strand','id'. Can optionally include collumns detailing the expression of each promoter in different conditions but they must come AFTER the first five",
                        default = "data/raw/transcripts.tsv",
                        type=str)
    parser.add_argument('-o',
                        '--outpath',
                        default="data/processed/", 
                        help="Path to save candidate element information in")
    parser.add_argument('-k',
                        '--k4',
                        default=None, 
                        nargs = "+",
                        help="Path to save candidate element information in")
    parser.add_argument("-chr","--chromosomes",
                        nargs = "+",
                        help="Chromsomes to retrieve",
                        default = CHROMS,
                        type = str)
    
    args = parser.parse_args()
    
    args.chromosomes = [name_chr(chrom) for chrom in args.chromosomes]
            
    if isinstance(args.k4,str):
        args.k4 = [args.k4]
    
    print("Calculating promoter regions")
    TSS = pd.read_csv(args.transcripts, sep = "\t",
                              na_filter = False,
                              dtype = {"chromosome":"str",
                                       "start":"int64",
                                       "end":"int64",
                                       "strand":"int",
                                       "id": "str"
                                      }
                             )
    if TSS.shape[1]>5:
        exp = TSS.values[:,5:]
        expnames = TSS.columns.values[5:]
    

    TSS = TSS[[item in args.chromosomes for item in TSS['chromosome'].values]]
    
    TSS_vec = np.zeros(TSS.shape[0])
    TSS_vec[TSS['strand'].values==1] = TSS['start'].values[TSS['strand'].values==1]
    TSS_vec[TSS['strand'].values==-1] = TSS['end'].values[TSS['strand'].values==-1]
     
    TSSregs = buffer_vec(TSS_vec, buffer = 1000)
    TSS['promoter start'] = TSSregs[:,0]
    TSS['promoter end'] = TSSregs[:,1]
    
    df = TSS[["chromosome",
              "promoter start",
              "promoter end",
              "id"]].copy()
    df.columns = ['chrom','start','end','id']
    for idx,name in enumerate(expnames):
        df[name] = exp[:,idx]
    
    k4 = [dtrvp('k4').from_bed(peakfile,
                              chrom_col = 0,
                              region_cols = (1,2),
                              allowed_chroms = args.chromosomes
                             ) for peakfile in args.k4]
    
    promk4 = np.zeros(df.shape[0])
    for chrom in df['chrom'].unique():
        idxs = df['chrom'] == chrom
        for track in k4:
            promk4[idxs] += track.stats(chrom, 
                                        df.loc[idxs,['start','end']].values.astype('int32'),
                                        stats_type='coverage')
            
    
    
    df = df.loc[promk4>0]
    df.to_csv(os.path.join(args.outpath,f"promoter_info.csv"),
              sep="\t",
              index = False)
    print("Done. Written promoters to: {}".format(os.path.join(args.outpath,f"promoter_info.csv")))
