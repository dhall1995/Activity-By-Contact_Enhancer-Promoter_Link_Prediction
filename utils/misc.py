import numpy as np
import pandas as pd
import itertools
import pickle

######################################################################
def save_obj(obj, out_path):
    with open(out_path +'.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(in_path):
    with open(in_path + '.pkl', 'rb') as f:
        return pickle.load(f)

def split_by_chr(df):
    return {k1: np.concatenate([item[None,1:] for item in list(g1)],axis = 0) for k1,g1 in itertools.groupby(sorted(df.values, key = lambda x:x[0]),lambda x: x[0])}

def process_regions(path, **kwargs):
    regs = pd.read_csv(path,
                     sep = "\t",
                     na_filter = False,
                     dtype = {"chr":"str",
                              "start":"int64",
                              "end":"int64",
                              "id": "str"
                             },
                       **kwargs
                    )
    
    return regs

def parse_bed(bed,
              allowed_chroms = [str(i+1) for i in np.arange(20)] + ['X'],
              **kwargs
             ):
    
    df = pd.read_csv(bed,
                     sep = "\t",
                     na_filter = False,
                     **kwargs
                    )[['chrom','start','end']]
    
    x1 = split_by_chr(df)
    
    return x1

def buffer_vec(vec, buffer = int(1e3)):
    out = np.empty((vec.shape[0],2))
    
    out[:,0] = vec - buffer
    out[:,1] = vec + buffer
    
    return out