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

def split_by_chr(df, chrom_col=0, accepted_cols=None):
    if accepted_cols is not None:
        return {k1: np.concatenate([item[None,accepted_cols] for item in list(g1)],axis = 0) for k1,g1 in itertools.groupby(sorted(df.values, key = lambda x:x[chrom_col]),lambda x: x[chrom_col])}
    else:
        return {k1: np.concatenate([item[None,1:] for item in list(g1)],axis = 0) for k1,g1 in itertools.groupby(sorted(df.values, key = lambda x:x[chrom_col]),lambda x: x[chrom_col])}

def process_regions(path, **kwargs):
    regs = pd.read_csv(path,
                     sep = "\t",
                     na_filter = False,
                     dtype = {"chrom":"str",
                              "start":"int64",
                              "end":"int64",
                              "id": "str"
                             },
                       **kwargs
                    )
    regs['chrom'] = regs['chrom'].astype('str')
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
    
    df['chrom'] = df['chrom'].astype('str')
    
    x1 = split_by_chr(df)
    
    return x1

def buffer_vec(vec, buffer = int(1e3), dtype = 'int32'):
    out = np.zeros((vec.shape[0],2)).astype(dtype)
    out[:,0] = vec - buffer
    out[:,1] = vec + buffer
    
    return out

def safe_divide(arr1, arr2):
    out = np.zeros(arr1.shape)
    out[arr2!=0] = np.divide(arr1[arr2!=0], arr2[arr2!=0])
    
    return out

def quantile_normalise(vec, ref):
    vec_idxs = np.argsort(vec)
    ref_idxs = np.argsort(ref)
    
    out = np.zeros(vec.shape)
    out[vec_idxs] = ref[ref_idxs]
    
    return out

def to_bool(vec, shape):
    out = np.array([False]*shape)
    out[vec] = True
    
    return out
