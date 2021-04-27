import numpy as np
import pandas as pd
import itertools
import pickle
from tqdm import tqdm

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

def name_chr(chrom):
    if 'chr' in chrom:
        return chrom[3:]
    else:
        return chrom


def make_link_df(rinfo_path, pinfo_path, link_path, binsize):
    links = np.load(link_path, allow_pickle=True)
    scores = {key.split("_")[0]: links[key] for key in links if "scores" in key}
    links = {key.split("_")[0]: links[key] for key in links if "links" in key}

    rinfo = pd.read_csv(rinfo_path,sep="\t",
                        na_filter = False,
                        dtype = {"chrom":"str",
                                 "start":"int64",
                                 "end":"int64",
                                 "strength":"float",
                                 "id": "str"
                                })
    
    rinfo["id"].loc[rinfo["id"]=='intergenic'] = ["intergenic_{}".format(idx) for idx, inter in enumerate(rinfo["id"].loc[rinfo["id"]=="intergenic"])]
    pinfo = pd.read_csv(pinfo_path,sep="\t",
                        na_filter = False,
                        dtype = {"chrom":"str",
                                 "start":"int64",
                                 "end":"int64",
                                 "id": "str"
                                })
    data = {chrom: [] for chrom in pinfo['chrom'].unique()}
    
    promoter = {'regions': split_by_chr(pinfo, accepted_cols = np.array([1,2])),
            'ids': split_by_chr(pinfo, accepted_cols = np.array(3))
           }
    regel = {'regions':split_by_chr(rinfo, accepted_cols = np.array([1,2])),
            'ids': split_by_chr(rinfo, accepted_cols = np.array(4))}
    

    for chrom in data:
        data[chrom] = np.array([[chrom,
                        promoter['regions'][chrom][link[0],0],promoter['regions'][chrom][link[0],1],
                        promoter['ids'][chrom][link[0]],
                        regel['regions'][chrom][link[1],0],regel['regions'][chrom][link[1],1],
                        regel['ids'][chrom][link[1]],
                        scores[chrom][lidx]
                       ] for lidx, link in enumerate(links[chrom])
                      ])
    
    df = pd.DataFrame(data = np.concatenate([data[chrom] for chrom in data], axis = 0),
                      columns = ['chrom','promoter_start','promoter_end','promoter_id',
                                 'regulatory_element_start','regulatory_element_end','regulatory_element_id',
                                 'ABC_score'])
    df['chrom'] = df['chrom'].astype('category')
    df['promoter_start'] = df['promoter_start'].astype('int64')
    df['promoter_end'] = df['promoter_end'].astype('int64')
    df['promoter_id'] = df['promoter_id'].astype('category')
    df['regulatory_element_start'] = df['regulatory_element_start'].astype('int64')
    df['regulatory_element_end'] = df['regulatory_element_end'].astype('int64')
    df['regulatory_element_id'] = df['regulatory_element_id'].astype('category')
    df['ABC_score'] = df['ABC_score'].astype('float')
    
    df['pbins'] = (df[['promoter_start','promoter_end']].mean(axis = 1)/binsize).astype('int32')
    df['equil_links'] = ["{}-{}-{}".format(item[0],
                                              item[1],
                                              item[2]) for item in df[['pbins',
                                                                       'promoter_id',
                                                                       'regulatory_element_id']].values]
    uids = df['equil_links'].unique()
    vals = []
    
    num = uids.shape[0]
    print("Filtering duplicate links")
    for uid in tqdm(uids, ascii=True, desc = "Processed Equivalent Links", total = num):
        eq_links = df[df['equil_links'] == uid][['chrom','promoter_start','promoter_end','promoter_id',
                                 'regulatory_element_start','regulatory_element_end','regulatory_element_id',
                                 'ABC_score']]
        
        max_score = np.argmax(eq_links['ABC_score'].values)
        
        vals.append(eq_links.values[max_score,:])
        
    newdf = pd.DataFrame(data = np.array(vals), columns = ['chrom','promoter_start','promoter_end','promoter_id',
                                 'regulatory_element_start','regulatory_element_end','regulatory_element_id',
                                 'ABC_score'])
    
    return newdf
