import numpy as np
from .misc import safe_divide
from scipy.sparse import coo_matrix
from scipy import sparse
import pandas as pd
import math
import numpy as np
from numpy import int32
import itertools
from scipy.sparse import coo_matrix as coo

CHR_KEY_SEP = ' '
DEFAULT_BIN_SIZE = 5
DEFAULT_MIN_BINS = 2
MIN_TRANS_COUNT = 5
CHROMS = ["chr"+str(i+1) for i in np.arange(19)] + ['chrX']

def load_npz_contacts(file_path, 
                      store_sparse=False,
                      display_counts=False,
                      normalize = False,
                      cut_centromeres = True,
                      chrompairs = [pair for pair in itertools.permutations(CHROMS,2)] + [(chrom, chrom) for chrom in CHROMS]
                     ):
    '''
    Utility function to load a .npz file containing contact information from a Hi-C experiment. 
    
    Arguments:
    
    - file_path: A .npz file generated using the nuc_tools ncc_bin tool. The function assumes a
                 File of this format
    - store_sparse: Boolean determining whether to return the contact matrices in sparse format
    - display_counts: Boolean determining whether to display summary plots of Hi-C counts
    - normalize: Boolean determining whether to normalise all matrix elements to lie between zero or one.
                 If False then raw contact counts are returned instead
    - cut_centromeres: Boolean determining whether to cut out the centromeres from the beginning of each
                       chromosome. Since the centromeres contain repetitive elements, they can't currently
                       be mapped by Hi-C so these rows and columns should be void of Hi-C contacts. This 
                       does affect indexing later on but other functions in this package should accommodate
                       for that
                       
    Returns:
    
    - bin_size: The size of each chromatin bin in basepairs.
    - chromo_limits: Dictionary detailing the start and end basepairs of each chromosome
                     in the contact dictionary. NOTE: the chromosome limits are inclusive
                     i.e. for each CHR_A we should have chromo_limits[CHR_A] = (start_A,
                     end_A) where all basepairs b on this chromsome satisfy:
                                 start_A <= b <= end_A
    - contacts: Dictionary of matrices detailing contacts between chromosome pairs
    '''
    for idx, pair in enumerate(chrompairs):
        c1, c2 = pair
        if "chr" not in c1:
            c1 = "chr"+c1
        if "chr" not in c2:
            c2 = "chr"+c2
        
        chrompairs[idx] = (c1,c2)
        
    allchroms = list(set([pair[0] for pair in chrompairs] + [pair[1] for pair in chrompairs]))
            
    file_dict = np.load(file_path, allow_pickle=True, encoding = 'bytes')
  
    chromo_limits = {}
    contacts = {}
    bin_size, min_bins = file_dict['params']
    bin_size = int(bin_size*1e3)
  
    chromo_hists = {}
    cis_chromo_hists = {}

    pair_keys = [key for key in file_dict.keys() if "cdata" in key and (key.split(CHR_KEY_SEP)[0],
                                                                        key.split(CHR_KEY_SEP)[1]) in chrompairs]
    nonpair_keys = [key for key in file_dict.keys() if (CHR_KEY_SEP not in key) and (key != 'params') and key in allchroms]
  
    for key in nonpair_keys:
        offset, count = file_dict[key]
        chromo_limits[key] = offset*bin_size, (offset+count)*bin_size
        chromo_hists[key] = np.zeros(count)
        cis_chromo_hists[key] = np.zeros(count)

    maxc = 1
    if normalize:
        for key in sorted(pair_keys):
            maxc = np.maximum(maxc, np.max(file_dict[key]))

    for key in sorted(pair_keys):
        chr_a, chr_b, _ = key.split(CHR_KEY_SEP)
        
        shape = file_dict[chr_a + CHR_KEY_SEP + chr_b + CHR_KEY_SEP + "shape"]
        mtype = "CSR"
        try:
            indices = file_dict[chr_a + CHR_KEY_SEP + chr_b + CHR_KEY_SEP + "ind"]
            indptr = file_dict[chr_a + CHR_KEY_SEP + chr_b + CHR_KEY_SEP + "indptr"]
        except:
            mtype = "COO"
            row = file_dict[chr_a + CHR_KEY_SEP + chr_b + CHR_KEY_SEP + "row"]
            col = file_dict[chr_a + CHR_KEY_SEP + chr_b + CHR_KEY_SEP + "col"]

        if mtype == "CSR":
            mat = sparse.csr_matrix((file_dict[key]/maxc, indices, indptr), shape = shape)
        else:
            mat = sparse.coo_matrix((file_dict[key]/maxc, (row, col)), shape = shape)

        if not store_sparse:
            mat = mat.toarray()
          
        if chr_a == chr_b:
            a, b = mat.shape
            cols = np.arange(a-1)
            rows = cols-1

            if not np.all(mat[rows, cols] == mat[cols, rows]): # Not symmetric
                mat += mat.T
        contacts[(chr_a, chr_b)] = mat  
    
    #Chromosomes in our dataset
    chroms = chromo_limits.keys()
    if cut_centromeres:
    #Exclude centromeres of chromosomes where we don't have any contact data
        for chrom in chroms:
            chrmax = chromo_limits[chrom][-1]
            temp = contacts[(chrom, chrom)].indices
            chromo_limits[chrom] = (bin_size*np.min(temp[temp>0]), chrmax)
    
    for pair in contacts:
        s0, s1 = int(chromo_limits[pair[0]][0]/bin_size), int(chromo_limits[pair[1]][0]/bin_size)
        try:
            contacts[pair] = contacts[pair][s0:,s1:]
        except:
            contacts[pair] = contacts[pair].tocsr()[s0:, s1:].tocoo()           
    
    if display_counts:
        # A simple 1D overview of count densities
 
        from matplotlib import pyplot as plt

        for chr_a, chr_b in contacts:
            mat = contacts[(chr_a, chr_b)]
            chromo_hists[chr_a] += mat.sum(axis=1)
            chromo_hists[chr_b] += mat.sum(axis=0)
 
            if chr_a == chr_b:
                cis_chromo_hists[chr_a] += mat.sum(axis=1)
                cis_chromo_hists[chr_b] += mat.sum(axis=0)
    
        all_sums = np.concatenate([chromo_hists[ch] for ch in chromo_hists])
        cis_sums = np.concatenate([cis_chromo_hists[ch] for ch in chromo_hists])
 
        fig, ax = plt.subplots()
 
        hist, edges = np.histogram(all_sums, bins=25, normed=False, range=(0, 500))
        ax.plot(edges[1:], hist, color='#0080FF', alpha=0.5, label='Whole genome (median=%d)' % np.median(all_sums))

        hist, edges = np.histogram(cis_sums, bins=25, normed=False, range=(0, 500))
        ax.plot(edges[1:], hist, color='#FF4000', alpha=0.5, label='Intra-chromo/contig (median=%d)' % np.median(cis_sums))
 
        ax.set_xlabel('Number of Hi-C RE fragment ends (%d kb region)' % (bin_size/1e3))
        ax.set_ylabel('Count')
 
        ax.legend()
 
        plt.show()

    return bin_size, chromo_limits, contacts


def save_npz_contacts(contacts, outpath, 
                      binsize = DEFAULT_BIN_SIZE,
                      minbincount= DEFAULT_MIN_BINS,
                      chr_key_sep = CHR_KEY_SEP):
    out = {}
    for chromo_key in contacts:
        chr_a, chr_b = chromo_key
        if "chr" not in chr_a:
            chr_a = "chr" + chr_a
        if "chr" not in chr_b:
            chr_b = "chr" + chr_b
        data_key = chr_a + chr_key_sep + chr_b + chr_key_sep + "cdata"    
      
        if chr_a == chr_b:
            indices_key = chr_a + chr_key_sep + chr_b + chr_key_sep + "ind"
            indptr_key = chr_a + chr_key_sep + chr_b + chr_key_sep + "indptr"
            shape_key = chr_a + chr_key_sep + chr_b + chr_key_sep + "shape"

            cdata = contacts[chromo_key].tocsr()
            out[data_key] = cdata.data.astype('float')
            out[indices_key] = cdata.indices.astype('int32')
            out[indptr_key] = cdata.indptr.astype('int32')
            out[shape_key] = cdata.shape

        else:
            row_key = chr_a + chr_key_sep + chr_b + chr_key_sep + "row"
            col_key = chr_a + chr_key_sep + chr_b + chr_key_sep + "col"
            shape_key = chr_a + chr_key_sep + chr_b + chr_key_sep + "shape"

            cdata = contacts[chromo_key]
            out[data_key] = cdata.data.astype('float')
            out[row_key] = cdata.row.astype('int32')
            out[col_key] = cdata.col.astype('int32')
            out[shape_key] = cdata.shape
        
        # Store bin offsets and spans
        out[chr_a] = np.array([0, cdata.shape[0]])
        out[chr_b] = np.array([0, cdata.shape[1]])
    
    out['params'] = np.array([binsize, minbincount])
        
     
    np.savez_compressed(outpath, 
                        **out)
                      

def row_col_sums(contacts, chroms = CHROMS):
    chromshapes = {chrom: contacts[(chrom,chrom)].shape[0] for chrom in chroms}
    
    sums = {chrom: np.zeros((chromshapes[chrom],1)) for chrom in chroms}

    for key in contacts:
        c1,c2 = key
        if c1 in chroms and c2 in chroms:
            sums[c1] += np.sum(contacts[key],axis =1)
            sums[c2] += np.sum(contacts[key],axis = 0).T
        if c1 not in chroms:
            sums[c2] += np.sum(contacts[key],axis =0).T
        if c2 not in chroms:
            sums[c1] += np.sum(contacts[key],axis = 1) 

    return sums

def vc_sqrt_norm(contact_path, sums, key):
    c1,c2 = key
    print('Starting chromosome pair {}-{}'.format(c1,c2))
    _, _, contacts = load_npz_contacts(contact_path,store_sparse=True,
                                      display_counts=False,
                                      normalize = False,
                                      cut_centromeres = False,
                                      chrompairs = [key]
                                     )
    contacts = contacts[key].tocoo()
    
    norms = np.dot(sums[c1], sums[c2].T)
    #goodcoverage = {c: np.where(np.array(sums[c]>0)[0,:])[0] for c in [c1,c2]}
    #norms[np.ix_(goodcoverage[c1],goodcoverage[c2])] = np.sqrt(1./norms[np.ix_(goodcoverage[c1],goodcoverage[c2])])
    norms[norms>0] = np.sqrt(1./norms[norms>0])
    
    print("Caculated row-sum coverage norms, chromosome pair {}-{}".format(c1,c2))
    newdata = np.multiply(contacts.data, [norms[contacts.row[idx],
                                                contacts.col[idx]] for idx in np.arange(contacts.row.shape[0])])

    
    print("Adjusted raw data, chromosome pair {}-{}. Converting to COO format".format(c1,c2))
    mat = coo((newdata, (contacts.row, contacts.col)), shape=contacts.shape)
    print("Converted VC normed data to COO format, chromosome pair {}-{}".format(c1,c2))
    print("Done chromosome pair {}-{}".format(c1,c2))
    return key, mat


from statsmodels import regression
from statsmodels.api import add_constant
from multiprocessing import Pool
from functools import partial

def get_mean_contact_strength_at_offset_(contact_path, offsets, chrom):
    _, _, contacts = load_npz_contacts(contact_path,store_sparse=True,
                                      display_counts=False,
                                      normalize = False,
                                      cut_centromeres = False,
                                      chrompairs = [(chrom,chrom)]
                                     )
    tempmat = contacts[(chrom,chrom)].tocoo().todense()
    return np.array([np.mean(np.diagonal(tempmat, o_idx)) for o_idx in offsets])

def infer_gamma(contact_path, binsize =5e3, thresh = 1e6, chroms = CHROMS, verbose = True):
    
    dists = (np.arange(binsize, thresh, binsize)/binsize).astype('int')
    strengths = np.zeros(dists.shape)
    p = Pool()
    fn = partial(get_mean_contact_strength_at_offset_, contact_path, dists)
    t_outs = p.imap(fn, (chrom for chrom in chroms))
    for t_out in t_outs:
        strengths += (1/len(chroms))*t_out
    
    X = np.log(dists)[1:]
    y = np.log(strengths)[1:]
    
    # Running the linear regression
    model = regression.linear_model.OLS(y, add_constant(X.reshape(-1,1))).fit()
    a = model.params[0]
    b = model.params[1]
    
    if verbose:
        print("Gamma regression summary:")
        print(model.summary())
        
    return b

def clip_contacts(prom_bins, 
                  enh_bins,
                  chr_contacts,
                  gamma,
                  chrom
                 ):
                  
    #Make sure we have any contact data for our promoter and enhancer bins 
    prom_bins[prom_bins>=chr_contacts.shape[0]] = 0
    enh_bins[enh_bins>=chr_contacts.shape[0]] = 0
    
    #row-wise rescale contact data to lie between 0-100
    densecontacts = np.array(chr_contacts.todense())
    rowmax = np.max(densecontacts, axis = 1)
    densecontacts = 100*safe_divide(densecontacts, np.repeat(rowmax[:,None], densecontacts.shape[1], axis = 1))
    
    #clip contacts to just contain promoter rows and enhancer collumns
    clipped_contacts = densecontacts[prom_bins,:]
    clipped_contacts = clipped_contacts[:,enh_bins]
    
    #For relatively sparse data we want to infer contact strengths from the theoretical
    #contact strengths given some power law relationship
    c_dists = abs(prom_bins[:,None] - enh_bins[None,:])
    
    #expected contact strengths
    ecs = np.zeros(c_dists.shape)
    diag = np.diagonal(chr_contacts.todense(),0)
    diag = np.median(diag[diag>0])
    ecs[(c_dists >= 200)&(c_dists <= 1000)] = diag*c_dists[(c_dists >= 200)&(c_dists <= 1000)]**gamma
    ecs[(c_dists < 200)] = diag*(200**gamma)

    #rescale expected contact frequencies to lie (theoretically) in the same scale as the actual contacts. For our
    #purposes the maximum theoretical contact frequency we should observe is just the median of the diagonal elements
    #in our true contact matrix. However, as per the methods outlined in the paper, for those e-p pairs closer than
    #1Mb we simply put their expected contact frequency as if they were 1Mb apart. In order to have the same scale, 
    #we don't really care about the 1Mb apart contact frequency, we only require that were a frequency as large as
    #the observed median diagonal frequency then it would be normalised to a score of 100
    ecs *= 100/diag
    
    
    #replace the zeros with inferred/theoretical contact strengths
    print("Adjusting contact strength for sparse portions of contact data in chromosome {}".format(chrom))
    outprint = ["Summary stats for clipped contacts chromosome {}:".format(chrom),
                "Number of zeros: {}".format(np.sum(clipped_contacts == 0)),
                "Median nonzero value: {}".format(np.median(clipped_contacts[clipped_contacts>0])),
                "Mean nonzero value: {}".format(np.mean(clipped_contacts[clipped_contacts>0])),
                "Maximum value: {}".format(np.max(clipped_contacts)),
                "Minimum nonzero value {}".format(np.min(clipped_contacts[clipped_contacts >0])),
                "######################################",
                "Summary stats for expected contact frequencies chromosome {}:".format(chrom),
                "Number of zeros: {}".format(np.sum(ecs == 0)),
                "Median nonzero value: {}".format(np.median(ecs[ecs>0])),
                "Mean nonzero value: {}".format(np.mean(ecs[ecs>0])),
                "Maximum value: {}".format(np.max(ecs)),
                "Minimum nonzero value {}".format(np.min(ecs[ecs >0])),
                "######################################",
                "######################################"
               ]
    logs = "\n".join(outprint)
          
    clipped_contacts += ecs
    clipped_contacts[clipped_contacts>100] = 100
    
    return clipped_contacts, logs

def theoretical_contact_strengths(prom_bins,
                                  enh_bins,
                                  gamma,
                                  chrom
                                 ):
    #infer contact strengths from the theoretical
    #contact strengths given some power law relationship
    c_dists = abs(prom_bins[:,None] - enh_bins[None,:])
    
    #expected contact strengths
    ecs = np.zeros(c_dists.shape)
    #for E-P pairs less than 5 megabase we use the theoretical contact strength
    ecs[(c_dists <= 1000)&(c_dists>0)] = c_dists[(c_dists <= 1000)&(c_dists>0)].astype('float')**float(gamma)
    ecs[c_dists == 0] = 1
    
    #rescale expected contact frequencies to lie in the range [0,100]
    ecs *= 100
    
    outprint = ["Summary stats for expected contact frequencies chromosome {}:".format(chrom),
                "Number of zeros: {}".format(np.sum(ecs == 0)),
                "Median nonzero value: {}".format(np.median(ecs[ecs>0])),
                "Mean nonzero value: {}".format(np.mean(ecs[ecs>0])),
                "Maximum value: {}".format(np.max(ecs)),
                "Minimum nonzero value {}".format(np.min(ecs[ecs >0])),
                "######################################",
                "######################################"
               ]
    logs = "\n".join(outprint)
    
    return ecs, logs
    