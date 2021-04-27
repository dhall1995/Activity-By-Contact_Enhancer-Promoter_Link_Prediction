import numpy as np
from .misc import safe_divide
import math
import numpy as np
from numpy import int32
from scipy.sparse import coo_matrix as coo

import straw
CHROMS = [str(i+1) for i in np.arange(19)] + ['X']

def straw_extract_chr_matrix(hic_file,
                             chrom,
                             binsize = 5000,
                             normalisation = 'KR',
                             trim = False
                            ):
    result = straw.straw(normalisation,hic_file,str(chrom), str(chrom),'BP',binsize)
    
    data = np.array(result[2] + result[2])
    goodidxs = ~np.isnan(data)

    if np.sum(goodidxs) == 0:
        print("Unable to get good data (no valid contact pairs) for requested normalisation on chromosome {}, using SQRT_VC instead".format(chrom))
        result = straw.straw('VC_SQRT',hic_file,str(chrom), str(chrom),'BP',binsize)
        data = np.array(result[2] + result[2])
        goodidxs = ~np.isnan(data)
        print("\t{} pairs covered with SQRT_VC normalisation".format(np.sum(goodidxs)))
        
    row  = np.array(result[0] + result[1])
    row = row[goodidxs]
    rowmin = np.min(row)
    rowmax = np.max(row)

    if trim:
        row = (row - rowmin)/binsize
    else:
        row = row/binsize
    row = row.astype('int32')

    col  = np.array(result[1] + result[0])
    col = col[goodidxs]
    if trim:
        col = (col - rowmin)/binsize
    else:
        col = col/binsize 
        
    col = col.astype('int32')

    mycoo = coo((data[goodidxs], (row, col)), shape=(np.max(row)+1, np.max(row)+1))
    
    return mycoo
    
    
from statsmodels import regression
from statsmodels.api import add_constant
from multiprocessing import Pool
from functools import partial

def get_mean_contact_strength_at_offset_(hic_file, offsets, chrom, binsize = 5000, normalisation='KR'):
    tempmat = straw_extract_chr_matrix(hic_file,
                                       chrom,
                                       binsize=binsize,
                                       trim = True,
                                       normalisation=normalisation)
    
    return np.array([np.mean(np.diagonal(tempmat.todense(), o_idx)) for o_idx in offsets])

def infer_gamma(hic_file, binsize =5e3, thresh = 1e6, normalisation = 'KR',
                chroms = CHROMS, verbose = True):
    
    dists = (np.arange(binsize, thresh, binsize)/binsize).astype('int')
    strengths = np.zeros(dists.shape)
    p = Pool()
    fn = partial(get_mean_contact_strength_at_offset_, hic_file, dists, binsize = binsize, normalisation=normalisation)
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
                  chrom,
                  thresh = 5e6,
                  binsize = 5000
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
    ecs[(c_dists >= thresh/(binsize*5))&(c_dists <= thresh/binsize)] = diag*c_dists[(c_dists >= thresh/(binsize*5))&(c_dists <= thresh/binsize)]**gamma
    ecs[(c_dists < thresh/(binsize*5))] = diag*((thresh/(binsize*5))**gamma)

    #rescale expected contact frequencies to lie (theoretically) in the same scale as the actual contacts. For our
    #purposes the maximum theoretical contact frequency we should observe is just the median of the diagonal elements
    #in our true contact matrix. However, as per the methods outlined in the paper, for those e-p pairs closer than
    #1Mb we simply put their expected contact frequency as if they were 1Mb apart. In order to have the same scale, 
    #we don't really care about the 1Mb apart contact frequency, we only require that were a frequency as large as
    #the observed median diagonal frequency then it would be normalised to a score of 100
    ecs *= 100/diag
    
    
    #replace the zeros with inferred/theoretical contact strengths
    try:
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
    except:
        logs = ""
          
    clipped_contacts += ecs
    clipped_contacts[clipped_contacts>100] = 100
    
    return clipped_contacts, logs

def theoretical_contact_strengths(prom_bins,
                                  enh_bins,
                                  gamma,
                                  chrom,
                                  binsize = 5000,
                                  thresh = 5e6
                                 ):
    #infer contact strengths from the theoretical
    #contact strengths given some power law relationship
    c_dists = abs(prom_bins[:,None] - enh_bins[None,:])
    
    #expected contact strengths
    ecs = np.zeros(c_dists.shape)
    #for E-P pairs less than 5 megabase we use the theoretical contact strength
    ecs[(c_dists <= thresh/binsize)&(c_dists>0)] = c_dists[(c_dists <= thresh/binsize)&(c_dists>0)].astype('float')**float(gamma)
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
    
