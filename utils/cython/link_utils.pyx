from libc.math cimport abs, sqrt, ceil, floor, log, log2, acos, cos
from numpy cimport ndarray
import numpy as np
from numpy import ones, zeros, int32, float32, uint8, fromstring
from numpy import sort, empty, array, arange, concatenate, searchsorted
from numpy import minimum, maximum,divide, std, mean
from numpy import min as nmin
from numpy import max as nmax
import h5py as h

def link_parent_and_child_multi_regions(ndarray[int, ndim=3] pairs,
                                        ndarray[int, ndim=3] regions,
                                        int cutoff = -1,
                                        allow_partial=False):
    '''
    Given an list of pairs of the form [a,b] where a<b and regions of the form [c,d] we want to work pairwise to
    return the region indices of the regions which overlap with each pair.
     
    Arguments:
    
    - pairs: (N,2) shape array where each row details a pair.
    - regions: (M,2) shape array where each row details a region.
    - exclude: Boolean. If true, then return the indices of pairs which don't overlap with any of the regions.
    - allow_partial: Boolean. If true then include pair indices which only partially overlap with a given region.
                     Here we consider the 'partial' to apply to the pair. That is, if allow_partial == False then
                     we must satisfy either:
                         - the entire pair is included in a region if exclude == False
                         - The entire pair must not be included in a region if exclude == True
                     If allow_partial == True then we must satisfy:
                         - At least part of the pair touches at least one region if exclude == False
                         - At least part of the pair touches no regions if exclude == True
                      
    Returns:
    
    - out: A maximum of an (N*M,2) shape array where each column details a pair group index and a regions group
           index where than pair group and that region group overlap. This is essentially gonna be in COO
           format
    '''
    cdef int i, j, k, np, nr, nl = 0
    cdef int overlap
    cdef int maxnp = len(pairs[0,:,0]) 
    cdef int maxnr = len(regions[0,:,0])
    cdef int npairgroups = len(pairs)
    cdef int nreggroups = len(regions)
    cdef ndarray[int, ndim=2] out = zeros((npairgroups*nreggroups,2), int32)
    
    
    for i in range(npairgroups):
        np = maxnp
        for k in range(maxnp):
            if pairs[i,k,0] == cutoff:
                np = k
                break
        for j in range(nreggroups):
            nr = maxnr
            for k in range(maxnr):
                if regions[j,k,0] == cutoff:
                    nr = k
                    break

            overlap = link_parent_and_child_regions(pairs[i,:np,:],
                                                    regions[j,:nr,:],
                                                    allow_partial = allow_partial,
                                                    bool_out = True)
            if overlap == 1:
                out[nl,0] = i
                out[nl,1] = j
                nl +=1
    
    return out[:nl,:]
                


def link_parent_and_child_regions(ndarray[int, ndim=2] parent_regs,
                                  ndarray[int, ndim=2] child_regs,
                                  allow_partial=False, 
                                  exclude = False,
                                  bool_out = False):
    '''
    Given some pairs and regions we want to return a COO matrix detailing overlaps between the parent
    and child regions.
    
    Arguments:
    
    - pairs: (N,2) shape array where each row details a pair.
    - regions: (M,2) shape array where each row details a region.
    - exclude: Boolean. If true then we retrurn a COO matrix of regions and pairs which don't satisfy our
               constraint.
    - allow_partial: Boolean. If true then include pair indices which only partially overlap with a given region.
                     Here we consider the 'partial' to apply to the pair. That is, if allow_partial == False then
                     we must satisfy either:
                         - the entire pair is included in a region if exclude == False
                         - The entire pair must not be included in a region if exclude == True
                     If allow_partial == True then we must satisfy:
                         - At least part of the pair touches at least one region if exclude == False
                         - At least part of the pair touches no regions if exclude == True
                      
    Returns:
    
    - out: A maximum of an (N*M,2) shape array where each column details a pair group index and a regions group
           index where than pair group and that region group overlap. This is essentially gonna be in COO
           format
    '''
    cdef int i, j, k
    cdef int exc = int(exclude)
    cdef int partial = int(allow_partial)
    cdef int np = len(parent_regs)
    cdef int nr = len(child_regs)
    cdef int bool_o = int(bool_out)
    cdef int nl = 0
    cdef ndarray[int, ndim=2] out = empty((np*nr,2), int32)

    
    for i in range(np):
        for j in range(nr):
            if (child_regs[j,1] > parent_regs[i,0]) and (parent_regs[i,1] > child_regs[j,0]):
                #Pair is at least partially overlapping with the region
                if partial:
                    if bool_o:
                        return 1
                    out[nl,0] = i
                    out[nl,1] = j
                    nl += 1
                    continue
                elif (child_regs[j,1] <= parent_regs[i,1]) and (child_regs[j,0] >= parent_regs[i,0]):
                    if exc == 0:
                        if bool_o:
                            return 1
                        out[nl,0] = i
                        out[nl,1] = j
                        nl += 1
            elif exc == 1:
                out[nl,0] = i
                out[nl,1] = j
                nl+=1
    if bool_o:
        return 0
    
    return out[:nl,:]