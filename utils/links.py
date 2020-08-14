from .cython.link_utils import link_parent_and_child_regions, link_parent_and_child_multi_regions
import numpy as np

def link_features(parent_regions,
                  child_regions):
    '''
    Given some list of regions associated with each parent and child, this returns an (N,2) shape
    array detailing the links between parents and children. 
    Arguments:
    parent_regions: A list of (M,2) shape arrays (one per parent) which need not be the same length
                    but which detail the regions associated with each parent.
    child_regions: A list of (M,2) shape arrays (one per child) which need not be the same length
                   but which detail the regions associated with each child.
    
    Returns:
    Links: (N,2) shape array (essentially a COO format sparse matrix) where each row details 
           a link between a parent and a child. 
    '''
    if len(parent_regions) == 0 or len(child_regions) == 0:
        return np.empty((0,2)).astype('int32')
    
    biggest_p = np.max([region.shape[0] for region in parent_regions])
    biggest_c = np.max([region.shape[0] for region in child_regions])
    
    smallest_p = np.min([np.min(region) for region in parent_regions])
    smallest_c = np.min([np.min(region) for region in child_regions])
    
    minval = np.minimum(smallest_p, smallest_c)
    
    pregions = np.full((len(parent_regions), biggest_p, 2), minval - 1)
    cregions = np.full((len(child_regions), biggest_c, 2), minval - 1)
    
    for idx in np.arange(len(parent_regions)):
        pregion = parent_regions[idx]
        pregions[idx,:pregion.shape[0],:] = pregion

    for idx in np.arange(len(child_regions)):
        cregion = child_regions[idx]
        cregions[idx,:cregion.shape[0],:] = cregion

    if pregions.shape[1] == 1 and cregions.shape[1] == 1:
        links = link_parent_and_child_regions(pregions[:,0,:].astype('int32'),
                                              cregions[:,0,:].astype('int32'),
                                              allow_partial = True
                                                  )
    else:
        links = link_parent_and_child_multi_regions(pregions.astype('int32'),
                                                        cregions.astype('int32'),
                                                        cutoff = minval-1,
                                                        allow_partial = True
                                                         )
    
    return links
    