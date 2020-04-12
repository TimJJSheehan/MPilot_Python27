import numpy as np
import copy as cp

# Functions to create congruent arrays from arrays that
# that are not congruent due to different transposition
# and/or one array having only some of the dimensions
# of the other array.
#
# This whole thing is based on dimension names, which
# is no surprise, since my intention to use it with
# arrays read from NetCDF files.
#
# In functions where it matters, the first array
# is the one that is intended to be reshaped.
#
# These functions assume dimensions have associated
# names.

def array_is_congruent_to_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):
    
    return src_dim_nms == tgt_dim_nms and src_arr.shape == tgt_arr.shape
    
# def array_is_congruent_to_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):

def can_expand_to_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):

    rtrn_val = True

    if not set(src_dim_nms).issubset(tgt_dim_nms):
        rtrn_val = False
    else:
        for dim_nm in src_dim_nms:
            if src_arr.shape[src_dim_nms.index(dim_nm)] != tgt_arr.shape[tgt_dim_nms.index(dim_nm)]:
                rtrn_val = False
                break

    return rtrn_val

# def can_expand_to_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):
    
def can_transpose_to_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):

    rtrn_val = True
    
    if set(src_dim_nms) != set(tgt_dim_nms):
        rtrn_val = False
    else:
        for dim_nm in src_dim_nms:
            if src_arr.shape[src_dim_nms.index(dim_nm)] != tgt_arr.shape[tgt_dim_nms.index(dim_nm)]:
                rtrn_val = False
                break

    return rtrn_val
        
# def can_transpose_to_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):

def expand_arr_to_match_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):
    
    src_dim_nms = list(src_dim_nms)
    tgt_dim_nms = list(tgt_dim_nms)

    if not can_expand_to_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):
        print 'ERROR: Target array not expandable to target array.'
        raise Exception('data mismatch')

    # dim names ordered for the expanded array
    rtrn_dim_nms = cp.copy(src_dim_nms)
    rtrnDims = list(src_arr.shape)

    # This expands the the target array into a new array
    # and transposes it.
    tgt_dim_nmsTmp = list(tgt_dim_nms)
    tgt_dim_nmsTmp.reverse()
    for dim_nm in tgt_dim_nmsTmp:
        if dim_nm not in src_dim_nms:
            rtrn_dim_nms.insert(0,dim_nm)
            rtrnDims.insert(0,tgt_arr.shape[tgt_dim_nms.index(dim_nm)])

    # build the transposition index list
    transposeLst = []
    for dim_nm in tgt_dim_nms:
        transposeLst.append(rtrn_dim_nms.index(dim_nm))

    if isinstance(src_arr,np.ma.MaskedArray):
        rtrnArr = np.ma.transpose((np.ma.resize(src_arr,rtrnDims)),transposeLst)
    else:
        rtrnArr = np.transpose((np.resize(src_arr,rtrnDims)),transposeLst)

    return rtrnArr

# def expand_arr_to_match_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):

def transpose_arr_to_match_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms,doCopy=False):

    if not can_transpose_to_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms):
        print 'ERROR: Target array not transposable to target array.'
        raise Exception('data mismatch')

    # build the transposition index list
    transposeLst = []
    for dim_nm in tgt_dim_nms:
        transposeLst.append(rtrn_dim_nms.index(dim_nm))

    if isinstance(rtrnArr,np.ma.MaskedArray):
        rtrnArr = np.ma.transpose((np.ma.resize(rtrnArr,rtrnDims)),transposeLst)
    else:
        rtrnArr = np.transpose((np.resize(rtrnArr,rtrnDims)),transposeLst)

    if doCopy: rtrnArr = cp.deepcopy(rtrnArr)

    return rtrnArr

# def transpose_arr_to_match_tgt(src_arr,src_dim_nms,tgt_arr,tgt_dim_nms,doCopy=False):

