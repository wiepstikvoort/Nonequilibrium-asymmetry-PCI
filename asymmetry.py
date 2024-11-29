# -*- coding: utf-8 -*-
"""
Goal: Calculate the asymmetry of a connectivity matrix

"""

#%% Import necessary packages

import numpy as np

#%% Calculate asymmetry

def calc_symm_pairs(EC, threshold = True, threshold_mean = False, thres_value = .012):
    """
    Parameters
    ----------
    EC : 
        generative Effective Connectivity, in format (n_roi,n_roi)
    threshold : 
        can be set to False if you don't want to use a threshold. Default = True.
    threshold_mean : 
        sets the threshold to the average of EC. Default = False
    thres_value : 
        threshold, float. Default = 0.012

    Returns
    -------
    sum_nr_pairs : 
        the number of pairs that surpass the asymmetry threshold

    """
    
    diff_EC         = EC - EC.T
    if threshold == True: 
        if threshold_mean == True:
            thres = np.mean(abs(diff_EC))
        else:
            thres = thres_value
        
        diff_EC = abs(diff_EC) > thres
        
    nr_of_pairs     = [diff_EC != 0]
    sum_nr_pairs    = np.sum(nr_of_pairs)/2
    
    return int(sum_nr_pairs)