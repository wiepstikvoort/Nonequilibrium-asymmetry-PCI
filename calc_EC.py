# -*- coding: utf-8 -*-
"""

Goal: to calculate time-lagged covariance, and to update a connectivity matrix
    when fitting the (linearised) hopf model
    
"""

#%% Import necessary packages

import numpy as np
import scipy.signal as ssl

#%% time lagged covariance without SC
def calc_EC_wo(tss, timelag = 1):
    """
    
    Parameters
    ----------
    tss : 
        non-perturbed timeseries, in format (n_roi, n_timesteps)
    timelag : 
        the number of timesteps of your timelag, default = 1
    
    Returns
    -------
    EC :
        time-lagged cov matrix in format(n_roi, n_roi)

    """

    n_roi = tss.shape[0]
    EC    = np.zeros((n_roi,n_roi))
    for i in range(n_roi):
        for j in range(n_roi):
            correlation = ssl.correlate(tss[i,:] - tss[i,:].mean(), tss[j,:] - tss[j,:].mean(), mode = 'same')
            lags        = ssl.correlation_lags(tss[i,:].shape[0], tss[j,:].shape[0], mode = 'same')
            EC[i,j]     = correlation[lags == timelag] / tss.shape[1]
                
            
    return EC

#%% 

def update_EC(eps_fc, eps_cov, FCemp, FCsim, covemp, covsim, SC, only_pos = True, norm_term = 0.2):
    """
    Parameters
    ----------
    eps_fc : 
        parameter, float
    eps_cov : 
        parameter, float
    FCemp : 
        empirical functional connectivity, format (n_roi, n_roi)
    FCsim : 
        simulated functional connectivity, format (n_roi, n_roi)
    covemp : 
        empirical effective connectivity, format (n_roi, n_roi)
    covsim : 
        simulated effective connectivity, format (n_roi, n_roi)
    SC : 
        structural connectivity, format (n_roi, n_roi)
    only_pos : 
        to keep the update of the SC in positive values. Default = True
            
    Returns
    -------
    An updated SC, format (n_roi, n_roi)

    """
    n_roi = SC.shape[0]
                    
    SCnew = SC + eps_fc * (FCemp - FCsim) + eps_cov * (covemp - covsim)
    for i in range(n_roi):
        for j in range(n_roi):
            if SC[i,j] == 0:
                SCnew[i,j] = 0
                
    if only_pos == True:
        SCnew[SCnew < 0] = 0
    SCnew /= np.max(abs(SCnew))
    SCnew *= norm_term
    
    return SCnew

#%% 

def calc_sigratio(covsim):
    """
    The calc_sigratio function calculates the normalization factor for the 
    time-lagged covariance matrix. This is used so that the FC, which is a 
    covariance normalized by the standard deviations of the two parts, and the 
    tauCOV are in the same space, dimensionless, to calculate the error in the 
    iterative fitting process. 
    
    Parameters
    ----------
    covsim : simulated tss put through calc_EC, format (n_roi,n_roi)

    Returns
    -------
    sigratios in format (n_roi,n_roi)

    """     
    sr = np.zeros((covsim.shape))        
    for i in range(covsim.shape[0]):
        for j in range(covsim.shape[1]):

            sr[i,j] = 1/np.sqrt(abs(covsim[i,i]))/np.sqrt(abs(covsim[j,j]))
    
    return sr
