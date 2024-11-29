# -*- coding: utf-8 -*-
"""
For details, please refer to:
Deco G, Sanz Perl Y, Bocaccio H, Tagliazucchi E, Kringelbach ML. 
The INSIDEOUT framework provides 478 precise signatures of the balance of intrinsic and extrinsic dynamics in brain states. 
Commun Biol. 479 2022 Jun 10;5(1):572.

"""

#%% Import necessary packages

import numpy as np

#%% Function insideout

def INSIDEOUT(ts):
    """
    Parameters
    ----------
    ts : 
        timeseries, in format (n_roi, time)

    Returns
    -------
    fowrev : 
        irreversibility, float

    """
    nt          = ts.shape[1]
    A           = ts[:,:nt-1].T
    B           = ts[:,1:].T
    A           = (A - A.mean(axis=0))/A.std(axis=0)
    B           = (B - B.mean(axis=0))/B.std(axis=0)

    FCtf        = np.dot(A.T, B)/A.shape[0]
    FCtr        = np.dot(B.T, A)/B.shape[0]

    Itauf       = -0.5 * np.log(1-np.multiply(FCtf, FCtf))
    Itaur       = -0.5 * np.log(1-np.multiply(FCtr, FCtr))

    reference   = ((Itauf - Itaur))**2

    index       = np.argwhere(reference>np.min(reference))
    fowrev      = np.nanmean(reference[index])
    
    return fowrev


