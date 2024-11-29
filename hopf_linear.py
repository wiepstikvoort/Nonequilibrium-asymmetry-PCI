# -*- coding: utf-8 -*-
"""
For details refer to:
Ponce-Alvarez A, Deco G. 
The Hopf whole-brain model and its linear approximation. 
Sci Rep. 2024 Jan 558 31;14(1):2615.
    
"""

#%% import necessary packages
import numpy as np
from scipy import linalg

#%% Functions

def correlation_from_covariance(covariance):
    """
    Parameters
    ----------
    covariance : matrix with covariance values, format (n_roi, n_roi)

    Returns
    -------
    correlation : matrix with correlation values, format (n_roi, n_roi)

    Computes the correlation matrix from a covariance matrix.
    code from: https://gist.github.com/wiso/ce2a9919ded228838703c1c7c7dad13b
    Replaces MATLAB's corrcov function.
    """
    v             = np.sqrt(np.diag(covariance))
    outer_v       = np.outer(v, v)
    correlation   = covariance / outer_v
    correlation[covariance == 0] = 0
    return correlation

def hopf_int(gC, f_diff, sigma, a = -0.02):
    """
    Parameters
    ----------
    gC : 
        (generative) SC, format (n_roi, n_roi)
    f_diff : 
        Hopf frequencies, format (1, n_roi)
    sigma : 
        noise amplitude, float
    a : 
        bifurcation parameter, default = -0.02

    Returns
    -------
    FC : 
        functional connectivity matrix, format (n_roi, n_roi)
    CV : 
        time-lagged covariance, format (n_roi, n_roi)
    """

    N           = len(gC)
    wo          = np.atleast_2d(f_diff).T.conj() * (2 * np.pi)

    Cvth        = np.zeros((2*N, 2*N))

    # Jacobian
    s           = np.sum(gC, axis=1)
    B           = np.diag(s)

    Axx         = a * np.eye(N) - B + gC
    Ayy         = Axx
    Axy         = -np.diag(wo[:,0])
    Ayx         = np.diag(wo[:,0])

    A           = np.block([[Axx, Axy], [Ayx, Ayy]])
    Qn          = (sigma**2) * np.eye(2*N)

    # Solve Sylvester equation
    Aconjtrans  = np.atleast_2d(A).T.conj()

    Cvth        = linalg.solve_sylvester(A, Aconjtrans, -Qn)
    FCth        = correlation_from_covariance(Cvth)
    FC          = FCth[0:N, 0:N]
    CV          = Cvth[0:N, 0:N]

    return FC, CV, Cvth, A