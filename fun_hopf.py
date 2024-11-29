# -*- coding: utf-8 -*-
"""

Goal: Create a function that runs the hopf model with perturbations by changing
    bifurcation parameter a 
    
For details of the Hopf model, please refer to:
Deco G, Kringelbach ML, Jirsa VK, Ritter P. 
The dynamics of resting fluctuations in the brain: 492 metastability and its 
dynamical cortical core. 
Sci Rep. 2017 Jun 8;7(1):3095.

"""

#%% Import necessar packages

import numpy as np
 
#%% Function hopf model including perturbations
def perturb_hopf(hopf_frequencies, SC, G, all_nodes = False, pert_mask = None, 
                 length_pert = 20, pert_ampl = 0.05, **pars):
    """
    Parameters
    ----------
    hopf_frequencies : 
        hopf frequencies, format (1,n_roi)
    SC : 
        structural connectivity (or any connectivity matrix that you want to use),
        gives the strengths of the weights between nodes, format (n_roi,n_roi)
    G : 
        global coupling factor, float
    all_nodes : 
        all nodes are set to be perturbed. Default = False
    pert_mask : 
        mask, nodes that are not perturbed set to False or 0, nodes that are
        perturbed set to True or 1, format (2,n_roi)
    length_pert : 
        length of perturbation, number of timesteps that you want the perturbation
        to last, float. Default = 20
    pert_ampl :  
        perturbation amplitude, the value for a during the perturbation for the
        perturbed node(s). Default = .05
    **pars : 
        Example of pars:
        pars = {'n_runs':10, 'n_roi':90, 'pre_pert':50, 'post_pert':50, 
                'TR':0.1, 'dt':0.025}
        Note: when fitting the model, TR is used to match the TR of the empirical
        fMRI data. Here it is used to capture more timepoints.
        
    Raises
    ------
    Exception
        when all_nodes = True while also giving a matrix to pert_mask,
        the function will not run

    Returns
    -------
    tss: 
        simulated timeseries, format (n_runs, n_roi, time)

    """

    #%% perturb_hopf
    if not pert_mask is None and all_nodes == True:
        raise Exception('You parsed a pert_mask while all_nodes was set to True, change all_nodes to False, or remove pert_mask !')

    n_runs      = pars['n_runs']         
    n_roi       = pars['n_roi']            
    pre_pert    = pars['pre_pert']     
    post_pert   = pars['post_pert'] 
    TR          = pars['TR']                 
    dt          = pars['dt']                 
    
    pre_p       = int(pre_pert/TR)
    length_p    = int(length_pert/TR) 
    post_p      = int(post_pert/TR)
    TIME_POINTS = pre_p + length_p + post_p
    sig         = 6e-4
    dsig        = np.sqrt(dt) * sig      

    tss = np.zeros((n_runs, n_roi, TIME_POINTS))
    
    omega = np.matlib.repmat(2 * np.pi * hopf_frequencies.T, 1, 2)     
    omega[:, 0] *= -1                                               
    weighted_conn = G * SC                                          
    sum_conn = np.matlib.repmat(weighted_conn.sum(1, keepdims=True), 1, 2)          
    for subsim in range(n_runs):
        # Initialize variables
        z = 0.1 * np.ones((n_roi, 2))        
        xs = np.zeros((TIME_POINTS, n_roi))   
        nn = 0                                
    
        # Discard the first 2k (transient)
        for t in np.arange(0, 2000+dt, dt):
            a           = -0.02 * np.ones((n_roi, 2))  
            zz          = z[:, ::-1]                       
            interaction = weighted_conn @ z - \
                sum_conn * z                            
            bifur_freq  = a * z + zz * omega            
            intra_terms = z * (z*z + zz*zz)
            
            # Gaussian noise
            noise       = dsig * np.random.normal(size=(n_roi, 2))
            z           = z + dt * (bifur_freq - intra_terms + interaction) + noise
    

        iter0 = 0
        while nn < TIME_POINTS:
            zz = z[:, ::-1]  
            interaction = weighted_conn @ z - \
                sum_conn * z  
            intra_terms = z * (z*z + zz*zz)

            # Gaussian noise
            noise = dsig * np.random.normal(size=(n_roi, 2))
    
            # Integrative step with perturbation
            if nn >= pre_p and nn < int(pre_p + length_p): 
                if all_nodes:
                    a           = pert_ampl * np.ones((n_roi, 2))
                else:
                    a           = -0.02 * np.ones((n_roi, 2)) 
                    a_pert      = (pert_ampl + 0.02) * pert_mask
                    a           = a + a_pert
                    
            else:    
            # Integrative step without perturbation
                a           = -0.02 * np.ones((n_roi, 2))  
    
            bifur_freq = a * z + zz * omega  
            z = z + dt * (bifur_freq - intra_terms + interaction) + noise
            iter0 += 1
            
            if iter0 >= TR/dt:      
                iter0 = 0
                xs[nn, :] = z[:, 0].T   
                nn += 1
    
        tss[subsim, :, :] = xs.T
    return tss

