# Noneq-asym-PCI-private

In this repository you can find the code that was used for the results in 'Nonequilibrium dynamics elicited as the origin of perturbative complexity'. 
Some functions are not available with Scipy versions older than v1.13.1

## Data
- gEC_W.npy and gEC_n3.npy
  - the generative effective connectivity matrices for the sleep dataset (n = 18, W = wakefulness, N3 = deep sleep)
  - The Sleep data set is publicly available at https://github.com/yonisanzperl/Perturbation_in_dynamical_models

## Functions
- fun_hopf.py
  - Simulates a timeseries using the Hopf model, including a perturbation by adjusting the bifurcation parameter for the nodes that you want to stimulate.
  - Deco G, Kringelbach ML, Jirsa VK, Ritter P. The dynamics of resting fluctuations in the brain: 492 metastability and its dynamical cortical core. Sci Rep. 2017 Jun 8;7(1):3095.
- insideout.py
  - Calculates the irreversibility of a timeseries.
  - Deco G, Sanz Perl Y, Bocaccio H, Tagliazucchi E, Kringelbach ML. The INSIDEOUT framework provides 478 precise signatures of the balance of intrinsic and extrinsic dynamics in brain states. Commun Biol. 479 2022 Jun 10;5(1):572.
- PCI_st.py
  - Calculates the State Transition Perturbational Complexity Index.
  - Comolatti R, Pigorini A, Casarotto S, Fecchio M, Faria G, Sarasso S, et al. A fast and general method to 516 empirically estimate the complexity of brain responses to transcranial and intracranial stimulations. 517 Brain Stimulat. 2019 Sep;12(5):1280–9. https://doi.org/10.1016/j.brs.2019.05.013
  - https://github.com/renzocom/PCIst
- asymmetry.py:
  - Calculates the number of asymmetric pairs of nodes in the connectivity matrices. 
- hopf_linear.py
  - Calculates the FC and time-lagged FC based on an SC, hopf frequencies, and a noise parameter
  - Ponce-Alvarez A, Deco G. The Hopf whole-brain model and its linear approximation. Sci Rep. 2024 Jan 558 31;14(1):2615.
- calc_EC.py
  - Contains code to calculated time-lagged covariance matrices, and to update your connectivity matrices when you are fitting a Hopf model with empirical fMRI data.
 
