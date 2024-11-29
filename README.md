# Noneq-asym-PCI-private

In this repository you can find the code that was used for the results in 'Nonequilibrium dynamics elicited as the origin of perturbative complexity'. 
Some functions are not available with Scipy versions older than v1.13.1

- fun_hopf.py
  - Simulates a timeseries using the Hopf model, including a perturbation by adjusting the bifurcation parameter for the nodes that you want to stimulate.
- hopf_linear.py
  - Calculates the FC and time-lagged FC based on an SC, hopf frequencies, and a noise parameter
- insideout.py
  - Calculates the irreversibility of a timeseries.
- PCI_st.py
  - Calculates the State Transition Perturbational Complexity Index. 
- asymmetry.py:
  - Calculates the number of asymmetric pairs of nodes in the connectivity matrices. 
- calc_EC.py
  - Contains code to calculated time-lagged covariance matrices, and to update your connectivity matrices when you are fitting a Hopf model with empirical fMRI data.
