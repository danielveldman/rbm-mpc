# RBM-MPC

This is the code we used for the numerical experiments with RBM-MPC, a combination of Model Predictive Control (MPC) and Random Batch Methods (RBMs). 
There are three critical parameters that need to be tuned in RBM-MPC, the prediction horizon $T$, the control horizon $\tau$, and the time scale of the randomization $h$. 
The RBM-MPC algorithm is applied to control a heat equation with periodic boundary conditions with is discretized by finite differences with $n$ spatial grid points. 
The code in this repository has been used to generate the numerical examples in the paper TODO. 

<p align="center">
<img src="https://github.com/danielveldman/rbm-mpc/blob/main/Figure.PNG" width="80%" height="80%" >
</p>

The file [C0_basic_code](C0_basic_code) computes one realization of an RBM-MPC trajectory. 
The files [C1_multiple_realizations](C1_multiple_realizations), [C2_varyh](C2_varyh), [C3_varyT](C3_varyT), [C4_varytau](C4_varytau), 
and [C5_varyn](C5_varyn) compute multiple realizations and show the influence of the parameters that define the RBM-MPC strategy. 
