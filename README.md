* Organization
  - There are four directories, each correspond to a subsection in the experiments section of the paper: lme (linear mixed effects modeling), mixtures (mixture modeling), ml (MovieLens data), and parafac (probabilistic parafac). 
  - Each directory has four sub directories: code, data, qsub, and result.
  - Directory 'code' has files of all the code (Matlab and R source code) that was used in the analysis. 
  - Directory 'data' has (if any) simulated data that was used in the analysis. This directory may be empty or absent.
  - Directory 'qsub' has SGE files (.q) that were used to submit jobs on a SGE cluster. 
  - Directory 'result' has a sub directory 'img' and stores the result (if any) produced in the analysis. This directory may be empty or absent.

* Files
  - 'simulate_data.R' contains the code to simulate and partition the data. 
  - 'analyze_result.R' contains the code for analyzing the results of MCMC, WASP, and competing methods and making plots/tables.
  - 'mcmc_sampler.R' contains the  code for the known/standard MCMC/Gibbs sampler for the model.
  - 'wasp_sampler.R' contains the  code for the MCMC/Gibbs sampler of a subset posterior distribution. This is a modified version of the code in 'mcmc_sampler.R' using stochastic approximation.
  - 'comp_sampler.R' contains the  code for the MCMC/Gibbs sampler of a subset posterior distribution in Consensus Monte Carlo (CMC) or Semiparametric Density Product (SDP). This is a modified version of the code in 'mcmc_sampler.R' by raising the prior to a power of '1/k', where 'k' is the number of subsets.
  - 'variational_bayes.R' contains the  code for the variational Bayes approach.
  - 'submit.R' contains the  code for the R code for submitting a job on the cluster. The files in 'qsub' directory use this file for running simulations.  

* Citation:
  If you use the code, then please cite the following three papers:
  - Srivastava, S., Li, C. and Dunson, D. B. (2017+). Scalable Bayes via barycenter in Wasserstein space. [<https://arxiv.org/abs/1508.05880>]
  - Li, C., Srivastava, S. and Dunson, D. B. (2017). Simple, scalable and accurate posterior interval estimation.  Biometrika 104: 665-680. [<https://arxiv.org/abs/1605.04029>]
  - Srivastava, S., Cevher, V., Tran-Dinh, Q. and Dunson, D. B. (2015). WASP: Scalable Bayes via barycenters of subset posteriors. Artificial Intelligence and Statistics: 912-920.
   
* Contact:
  Please email Cheng Li (<stalic@nus.edu.sg>) or Sanvesh Srivastava (<sanvesh-srivastava@uiowa.edu>) if you have any questions related to the code.

* Acknowledgment
  - The files 'callLpSolver.m', 'recoverSolution.m', and the *.m files for computing the WASP are based on an algorithm due to Volkan Cevher (<https://lions.epfl.ch/>) and Quoc Tran-Dinh (<http://trandinhquoc.com/>). The algorithm can be found in Srivastava et al. (2015).
  - Some code for MovieLens data analysis and linear mixed effects modeling has been borrowed from Patrick O. Perry (<http://ptrckprry.com/code/>).
  - Please email us if you think that we have missed citations to your paper/work. 

