Implementation for manuscript: ``A Hierarchical Bayesian Model for Studying the Impact of Stroke on Brain Motor Function''

FILES:
* func_hie_bvar_fmri_multisub.r
The file contains the key functions to compute the MCMC posterior samples for the hierarchical Bayesian model for multi-subject fMRI data. 

* main_code.cpp
The file contains the main cpp code for updating the subject-level parameters, which is called by func_hie_bvar_fmri_multisub.r

* side_code.cpp
The file contains the required functions for main_code.cpp

* func_prelim_est.r 
The file contains the functions to perform preliminary analysis (non-Bayesian) for fMRI data, including BOLD amplitude and HRF estimation, and BVAR model fitting for replicated time series. 

* test_toy.r
The file contains a toy example of how to use the provided function to estimate the model. 

* hrfbasis_J5 
The file contains the constrained linear basis for HRF

* data/toyData.RData
The file is a R workspace containing a simulated dataset
 

THE KEY FUNCTION:
* hbvar_fmri()
 - Description: 
	the function to perform the hierarchical Bayesian model for multi-subject fMRI data. 
 - Arguments:
	data: 		fMRI data; a list with each element being a matrix of fMRI data for a subject, with each column being the mean fMRI time series (concatenated across sessions) for a ROI. 
	P: 			the number of ROIs (should be equal to the number of columns of the data matrix for each subject)
	T_singleSession: the length of time series for a single session
	R: 			the number of sessions per subject (T_singleSession * R should be equal to the number of rows of the data matrix for each subject)
	sti: 		a binary vector of indicators for whether the stimulus (task) condition is on (1) or off (0) for each TR. Its length should be equal to T_singleSession * R
	TR: 		repetition time of the fMRI scan
	L: 			the lag of VAR model
	prior: 		a list of prior constants. A list containing 6 elements: beta, phi, pi and sigma_beta, sigma_phi, omega. Default is list(phi=list(mean=rep(0, L*P^2*2), sigma2=rep(1, L*P^2*2)), beta = list(mean=rep(0, 3*P), sigma2=rep(1, 3*P)), pi=c(0.1, 0.1), sigma_phi=c(a=1, b=1), sigma_beta=c(a=1,b=1), omega=c(nu=P+1)). For beta and phi, it is a list with the 1st element being the vector of the mean and the 2nd element being the vector of the variance for each element. For pi, it is the argument (a vector of length L+1) for the Dirichlet distribution for pi. For sigma_phi and sigma_beta, it is a vector the two arguments for the gamma distribution for the variances for each element of phi (or beta). For omega, it is the DF argument of the Wishart distribution for Omega. 
	MCMC_setting: specification of the MCMC size; a list containing 3 elements: chain_size (size of the posterior chain), thin (value for thinning, i.e., only every 'thin'-th iterations are recorded in the posterior chain), and burn (burn*thin is the total size of the burn-in period). Default is list(chain_size=1000L, thin=5L, burn=1000L)
	verbose: 	print iteration information during formal estimation. Default is TRUE
	pilot_prior: a list of prior constants used for the pilot estimation. A list containing 4 elements: beta, phi, pi and omega. Default is list(phi=list(mean=rep(0, 2*(P*P)*L), var_cov=diag(rep(0.5^2, 2*(P*P)*L))), beta = list(mean=rep(0, 3*P), var_cov=diag(rep(50^2, 3*P))), pi=matrix(rep(1/(1+L),(1+L)*P*P),1+L), omega=c(nu=P+1, psi_omega = diag(P))). For beta and phi, it is a list with the 1st element being the vector of the mean and the 2nd element being the variance-covariance matrix for the vector phi (or beta). For pi, it is a 2 by (2*L*P^2) matrix, with each column being the probability vector for the multinomial distribution for each element of xi. For omega, it is a list with the 1st element being the DF argument and the 2nd element being the scale matrix argument of the Wishart distribution for Omega. 
	pilot_MCMC_setting: specification of the MCMC size for the pilot estimation, similar to MCMC_setting. Default: list(chain_size = 1000L, thin=5L, burn=1000L)
	pilot_verbose: print iteration information during the pilot estimation
	mc.cores: 	the number of CPUs to use for parallel computing. Default is: min(length(data), 60)
	microtime: 	whether to use a finer scale to compute the convolution. Default is FALSE
	pilot_PLOT: whether to generate a pdf file for results of the preliminary analysis. Default is FALSE
	pilot_pdfname: if pilot_PLOT=TRUE, this will be the file name of the generated pdf file. Default is "pilot_estimation"
	roinames: 	if pilot_PLOT=TRUE, this will be the names of the ROIs. Default is 1:P
	dataname: 	name of the data set to be put in the titles of the plots in the pdf file for the preliminary analysis

 - Values:
	a list of the MCMC posterior samples. Each element is a list containing one set of the posterior MCMC sample for all the subjects. In the list, the first N elements are the N lists of the posterior sample for the subject-level parameters for the N subjects, and the other elements are vphis_multi_mean (group mean for phi), phi_ome (precision of each element of subject-level phi), vbeta_multi_mean (group mean for beta), beta_ome (precision of each element of subject-level beta), vd_multi_mean (empirical group mean for d), d_ome (empirical precision of each element of subject-level d), psi_omega (group-level Omega), prob_ini (group-level probability for each xi). Each list of the subject parameters contains a posterior sample of vbeta (vectorized beta), vphis (vectorized phi), vd (vectorized d) and Omega. 
			
 - Notes (!IMPORTANT!): 
	The order of the phi elements is: phi_{p=1,q=1:P,l=1,k=1}, ..., phi_{p=P,q=1:P,l=1,k=1}, ..., phi_{p=P,q=1:P,l=L,k=1}, ..., phi_{p=P,q=1:P,l=L,k=2}. Here p is the row number and q is the column number. Therefore for phi_array = array(phi, c(P, P, L, 2)), phi_array[,,l,k] returns the *TRANSPOSE* of Phi_{pq}(l) in the manuscript. 
	The order of the beta elements is: beta_{[p=1:P, k=0}, beta_{p=1:P, k=1}, beta_{[p=1:P, k=2}
	The order of the d elements is: d_{[j=1:J, p=1}, ..., d_{j=1:J, p=P}


MORE DETAILS OF THE FILES:
* func_hie_bvar_fmri_multisub.r 
 - hbvar_fmri(): estimate the hierarchical Bayesian model for multi-subject fMRI data; returns the MCMC posterior samples
 - formal_estimation(): an internal function to obtain the MCMC posterior samples based on the given initial values; called by hbvar_fmri()
 - pilot_estimation(): perform pilot estimation to obtain good initial values for formal model estimation, assuming independence across subjects; called by hbvar_fmri()
 - prepare_hrfbasis(): an internal function to return the convolutions and the prior distributions for HRF basis
 - bvarhrf_singlesub_parallel(): an internal function to sample the subject-level parameters for all subjects simultaneous through parallel computing, but independently
 - bvarhrf_singlesub(): an internal wrapper function to sample the subject-level parameters with the function argument being a list of arguments
 - bvarhrf_singlesub_raw(): an internal R wrapper function to sample the subject-level parameters with the function, compiled from C++ code
 - conv_multitrial(): calculate the convolution between the hrf (basis functions) and the stimulus indicator with multiple fMRI sessions (trial)
 - conv_multitrial_fine(): calculate the convolution between the hrf (basis functions) and the stimulus indicator with multiple sessions with a fine grid
 
* func_prelim_est.r 
 - prelim_est(): estimate BOLD amplitude and region-specific HRF assuming independent Gaussian error; called by func_hie_bvar_fmri_multisub.r
 - VAR.multiTrial(): fit a VAR model for replicated time series; called by func_hie_bvar_fmri_multisub.r
 - prelim_Hessian(): an internal function to estimate the Hessian matrix of the point estimates, called by prelim_est()


* hrfbasis_J5
 An internal list containing H (a matrix of the basis vectors), t (time for the basis vectors, in sec), mean (distribution for HRF basis coefficients) and cov (covariance matrix for HRF basis functions) 

* data/toyData.RData
 contains two lists called "data" and "group_parameters". "data" contains 15 simulated data sets (see manuscript for the simulation setting) and the true parameter values. Each dataset is a list containing: y (fMRI data matrix), time (time for the fMRI data, in seconds), sti (a vector of indicator for whether the stimulus condition is on), mu (a matrix of the true mean of the data), Sigma (true variance matrix for the uncorrelated noise), phi1 (phi for the condition 1), phi2 (phi for condition 2), beta (true beta), sdmat (true d), hrf (true hrf), T_singleSession (number of time points for a single session), N_session (number of sessions). "group_parameters" contains true group parameters used to generate the data. 

 

