rm(list=ls())#
 
load("data/toyData.RData")

data_sub1 = data[[1]]
P = dim(data_sub1$y)[2]
T_singleSession = data_sub1$T_singleSession
R = data_sub1$N_session
sti = as.integer(data_sub1$sti)
TR = data_sub1$time[2] - data_sub1$time[1]

L = 1L
prior_phi_mu = rep(0, 2*(P*P)*L)
prior_phi_sig = diag(rep(0.5^2, 2*(P*P)*L))
prior_beta_mu = rep(0, 3*P)
prior_beta_sig = diag(rep(50^2, 3*P))
prob_ini = matrix(rep(1/(1+L),(1+L)*P*P),1+L)
nu = P+1
Omega_mu = diag(P)

pilot_prior = list(beta = list(prior_beta_mu, prior_beta_sig), phi=list(prior_phi_mu, prior_phi_sig), pi = prob_ini, omega = list(nu, Omega_mu))

prior = list(phi = list(rep(0, L*P^2*2), rep(1, L*P^2*2)), beta= list(rep(0, 3*P), rep(25, 3*P)), pi = rep(.1, L+1), sigma_phi=c(1,1), sigma_beta=c(1,1), omega = P+1)

data.fmrionly = lapply(data, function(x){x$y})
mc.cores = length(data)

source("func_hie_bvar_fmri_multisub.r")

test = hbvar_fmri(data.fmrionly, P, T_singleSession, R, sti, TR, L, 
	prior= prior, verbose=TRUE, MCMC_setting = list(chain_size = 100L, burn=100L), 
	pilot_prior=pilot_prior, pilot_MCMC_setting = list(chain_size = 10L, thin=5L, burn=10L), pilot_verbose=TRUE, 
	mc.cores = mc.cores, microtime=FALSE, 
	pilot_PLOT=FALSE)


