package_install = function(package){if(!package%in%installed.packages()[,1]) install.packages(package)}

package_install("signal")
package_install("Matrix")
package_install("mvtnorm")
package_install("parallel")
package_install("doMC")
package_install("Rcpp")
package_install("RcppArmadillo")
package_install("inline")

library(signal)
library(Matrix)
library(mvtnorm)
library(parallel)
library(doMC)
library(Rcpp)
library(RcppArmadillo)
library(inline)

source("func_prelim_est.r")

main_code = paste(readLines("main_code.cpp"), collapse="\r\n")
side_code = paste(readLines("side_code.cpp"), collapse="\r\n")

bvarhrf_singlesub_raw <- cxxfunction(signature(n_time="int", end="int", stimulus="int", lag="int", signal="numeric", Hconv="numeric", Hsum="numeric",
	prior_d_part="numeric", prior_phi_part="numeric", prior_beta_part="numeric", prior_d_ome="numeric", prior_phi_ome="numeric", 
	prior_beta_ome="numeric", prior_prob="numeric", Nu="integer", psi_omega="numeric",
	phiinits="numeric", vxicatinits="int", omegainits="numeric", vbeta="numeric", vd="numeric",
	noise="numeric", chain_size="int", thin="int", burn="int"), 
main_code, plugin="RcppArmadillo", includes=side_code )


bvarhrf_singlesub = function(paras){
	set.seed(paras$seed)
	return(
		bvarhrf_singlesub_raw(
		paras$T, paras$end, paras$sti, paras$L, paras$y, paras$Hc, paras$Hsum,
		paras$prior_d_part, paras$prior_phi_part, paras$prior_beta_part, paras$prior_d_ome, paras$prior_phi_ome,  
		paras$prior_beta_ome, paras$prob_ini, paras$nu, paras$psi_omega, 
		paras$phiinits, paras$vxicatinits, paras$omegainits, paras$vbeta, paras$vd,
		paras$noise, paras$chain_size, paras$thin, paras$burn)
	)
}

bvarhrf_singlesub_parallel = function(paras_ch, mc.cores = length(paras_ch)){
	return(mclapply(paras_ch, bvarhrf_singlesub, mc.cores=mc.cores))
}


conv_multiSession = function(hrf, sti, end){
#convolution for multiple trials
	T = length(sti)
	x1=rep(NA, T)
	n.ses = length(end) - 1
	for(ses in 1:n.ses){
		sti.ses = sti[(end[ses]+1):end[ses+1]]
		x1[(end[ses]+1):end[ses+1]]=conv(hrf, sti.ses)[1:length(sti.ses)]
	}
	return(x1)
}

conv_multiSession_fine = function(hrffine, sti, end, hrftfine, TR){
#convolution for multiple trials with microtime
	T = length(sti)
	x1=rep(NA, T)
	n.ses = length(end) - 1
	deltat = hrftfine[2] - hrftfine[1]
	nTR = TR/deltat
	for(ses in 1:n.ses){
		sti.ses = sti[(end[ses]+1):end[ses+1]]
		sti.ses.fine = rep(sti.ses, each = nTR)
		x1[(end[ses]+1):end[ses+1]]= deltat*conv(hrffine, sti.ses.fine)[(1:length(sti.ses))*nTR]
	}
	return(x1)
}

prepare_hrfbasis = function(time, end, sti, microtime, TR=time[2]-time[1]){
	load("hrfbasis_J5")
	J = 5
	H = hrfbasisfine$H
	if(!microtime){
		H = H[hrfbasisfine$t %in% c(0, time), ]
		HC1 = apply(H, 2, conv_multiSession, sti, end)
		HC2 = apply(H, 2, conv_multiSession, 1-sti, end)
		delta_t = 1
	}else{
		HC1 = apply(H, 2, conv_multiSession_fine, sti, end, hrfbasisfine$t, TR)
		HC2 = apply(H, 2, conv_multiSession_fine, 1-sti, end, hrfbasisfine$t, TR)
		delta_t = hrfbasisfine$t[2] - hrfbasisfine$t[1]
	}
	Hc = cbind(HC1, HC2)
	Hsum = c(apply(H,2,sum))*delta_t
	scale = c(Hsum%*%hrfbasisfine$mean)
	d_mean = hrfbasisfine$mean/ scale
	d_cov = hrfbasisfine$cov/(scale)^2
	return(list(Hc = Hc, Hsum = Hsum, mean = d_mean, cov = d_cov))
}



	
pilot_estimation = function(data, P, T_singleSession, R, sti, TR, L, prior=list(), MCMC_setting = list(chain_size = 1000L, thin=5L, burn=1000L), verbose, PLOT, pdfname ="pilot_estimation" , roinames=1:P, dataname = "", mc.cores = min(length(data), 60), microtime=FALSE){
# data is a list with each element being the T by P data matrix for a subject
	
	time = (1:(T_singleSession*R))*TR
	T = T_singleSession*R
	end = (0:R)*T_singleSession
	N = length(data) 
	mc.cores = min(mc.cores, N)
	hrfbasis = prepare_hrfbasis(time, end, sti, microtime, TR=TR)
	J = length(hrfbasis$mean)
	Hc = hrfbasis$Hc
	Hsum = hrfbasis$Hsum
	prior_names = names(prior)
		
	prior_d_mu = rep(hrfbasis$mean,P)
	d.i.sig.prior.p = d.i.sig.prior = solve(hrfbasis$cov)
	for(p in 2:P){
		d.i.sig.prior = bdiag(d.i.sig.prior, d.i.sig.prior.p) 
	}
	prior_d_ome = as.matrix(d.i.sig.prior)
	prior_d_part = c(prior_d_ome%*%prior_d_mu)
	
	if(!"phi" %in% prior_names){
		prior_phi_mu = rep(0, 2*(P*P)*L)
		prior_phi_ome = diag(1/rep(0.5^2, 2*(P*P)*L))
	}else{
		prior_phi_mu = prior$phi[[1]]
		prior_phi_ome = solve(prior$phi[[2]])	
	}
	prior_phi_part = c(prior_phi_ome%*%prior_phi_mu)
	
	if(!"beta" %in% prior_names){
		prior_beta_mu = rep(0, 3*P)
		prior_beta_ome = diag(1/rep(50^2, 3*P))
	}else{
		prior_beta_mu = prior$beta[[1]]
		prior_beta_ome = solve(prior$beta[[2]])		
	}
	prior_beta_part = c(prior_beta_ome%*%prior_beta_mu)

	if(!"pi" %in% prior_names){
		prob_ini = matrix(rep(1/(1+L),(1+L)*P*P),1+L)
	}else{
		prob_ini = prior$pi
	}
	
	if(!"omega" %in% prior_names){
		nu = P+1
		psi_omega = diag(P)		
	}else{
		nu = prior$omega[[1]]
		psi_omega = prior$omega[[2]]				
	}
	
	load("hrfbasis_J5")
	if(PLOT) pdf(paste(pdfname,".pdf",sep=""))

	# preliminary analysis for beta and d
	pre_glm = prelim_est(data, T, P, end, Hc, Hsum, J, hrfbasisfine, prior_d_mu, prior_d_ome, 
		roi = roinames, dataname=dataname, verbose = verbose, PLOT=PLOT)
	if(PLOT) dev.off()
	
	paras_ch = vector("list",N)
	seed_ch = sample(N*100, N)
	n_phi = length(prior_phi_mu)
	
	phi.reorder = apply(array(1:n_phi, c(P, P, L)), c(1,3), t)
	for(sub in 1:N){
		y = as.matrix(data[[sub]])#[data$time %in% time,])
		noise  = y - pre_glm$mu_all[,,sub]
		# preliminary analysis for phi and Omega
		phipre = VAR.multiSession(noise, end, L, prior_phi_part[phi.reorder], prior_phi_ome[phi.reorder, phi.reorder])
		phiinits = c(phipre$phi[phi.reorder])
		vxicatinits = matrix(abs(phiinits) > 0.1, P, P)
		omegainits = solve(phipre$Sigma)
		
		paras_ch[[sub]] = list(ID = sub, T = T, end=end, sti=sti, L=L, y=y, Hc=Hc, Hsum=Hsum,
			prior_d_part = prior_d_part, prior_phi_part= prior_phi_part, prior_beta_part = prior_beta_part, 
			prior_d_ome=prior_d_ome, prior_phi_ome=prior_phi_ome, prior_beta_ome=prior_beta_ome, nu=nu, psi_omega=psi_omega,
			prob_ini=prob_ini, 
			phiinits= rep(phiinits, 2), vxicatinits=vxicatinits, omegainits=omegainits, vbeta=c(pre_glm$beta_all[,,sub]), 
			vd=c(pre_glm$d_all[,,sub]), noise=noise, 
			chain_size=MCMC_setting$chain_size, thin=MCMC_setting$thin, burn=MCMC_setting$burn, seed = seed_ch[sub]
		)
	}

	result = bvarhrf_singlesub_parallel(paras_ch, mc.cores = mc.cores)
	return(list(result = result, time = time, T = T, end = end, N = N, mc.cores=mc.cores, hrfbasis=hrfbasis, J = J, Hc=Hc, Hsum = Hsum, hrfbasisfine = hrfbasisfine, MCMC_setting=MCMC_setting, prior_d_mu=prior_d_mu, prior_d_ome=prior_d_ome))

}


formal_estimation = function(pilot, data, P, L, sti, prior=list(), MCMC_setting=list(chain_size = 5000L, burn=3000L), verbose){

	result = pilot$result
	time = pilot$time
	T = pilot$T
	end = pilot$end
	N = pilot$N
	mc.cores=pilot$mc.cores
	hrfbasis=pilot$hrfbasis
	J = pilot$J
	Hc=pilot$Hc
	Hsum = pilot$Hsum
	hrfbasisfine = pilot$hrfbasisfine
	prior_d_mu = pilot$prior_d_mu
	prior_d_ome = pilot$prior_d_ome
	prior_d_part = c(prior_d_ome%*%prior_d_mu)
	
	pre_chain_size = pilot$MCMC_setting$chain_size
	multi_burn = MCMC_setting$burn
	multi_chain_size = MCMC_setting$chain_size
	chain_size = 1L
	thin = 5L

	prior_names = names(prior)	
	if(!"phi" %in% prior_names){
		prior_phi_ome = rep(1, L*P^2*2) 
		prior_phi_mu = rep(0, L*P^2*2)
	}else{
		prior_phi_mu = prior$phi[[1]]
		prior_phi_ome = 1/(prior$phi[[2]])	
	}
	prior_phi_part = c(prior_phi_ome%*%prior_phi_mu)
	
	if(!"beta" %in% prior_names){
		prior_beta_mu = rep(0, 3*P) 
		prior_beta_ome = rep(1/25, 3*P)
	}else{
		prior_beta_mu = prior$beta[[1]]
		prior_beta_ome = 1/(prior$beta[[2]])		
	}
	prior_beta_part = c(prior_beta_ome%*%prior_beta_mu)

	if(!"pi" %in% prior_names){
		a_pi = 0.1
		b_pi = 0.1
	}else{
		a_pi = prior$pi[[1]]
		b_pi = prior$pi[[2]]
	}
	if(!"sigma_phi" %in% prior_names){
		a_phi = 1
		b_phi = 1
	}else{
		a_phi = prior$sigma_phi[[1]]
		b_phi = prior$sigma_phi[[2]]
	}

	if(!"sigma_beta" %in% prior_names){
		a_beta = 1
		b_beta = 1
	}else{
		a_beta = prior$sigma_beta[[1]]
		b_beta = prior$sigma_beta[[2]]
	}
	a_d = b_d = 1
	
	if(!"omega" %in% prior_names){
		nu = P+1
	}else{
		nu = prior$omega[[1]]
	}

	# initialize the variables and preallocate space  
	sapply(result, function(x){x$vphiraw[pre_chain_size,]})
	phi_ome = (1/apply(sapply(result, function(x){x$vphiraw[pre_chain_size,]}),1,sd)^2)
	beta_ome = (1/apply(sapply(result, function(x){x$vbeta[pre_chain_size,]}),1,sd)^2)
	d_ome = diag(prior_d_ome)

	vphiraw_multi = matrix(NA, N, 2*P*P*L)
	vd_multi = matrix(NA, N, P*J)
	vbeta_multi = matrix(NA, N, P*3)
	vxicat_multi = array(NA, c(P,P,N))
	Omega_multi = array(NA, c(P, P, N))
	
	result_less = rep( list(list(vbeta = array(.1, c(chain_size, 3*P)), vphis = array(.1, c(chain_size, P*P*L*2)), vd = array(.1, c(chain_size, J*P)), Omega = array(.1, c(P,P,chain_size)))), N) 
	result_less$vphis_multi_mean = rep(.1, P*P*L*2)
	result_less$phi_ome = rep(1, P*P*L*2)*.1
	result_less$vbeta_multi_mean = rep(.1, 3*P)
	result_less$beta_ome = rep(1, 3*P)*.1
	result_less$vd_multi_mean = rep(.1, J*P)
	result_less$d_ome = rep(1, J*P)*.1
	result_less$psi_omega = diag(P)*.1
	result_less$prob_ini = matrix(.1, 2, P*P*L)				
	multi_result = rep(list(result_less), multi_chain_size)		

	# generate random variables in advance to reduce computation time
	random_B = multi_burn+multi_chain_size
	vphis_multi_mean_random = matrix(rnorm((random_B)*(2*P*P*L)), (random_B))
	vd_multi_mean_random = matrix(rnorm((random_B)*(J*P)), (random_B))
	vbeta_multi_mean_random = matrix(rnorm((random_B)*(3*P)), (random_B))
	phi_ome_random = matrix(rgamma(random_B*2*P*P*L, N/2+a_phi, 1), (random_B))
	d_ome_random = matrix(rgamma(random_B*J*P, N/2+a_d, 1), (random_B))
	beta_ome_random = matrix(rgamma(random_B*3*P, N/2+a_beta, 1), (random_B))
	psi_omega_random = rWishart(random_B, N*nu, diag(P))

	# generate seed number to use in "foreach" at each iteration
	seed.iter = matrix(sample((multi_burn+multi_chain_size)*N, N*(multi_burn+multi_chain_size), 
	replace=((multi_burn+multi_chain_size)*N>2e9)), multi_burn+multi_chain_size)

	options(warn=2)
	registerDoMC(mc.cores)
	getDoParWorkers()	
	
	chain_size = pre_chain_size
	for(iter in 1:(multi_burn + multi_chain_size)){

		if(iter==2){chain_size = 1}
		for(sub in 1:N){
			result_sub = result[[sub]]
			vphiraw_multi[sub, ] = result_sub$vphiraw[chain_size,]
			vd_multi[sub, ] = result_sub$vd[chain_size,]
			vbeta_multi[sub, ] = result_sub$vbeta[chain_size,]
			vxicat_multi[,, sub ] = result_sub$vxicat[,,chain_size]
			Omega_multi[,, sub ] = result_sub$Omega[,,chain_size]		
		}


		############### update for the group

		vphis_multi_mean = apply(vphiraw_multi, 2, mean)
		vphis_multi_sig_post = 1/((phi_ome)*N + prior_phi_ome)
		vphis_multi_mean_post = vphis_multi_sig_post*((phi_ome)*N*vphis_multi_mean + prior_phi_part)
		vphis_multi_mean = c(vphis_multi_mean_post + vphis_multi_mean_random[iter,] * sqrt(vphis_multi_sig_post))
		temp = vphiraw_multi - matrix(rep(vphis_multi_mean, each = N), N)
		phi_ome = (phi_ome_random[iter,] / ((colSums(temp*temp))/2+b_phi))
		phi_part = vphis_multi_mean*phi_ome

		vd_multi_mean_post = apply(vd_multi, 2, mean)
		vd_multi_sig_post = 1/((d_ome)*N)
		vd_multi_mean = c(vd_multi_mean_post + vd_multi_mean_random[iter,] * sqrt(vd_multi_sig_post))
		temp = vd_multi - matrix(rep(vd_multi_mean, each = N), N)
		d_ome = (d_ome_random[iter,] / ((colSums(temp*temp))/2+b_d))

		vbeta_multi_mean = apply(vbeta_multi, 2, mean)
		vbeta_multi_sig_post = 1/((beta_ome)*N + (prior_beta_ome))
		vbeta_multi_mean_post = vbeta_multi_sig_post*((beta_ome)*N*vbeta_multi_mean + prior_beta_part)
		vbeta_multi_mean = c(vbeta_multi_mean_post + vbeta_multi_mean_random[iter,] * sqrt(vbeta_multi_sig_post))
		temp = vbeta_multi - matrix(rep(vbeta_multi_mean, each = N), N)
		beta_ome = (beta_ome_random[iter,] / ((colSums(temp*temp))/2+b_beta))
		beta_part = vbeta_multi_mean*beta_ome

		Omega_sum = apply(Omega_multi, c(1,2), sum )
		psi_omega = Omega_sum%*%solve(psi_omega_random[,,iter])
				
		vxicat_count = apply(vxicat_multi, c(1,2), sum) 
		prob_post = rbeta(P*P, c(vxicat_count)+a_pi, N-c(vxicat_count)+b_pi )
		prob_post = rbind(1-prob_post, prob_post)

		seed = seed.iter[iter, ]
		result = foreach(sub = 1:N)%dopar%{
			set.seed(seed[sub])
			res = result[[sub]]
			return( 
				bvarhrf_singlesub_raw(
				T, end, sti, L, data[[sub]], Hc, Hsum,
				prior_d_part, phi_part, beta_part, prior_d_ome, diag(phi_ome),  
				diag(beta_ome), prob_post, nu, psi_omega, 
				res$vphis[chain_size,], res$vxicat[,,chain_size], res$Omega[,,chain_size], res$vbeta[chain_size,], res$vd[chain_size, ], res$U[,,chain_size],
				1L, thin, 0L)
						
			)
		}
		
		if(iter>multi_burn){
			result_less[1:N] = lapply(result, function(x){x[c("vbeta","vphis","vd","Omega")]})
			result_less$vphis_multi_mean = vphis_multi_mean
			result_less$phi_ome = (phi_ome)
			result_less$vbeta_multi_mean = vbeta_multi_mean
			result_less$beta_ome = (beta_ome)
			result_less$vd_multi_mean = vd_multi_mean
			result_less$d_ome = (d_ome)
			result_less$psi_omega = psi_omega
			result_less$prob_ini = prob_post
						
			multi_result[[iter-multi_burn]] = result_less		
		}

		if(verbose){
			if(iter%%100==0){
				print(iter)
				print(Sys.time())
				flush.console()
			}
		}

	}
	return(multi_result)
}


hbvar_fmri = function(data, P, T_singleSession, R, sti, TR, L , 
	prior= list(), MCMC_setting = list(chain_size = 5000L, burn=3000L), verbose=TRUE, 
	pilot_prior=list(), pilot_MCMC_setting = list(chain_size = 1000L, thin=5L, burn=1000L), pilot_verbose=FALSE, 
	mc.cores = length(data), microtime=FALSE, 
	pilot_PLOT=FALSE, pilot_pdfname ="pilot_estimation", roinames=1:P, dataname = ""){
	mc.cores = min(mc.cores, detectCores())
	pilot = pilot_estimation(data, P, T_singleSession, R, sti, TR, L, pilot_prior, pilot_MCMC_setting, pilot_verbose, pilot_PLOT, pilot_pdfname, roinames, dataname, mc.cores = mc.cores, microtime)
	print("variables initialized")

	return( formal_estimation(pilot, data, P, L, sti, prior, MCMC_setting, verbose) )

}

