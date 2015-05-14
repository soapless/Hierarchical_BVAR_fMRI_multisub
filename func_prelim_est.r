library(mvtnorm)


prelim_Hessian = function(y, d, Ups, HC, Hsum, sig, J, prior_d_ome){
# calculate Hessian matrix of -log(likelihood) 
	co0 = cbind(HC[,1:J]%*%d, HC[,1:J+J]%*%d)
	tmp1 = HC[,1:(J-1)] - HC[,J]%*%t(Hsum[1:(J-1)])/Hsum[J]
	tmp2 = HC[,J+1:(J-1)] - HC[,J+J]%*%t(Hsum[1:(J-1)])/Hsum[J]
	co = cbind((HC[,J]/Hsum[J] + (tmp1)%*%d[1:(J-1)]), 
		HC[,J+J]/Hsum[J] + (tmp2)%*%d[1:(J-1)])
	co2 = Ups[2]*(tmp1) + Ups[3]*(tmp2)

	res = c(y - cbind(1,co0)%*%Ups)
	Hes = matrix(NA, 3+J, 3+J)
	Hes[1,1] = length(y)/sig
	Hes[2:3,2:3] = t(co)%*%co/sig
	Hes[3+1:(J-1), 3+1:(J-1)] = t(co2)%*%co2/sig + prior_d_ome[1:(J-1), 1:(J-1)]
	Hes[1,2:3] = Hes[2:3, 1] = colSums(co)/sig; 
	Hes[1, 3+1:(J-1)] = Hes[3+1:(J-1), 1] = colSums(co2)/sig;  
	Hes[2,3+1:(J-1)] = Hes[3+1:(J-1), 2] = (c(-t(res)%*%tmp1) + c(t(co2)%*%co[,1]))/sig
	Hes[3,3+1:(J-1)] = Hes[3+1:(J-1), 3] = (c(-t(res)%*%tmp2) + c(t(co2)%*%co[,2]))/sig

	Hes[3+(J-1)+1,1] = 	Hes[1, 3+(J-1)+1] = sum(res)/sig^2
	Hes[3+(J-1)+1,2:3] = Hes[2:3, 3+(J-1)+1] = colSums(diag(res)%*%co)/sig^2
	Hes[3+(J-1)+1,3+(1:(J-1))] = Hes[3+(1:(J-1)), 3+(J-1)+1] = colSums(diag(res)%*%co2)/sig^2
	Hes[3+(J-1)+1, 3+(J-1)+1] = -T/2/sig^2 + 1/sig^3*sum(res**2)
	varnames = c(paste("beta",1:3,sep=""), paste("d",1:(J-1),sep=""), "sig2")
	dimnames(Hes) = list(varnames , varnames)
	return(Hes)
}

prelim_est = function(data, T, P, end, HC, Hsum, J, hrfbasisfine, prior_d_mu, prior_d_ome, 
	roi = 1:P, dataname="", verbose = TRUE, PLOT=FALSE){
# function to get rough estimates for beta, d and sig, assuming Y = beta0 + X1* beta1 + X2* beta2 + eps
# where Xi = HC_i %*% d
# assuming independence temporally and spatially
# informative prior, but sum(Hd)=1 is satisfied in the estimation

	if(class(data)=="list"){
		N = length(data)
	}else{
		N = 1
		data = list(data)
	}
	if(dim(data[[1]])[1]==P){
		stop(paste(dataname, ": data needs to be transposed to a T by P matrix"))
	}


	d.prior.part = c(prior_d_ome%*%prior_d_mu)
	mutmp = NULL
	isigtmp = array(NA, c(J-1,J-1,P))
	mutmp = matrix(prior_d_mu, J, P)[1:(J-1),]
	for(p in 1:P){
		isigtmp[,,p] = prior_d_ome[(1:(J-1)) + (p-1)*J, (1:(J-1)) + (p-1)*J]
	}

	HC_new = HC
	for(j in 1:(J-1)){
		HC_new[,j] = HC_new[,j] - Hsum[j]*HC_new[,J]/Hsum[J]
		HC_new[,j+J] = HC_new[,j+J] - Hsum[j]*HC_new[,J+J]/Hsum[J] 
	}
	HC_new[,c(J,2*J)] = HC_new[,c(J,2*J)]/Hsum[J]
			
	TOL = 1e-5	# tolerance of error to claim convergence	
	B = 2000L # maximum number of iterations
	M = 10L # number of estimations to search global minimum
		
	d_all = array(0, c(J,P,N))
	beta_all = array(0, c(P, 3, N))
	Hessian_eigen = array(0, c(3+J,P, N))
	sse = lik = array(NA, c(P,M))
	lik_all = array(NA, c(dim(sse), N))
	sig_all = array(NA, c(P, N))
	vd_temp = array(0, c(J,P,M))
	vbeta_temp = array(0, c(P,3,M))
	vbeta_vcov_all = array(0, c(3, 3, P, M))
	vd_vcov_all = array(0, c(J, J, P, M))
	vbeta_vcov_sub = array(0, c(3*P, 3*P, N))
	vd_vcov_sub = array(0, c(J*P, J*P, N))
	betavcov = dvcov = vector("list", P)
	mu_M = array(0, c(T, P, M))
	mu_all = array(0, c(T, P, N))

	for(sub in 1:N){
		y = data[[sub]]
		for(m in 1:M){
			converge = FALSE
			vbeta_ch = array(NA, c(B, 3*P))
			vd_ch = array(NA, c(B, J*P))
			X = Z = ytemp = list()
			vbeta = rep(NA, 3*P)
			vd = rep(NA, J*P)

			######## get initial values ########
			for(p in 1:P){
				covtmp = solve(prior_d_ome[(1:(J)) + (p-1)*J, (1:(J)) + (p-1)*J])
				meantmp = prior_d_mu[(1:J) + (p-1)*J] 
				dtemp = c(rmvnorm(1, meantmp, covtmp))
				scale = as.numeric(matrix(Hsum,1,J)%*%matrix(dtemp,J,1))
				dtemp = dtemp/scale
				X[[p]] = cbind(1, HC[,1:J]%*%dtemp, HC[,(1:J)+J]%*%dtemp)
				m1 = lm(y[,p]~ 0+X[[p]])
				btemp = coef(m1)
				sse[p,m] = sum(m1$res**2)
				vbeta[c(2*P+p, p, P+p)] = btemp
				Z[[p]] = HC[,1:J]*btemp[2] + HC[,(1:J)+J]*btemp[3]
				ytemp[[p]] = y[,p] - btemp[1]
				Z[[p]] = HC_new[,1:(J-1)]*btemp[2] + HC_new[,(1:(J-1))+J]*btemp[3]
				ytemp[[p]] = ytemp[[p]] - HC_new[,J]*btemp[2] - HC_new[,2*J]*btemp[3]
			}

			
			######## begin iterative estimation ########
			for(iter in 1:B){
				vbeta = rep(NA, 3*P)
				vd = rep(NA, J*P)
				
				for(p in 1:P){
					sig = sse[p,m]/T
					# estimate d
					dvcov[[p]] = solve(t(Z[[p]])%*%Z[[p]]/sig + prior_d_ome[(1:(J-1)) + (p-1)*J, (1:(J-1)) + (p-1)*J])
					dtemp = dvcov[[p]]%*%(t(Z[[p]])%*%ytemp[[p]]/sig + d.prior.part[(1:(J-1)) + (p-1)*J])
					dtemp = c(dtemp, (1 - Hsum[1:(J-1)]%*%dtemp)/Hsum[J])
					vd[(p-1)*J+1:J] = dtemp                                                                
					X[[p]] = cbind(1, HC[,1:J]%*%dtemp, HC[,(1:J)+J]%*%dtemp)

					# estimate beta
					m1 = lm(y[,p]~ 0+X[[p]])
					btemp = coef(m1)
					vbeta[c(2*P+p, p, P+p)] = btemp
					Z[[p]] = HC_new[,1:(J-1)]*btemp[2] + HC_new[,(1:(J-1))+J]*btemp[3]
					ytemp[[p]] = y[,p] - btemp[1] - HC_new[,J]*btemp[2] - HC_new[,2*J]*btemp[3]

					# calculate other related quantities
					msum = summary(m1)
					betavcov[[p]] = msum$cov[c(2,3,1), c(2,3,1)] * (msum$sigma**2)
					mu_M[,p,m] = m1$fit
					sse[p,m] = sum(m1$res**2)
					sig = sse[p,m]/T
					lik[p,m] = -T/2*log(sig)- .5*sse[p,m]/sig - .5*t(dtemp[-J]- mutmp[,p])%*%isigtmp[,,p]%*%(dtemp[-J]- mutmp[,p])						
				}

				vbeta_ch[iter, ] = vbeta
				vd_ch[iter, ] = vd

				if(iter>=2){
					if(max(abs(vbeta_ch[iter,]-vbeta_ch[iter-1, ]))<TOL & max(abs(vd_ch[iter,]-vd_ch[iter-1, ]))<TOL){
						converge = TRUE
						if(verbose)
						print( paste(dataname, 'subject', sub, ': converge, iter:', iter))			
						for(p in 1:P){
							vd_vcov_all[1:(J-1), 1:(J-1), p, m] = dvcov[[p]]
							vd_vcov_all[J, J, p, m] = t(Hsum[1:(J-1)]/Hsum[J])%*%dvcov[[p]]%*%(Hsum[1:(J-1)]/Hsum[J])
							vbeta_vcov_all[,,p, m] = betavcov[[p]]
						}
						break
					}
				}
				
				flush.console()
			}
			
			if(!converge){
				if(verbose)
				print(paste(dataname, 'subject', sub, ': not converge, iter:', iter))
			}
			vd_temp[,,m] = matrix(vd,J)
			vbeta_temp[,,m] = matrix(vbeta, P)

		}

		if(verbose){
			print("coefficient of variation of loglik across 10 runs:")
			print(apply(lik, 1, sd)/apply(lik, 1, mean))
		}
		select = apply(-lik, 1, order)[1, ]

		for(p in 1:P){
			vd_vcov_sub[((p-1)*J+1):(p*J), ((p-1)*J+1):(p*J),sub] = 
				vd_vcov_all[,,p,select[p]]
			vbeta_vcov_sub[c(p,p+P, p+2*P), c(p,p+P, p+2*P),sub] = 
				vbeta_vcov_all[,,p,select[p]]
			d_all[,p,sub] = vd_temp[,p,select[p]]
			beta_all[p,,sub] = vbeta_temp[p,,select[p]]
			sig = sig_all[p,sub] = sse[p,select[p]]/T
			mu_all[,p,sub] = mu_M[,p,select[p]]
			Hes = prelim_Hessian(y[,p], d_all[,p,sub], beta_all[p,,sub], HC, Hsum, sig, J, prior_d_ome[(1:J) + (p-1)*J, (1:J) + (p-1)*J])			
			eigenval = eigen(Hes)$val
			Hessian_eigen[,p,sub] = eigenval
		}
		
		lik_all[,,sub] = lik
	}


	####### begin plotting ######
	if(PLOT){
		Hfine = hrfbasisfine$H
		chrf = dgamma(hrfbasisfine$t, 6, 1) - 1/6*dgamma(hrfbasisfine$t, 16, 1)
		chrf = chrf/max(chrf)
		####### y vs mean fit plot and residual vs mean fit plot ######
		for(sub in 1:N){
			y = data[[sub]]
			par(mfrow=c(3,2))
			for(p in 1:P){
				plot(y[,p], type="l", col="black", main=paste("mean fit: sub",sub, " roi", roi[p]), ylab="", xlab="time (sec)")
				lines(mu_all[,p,sub], col="blue")
				# if("end" %in% ls())
				abline(v = end, lwd =1.5, col="green", lty = 2)
				plot(mu_all[,p,sub], (y[,p]-mu_all[,p,sub]), main=paste("sub",sub, " roi", roi[p]), ylab="residual", xlab="mean fit")	
			}
			
			#######plot all hrfs together######
			par(mfrow=c(1,1))
			hsub = Hfine%*% d_all[,,sub]
			hsub = hsub%*%diag(1/apply(hsub, 2, max))
			plot(hrfbasisfine$t, chrf, type="l", col="black", main=paste("HRF - sub",sub), ylab="",xlab="time (sec)", 
			ylim=range(c(chrf, hsub)), #ylim = , range(hsub),
			xlim=c(0, min(29,max(hrfbasisfine$t))))
			for(p in 1:P){
				lines(hrfbasisfine$t, hsub[,p], col=rainbow(P)[p])
			}
		}

		if(N>1){
			########plot hrfs of the same roi, from different subjects ##########
			par(mfrow= c(2,2))
			for(p in 1:P){
				hsub = hrfbasisfine$H%*%d_all[,p,]
				hsub = array(hsub, c(dim(hrfbasisfine$H)[1], N))
				hsub = hsub%*%diag(c(1/apply(hsub, 2, max)))
				plot(hrfbasisfine$t, chrf, type="l", col="black", main=paste("HRF -roi",roi[p]), ylab="",xlab="time (sec)", 
				ylim=c(-1, 1), #ylim = range(hsub)
				xlim=c(0, min(29,max(hrfbasisfine$t))), cex.lab=1.5, cex.axis=1.5)
				for(sub in 1:N){		
					lines(hrfbasisfine$t, hsub[,sub], col=rainbow(N)[sub], lwd=2)
					colval = col2rgb(rainbow(floor(N*1.2))[sub])/255
					col = rgb(colval[1],colval[2],colval[3],alpha=0.5) 
				}
			}

			########plot hrfs of the same subject, from different rois ##########
			par(mfrow= c(2,2))
			for(sub in 1:N){
				hsub = hrfbasisfine$H%*%d_all[,,sub]
				hsub = hsub%*%diag(1/apply(hsub, 2, max))
				plot(hrfbasisfine$t, chrf, type="l", col="black", main=paste("HRF - sub",sub), ylab="",xlab="time (sec)",
				ylim=c(-1,1), #ylim = range(hsub), 
				xlim=c(0,min(29,max(hrfbasisfine$t))), cex.lab=1.5, cex.axis=1.5)
				for(p in 1:P){
					lines(hrfbasisfine$t, hsub[,p], col=rainbow(P)[p], lwd=2)
					colval = col2rgb(rainbow(floor(P*1.2))[p])/255
					col = rgb(colval[1],colval[2],colval[3],alpha=0.5) 
				}	
			}
		}
	}

	return(list(beta_all=beta_all, d_all=d_all, sig_all = sig_all, 
	vd_vcov=vd_vcov_sub, vbeta_vcov=vbeta_vcov_sub, Hessian_eigen = Hessian_eigen, mu_all = mu_all))
}


VAR.multiSession = function(Y, end, L, prior_part, prior_isig){
	n.ses = length(end)-1
	P = dim(Y)[2]
	Sigma = diag(P)
	EPS = 1e-6
	loglik = -Inf
	T = dim(Y)[1]

	Z = ZZ = vy = Zr = list()
	for(i in 1:n.ses){
		T.ses = end[i+1] - end[i]
		t.ses = (end[i]+1):end[i+1]
		vy[[i]] = c(t(Y[t.ses[(L+1):T.ses], ]))
		Ztemp = NULL
		for(l in 1:L){
			Zl = t(Y[t.ses[(L-l+1):(T.ses-l)],])
			Ztemp = rbind(Ztemp, Zl)
		}
		Z[[i]] = Ztemp
		ZZ[[i]] = Z[[i]]%*%t(Z[[i]])
	}

	ZZ.all= 0
	for(i in 1:n.ses){
		ZZ.all = ZZ.all + ZZ[[i]]
	}
	i.ZZ.all = solve(ZZ.all)
	vb = 0
	for(i in 1:n.ses){
		Zr[[i]] = kronecker(i.ZZ.all%*%Z[[i]], diag(P))
		vb = vb + Zr[[i]]%*%vy[[i]]
	}
	for(iter in 1:500){

		phicov = kronecker(i.ZZ.all, Sigma)
		iphicov = solve(phicov)
		cov_post = solve(iphicov + prior_isig)
		vb_post = cov_post%*%(iphicov%*%vb + prior_part)
		
		B = matrix(vb_post, P)
		U = list()
		UU.all = 0
		for(i in 1:n.ses){
			U[[i]] = matrix(vy[[i]],P) - B%*%Z[[i]]		
			UU.all = UU.all + U[[i]]%*%t(U[[i]])
		}
		Sigma = UU.all/(T-L*(P+n.ses)) 
		
		loglik_new = -0.5*((det(Sigma))*(T-L*n.ses)  +  sum(diag(UU.all%*%Sigma)) )

		if(abs(loglik_new - loglik)< EPS) break

		loglik = loglik_new
	}

	return(list(phi=vb_post, phicov = cov_post, Sigma = Sigma, eps = U))
}