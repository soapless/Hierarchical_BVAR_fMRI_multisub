arma_rng::set_seed_random();
// ********** data ********* //
const int T = as<int>(n_time);
const uvec sti=as<uvec>(stimulus);
const uvec End = as<uvec>(end); // end is in terms of TR (number of obs) not time in sec.
const mat y=as<mat>(signal);
const mat Hc = as<mat>(Hconv); // convolution of hrf basis and condition indicator, T by 2P
const rowvec Hs = as<rowvec>(Hsum);
const int L=as<int>(lag);


// *********** priors ************* //
const rowvec Prior_d_part = as<rowvec>(prior_d_part);
const mat Prior_d_ome = as<mat>(prior_d_ome);
const rowvec Prior_phi_part = as<rowvec>(prior_phi_part),  Prior_beta_part = as<rowvec>(prior_beta_part);
const mat Prior_phi_ome = as<mat>(prior_phi_ome),  Prior_beta_ome = as<mat>(prior_beta_ome);
mat prob = as<mat>(prior_prob); 
const mat Psi_Omega = as<mat>(psi_omega);
const int nu = as<int>(Nu);

// ****** initial parameter values ******* //
rowvec vphis=as<rowvec>(phiinits); 
umat vxicat=as<umat>(vxicatinits); 
mat Omega=as<mat>(omegainits); 
rowvec vec_beta = as<rowvec>(vbeta);
rowvec vec_d = as<rowvec>(vd);
mat U = as<mat>(noise);
// mat y_center = as<mat>(y_nointer);
mat y_center = y;
y_center.each_row() -= vec_beta.subvec(2*y.n_cols, 3*y.n_cols-1);

// ****** chain size ******* ///
const int B = as<int>(chain_size), Thin= as<int>(thin);
const int burnin= as<int>(burn);
const int newB = B + burnin;

//const bool Return_residual = as<bool>(return_residual); 

// ****** dimension and auxiliary variables ******* ///
const int J = Hc.n_cols/2;
const int P = y.n_cols, JP = J*P;
int n_M = 3*P;
const int Psq = P*P, n_phi = 2*L*P*P;
int No_Run = End.n_elem-1;

const int wT = T - L*No_Run;
uvec WTIME(wT); // the index for the time points with the first L obs are discarded for each session
for(unsigned int run=1;run<End.n_elem;run++){
	WTIME.subvec(End[run-1]-(run-1)*L, End[run]-(run)*L-1)= linspace<uvec>(End[run-1]+L, End[run]-1, End[run]-End[run-1]-L);
}

uvec WTIME_past(wT); 
WTIME_past = WTIME; 
WTIME_past -= L;

uvec sti0_past = find(sti.elem(WTIME_past)==0); 
uvec sti1_past = find(sti.elem(WTIME_past)==1); 


int iter=0, p=0, q=0, jump=0, index=0; 


// ******* create model parameters *********** ///
mat wy = zeros(wT, P);
mat SSE(P,P), resid(wT,P), vx1_new(T,P), vx2_new(T,P), //vx1_old(T,P), vx2_old(T,P), 
 wU(wT,P), wmean(wT,P), wy_center(wT,P); 

mat ibeta_vcov(n_M,n_M), iphi_vcov(n_phi,n_phi), post_phi_ome(n_phi,n_phi), post_phi_sig(n_phi,n_phi), post_beta_ome(n_M,n_M), post_beta_sig(n_M,n_M);
mat id_vcov(JP,JP), post_d_ome(JP,JP), post_d_sig(JP,JP), dmat(J,P);

rowvec SSXy=zeros(1,n_M), post_beta_mu=zeros(1,n_M);
rowvec SSHXy = zeros(1,JP), post_d_mu = zeros(1, JP); 

rowvec vphi=zeros(1,n_phi), SSZU=zeros(1,n_phi), post_phi_mu=zeros(1,n_phi);

mat probpost(1+L,Psq);

cube matphi1=zeros(P,P,L), matphi2=zeros(P,P,L), wX_new(n_M,P,wT), wZ(n_phi,P,wT), wHX(JP, P, wT);

matphi1 = cube(vphis.begin(), P,P,L);
matphi2 = cube(vphis.begin()+ Psq*L, P,P,L);

uvec shuffle(Psq);
shuffle = linspace<uvec>(0, Psq-1, Psq);

rowvec SSZU1 = zeros(1,n_phi/2), SSZU2 = zeros(1,n_phi/2);
mat iphi_vcov1 = zeros(n_phi/2, n_phi/2), iphi_vcov2= zeros(n_phi/2, n_phi/2), post_phi_ome1= zeros(n_phi/2, n_phi/2), post_phi_ome2= zeros(n_phi/2, n_phi/2); 
cube wZ1=zeros(n_phi/2, P, sti1_past.n_elem), wZ2=zeros(n_phi/2, P, sti0_past.n_elem); 


// ************** output chain *************** ///
mat vbeta_ch(B, n_M);
mat vd_ch(B, JP), vphis_ch(B, n_phi), vphiraw_ch(B, n_phi);
cube Omega_ch(P,P,B);
ucube vxicat_ch(P,P,B);
cube U_ch(T,P,B); // , wy_ch(wT, P, B); //, vx1_ch(T,P,B), vx2_ch(T,P,B);

for(iter=0; iter<newB; iter++){
	for(jump=0; jump< Thin; jump++){

		// ********** update phi ************* ///
		if(L>0){		
			wU = U.rows(WTIME);
			phi_design(U,sti,L, vxicat, wZ, WTIME);
			if(L==1){
				wZ1 = noncont_slices(wZ, sti1_past, 0, 0, n_phi/2-1, P-1);
				wZ2 = noncont_slices(wZ, sti0_past, n_phi/2, 0, n_phi-1, P-1);
				lm2(wU.rows(sti1_past), wZ1, Omega, SSZU1, iphi_vcov1);
				lm2(wU.rows(sti0_past), wZ2, Omega, SSZU2, iphi_vcov2);
				post_phi_ome1 = Prior_phi_ome.submat(0,0,n_phi/2-1,n_phi/2-1) + iphi_vcov1;
				post_phi_ome2 = Prior_phi_ome.submat(n_phi/2,n_phi/2,n_phi-1,n_phi-1) + iphi_vcov2;
				block_diag(post_phi_sig, inv_sympd(post_phi_ome1), inv_sympd(post_phi_ome2));			
				post_phi_mu = (Prior_phi_part + join_rows(SSZU1, SSZU2))*post_phi_sig;
			}else{
				lm2(wU, wZ, Omega, SSZU, iphi_vcov);
				post_phi_ome = Prior_phi_ome + iphi_vcov ;
				post_phi_sig = inv_sympd(post_phi_ome);
				post_phi_mu = (Prior_phi_part + SSZU)*post_phi_sig;
			}
			vphi = rmnorm(post_phi_mu, post_phi_sig);
			//		cout<<"phi pass"<<endl;

			vphis = star(vphi, vxicat);
							
			//update probpost (calculate prob for lag = 0,1,...,L)
			vec lik(L+1), d(wT), sum_resid_Omega(wT), sum_resid_Omega_cat(wT); 

			double total;
			int cat=0; 
			std::random_shuffle( shuffle.begin(), shuffle.end());

			// ************ update xi ************* ///		
			for(int jj=0;jj<Psq;jj++){
				int j = shuffle[jj];
				lmMean( wZ, vphis, resid);
				resid = wU - resid;


				p = j/P; q = j - p*P;
				cat = vxicat.at(q,p);

				mat Uq=U.col(q);
				uvec which;

				// l = vxicat.at(q,p)
				sum_resid_Omega_cat = sum( resid*diagmat(Omega.col(p)), 1);
				lik[cat] = 0;

				// l < vxicat.at(q,p)
				sum_resid_Omega = sum_resid_Omega_cat;
				for(int l=cat-1;l>=0;l--){
					which = WTIME-l-1; 
					d = Uq.elem(which)%(conv_to< vec >::from(sti.elem(which))*vphi[Psq*l+j]+(1-conv_to< vec >::from(sti.elem(which)))*vphi[Psq*(L+l)+j]);
					lik[l] = sum(2*d%sum_resid_Omega+ (d%d)*Omega.at(p,p));
					lik[l] += lik[l+1]; 
					if(l>0){
						sum_resid_Omega += d*Omega.at(p,p);
					}
				}

				sum_resid_Omega = sum_resid_Omega_cat;
				for( int l=cat+1;l<=L;l++){
					which = WTIME-l;
					d = - Uq.elem(which)%(conv_to< vec >::from(sti.elem(which))*vphi[Psq*(l-1)+j]+(1-conv_to< vec >::from(sti.elem(which)))*vphi[Psq*(L+l-1)+j]);
								
					lik[l] = sum(2*d%sum_resid_Omega+ d%d*Omega.at(p,p));
					lik[l] += lik[l-1]; 
					if(l<L){
						sum_resid_Omega += d*Omega.at(p,p);
					}
				}

				lik = -lik/2;
				lik = lik - max(lik);
				lik = exp(lik);
				lik = lik % prob.col(j);
				total = sum(lik);
				probpost.col(j) = lik/total;
				vxicat.at(q,p) = sample(probpost.col(j));
				star(vphis, vphi, j, vxicat.at(q,p), Psq, L);

			}

			
			//transform vphis to matphi1 and matphi2
			matphi1 = cube(vphis.begin(), P,P,L);
			matphi2 = cube(vphis.begin()+Psq*L, P,P,L);
			
		}

		// prewhiten y
		white(y, y_center, sti, matphi1, matphi2, wy, wy_center, WTIME);
		//prewhiten Hc
		white(Hc, vec_beta, sti, matphi1, matphi2, wHX, WTIME);
		

// *********** calculate HRF basis coefs ************** /// 		
		lm2(wy_center, wHX, Omega, SSHXy, id_vcov);
		post_d_ome = Prior_d_ome + id_vcov ;
		mat post_d_sig = inv_sympd(post_d_ome);
		
		//			cout<<"d pass"<<endl;

		post_d_mu = (Prior_d_part + SSHXy)*post_d_sig;
		vec_d = rmnorm(post_d_mu, post_d_sig);
		dmat = reshape(vec_d, J, P);
		rowvec scale = Hs*dmat;
		dmat.each_row() /= scale;
		vec_d = reshape(dmat, 1, JP);		

		vx1_new = Hc.cols(0,J-1)*dmat;
		vx2_new = Hc.cols(J,2*J-1)*dmat; 
		
		
// *********** calculate beta ************** /// 		
		//prewhiten X
		white(vx1_new, vx2_new, sti,  matphi1, matphi2, wX_new, WTIME);
		
		lm2(wy, wX_new, Omega, SSXy, ibeta_vcov);		
		post_beta_ome = Prior_beta_ome + ibeta_vcov;
		post_beta_sig = inv_sympd(post_beta_ome);
		
		//			cout<<"beta pass"<<endl;

		post_beta_mu = (Prior_beta_part + SSXy)*post_beta_sig;
		vec_beta = rmnorm(post_beta_mu, post_beta_sig);
		
		y_center = y;
		y_center.each_row() -= vec_beta.subvec(2*y.n_cols, 3*y.n_cols-1);

// *********** calculate fit and noise ************** /// 		
		lmMean(vx1_new, vx2_new, vec_beta, U);
		// calculate U (noise)
		U = y - U;

		
// *********** calculate noise variance ************** /// 		
		// update Omega (sigma)
		lmMean(wX_new, vec_beta, wmean);
		resid = wy - wmean;
		SSE = resid.t()*resid;
		Omega =rwish(wT-1 + nu, inv_sympd(SSE+ Psi_Omega));
		
		//			cout<<"omega pass"<<endl;

	}
	
	if(iter>= burnin){
	
		vbeta_ch.row(index) = vec_beta;
		vxicat_ch.slice(index) = vxicat;
		Omega_ch.slice(index) = Omega;
		//resid_ch.slice(index) = resid;
		vphis_ch.row(index) = vphis;
		vphiraw_ch.row(index) = vphi;
		vd_ch.row(index) = vec_d;
		U_ch.slice(index) = U;
		
		index++;
	}
}
											
	return List::create( 
	Named("vbeta") = vbeta_ch, 
	Named("vphis") = vphis_ch, 
	Named("vphiraw") = vphiraw_ch, 
	Named("vd") = vd_ch, 
	Named("vxicat") = vxicat_ch, 
	Named("Omega") = Omega_ch, 
	Named("U") = U_ch 
	//, Named("y_center") = y_center
	);
