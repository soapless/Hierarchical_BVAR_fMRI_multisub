#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace arma;

RNGScope scope;

#define H_EPS 0.000001
#define BOUND 100
typedef arma::vec::iterator vec_iter;
typedef arma::mat::iterator mat_iter;
typedef arma::cube::iterator cub_iter;

void block_diag(mat & x, const mat& x1, const mat& x2){
	x.zeros();
	int n1 = x1.n_rows, m1 = x1.n_cols, n2= x2.n_rows, m2 = x2.n_cols; 
	if((x.n_rows != (unsigned) (n1+n2)) | (x.n_cols != (unsigned)(m1+m2)))
		x = zeros(n1+n2, m1+m2);
	x.submat(0,0,n1-1, m1-1) = x1;	
	x.submat(n1,m1, n1+n2-1, m1+m2-1) = x2;		
}

mat block_diag(const mat& x1, const mat& x2){
	int n1 = x1.n_rows, m1 = x1.n_cols, n2= x2.n_rows, m2 = x2.n_cols; 
	mat x = zeros(n1+n2, m1+m2);
	x.submat(0,0,n1-1, m1-1) = x1;	
	x.submat(n1,m1, n1+n2-1, m1+m2-1) = x2;		
	return x; 
}


cube noncont_slices(const cube & X, const uvec & index){
	cube Xsub(X.n_rows, X.n_cols, index.n_elem);
	for(unsigned int i=0; i<index.n_elem; i++){
		Xsub.slice(i) = X.slice(index[i]);
	}
	return Xsub; 
}

cube noncont_slices(const cube & X, const uvec & index, int r_start, int c_start, int r_stop, int c_stop){
	cube Xsub(r_stop-r_start+1, c_stop-c_start+1, index.n_elem);
	for(unsigned int i=0; i<index.n_elem; i++){
		Xsub.slice(i) = X.slice(index[i]).submat(r_start, c_start, r_stop, c_stop);
	}
	return Xsub; 
}


rowvec rmnorm(const rowvec& mean, const mat& var, bool inverse = false ){
	int m = mean.n_elem;
	rowvec rand(m);
	mat Lmat = chol(var);
	rand = randn(1,m);
	if(!inverse){
		rand = rand*Lmat + mean;
	}else{
		rand = trans(solve( trimatu(Lmat), rand.t() )) + mean;
	}
	return rand;
}

vec rchisq(const int v, const int p){
	vec res(p);
	for(int i=0;i<p;i++){
		res[i] = rgamma(1,(v-i)/2,2)[0];
	}
	return res;
}

mat rwish(const int v, const mat& S){
	int p = S.n_rows;
	mat CC = chol(S), Z=zeros(p,p);
	Z.diag() = sqrt(rchisq(v,p));
	if(p>1){
		for(int i=1;i<p;i++)
			Z(span(0,i-1),i) = randn(i);	
	}
	mat cross = Z*CC;
	cross = cross.t()*cross;
	return cross;
}


unsigned int sample(const vec& prob){
	unsigned int N = prob.n_elem;
	vec cdf = cumsum(prob);
	double u = as_scalar(randu(1));
	unsigned int i=0;
	for(i=0;i<N;i++){
		if(u<cdf[i]) break;
	}
	if(i==N) i=N-1; 
	return i;
}

/*
void phi_design(mat& U, const uvec& sti, const int L, umat& vxicat, cube& wZ, const uvec& End){
	int No_Run = End.n_elem-1, P=U.n_cols;
	int Psq=P*P, wt;
	wZ.zeros();
	for(int run=0;run<No_Run;run++){
		for(int t=End[run]+L;t<End[run+1];t++){
			wt = t-L*(run+1);
			for(int l=1;l<=L;l++){
				if(sti[t-l]){
					for(int p=0;p<P;p++){
						uvec cond= (vxicat.col(p) >= l);
						wZ.slice(wt).col(p).subvec((l-1)*Psq+ p*P,(l-1)*Psq+(p+1)*P-1) = (cond)%trans(U.row(t-l));
					}
					wZ.slice(t-L*(run+1)).rows((L+l-1)*Psq,(L+l)*Psq-1).zeros();
				}else{
					wZ.slice(t-L*(run+1)).rows((l-1)*Psq,l*Psq-1).zeros();
					for(int p=0;p<P;p++){
						uvec cond= (vxicat.col(p) >= l);
						wZ.slice(wt).col(p).subvec((L+l-1)*Psq+p*P,(L+l-1)*Psq+(p+1)*P-1)=cond%trans(U.row(t-l));
					}
				}
			}
		}
	}
}
*/

void phi_design(mat& U, const uvec& sti, const int L, umat& vxicat, cube& wZ, const uvec& WTIME){
	int P=U.n_cols, wT = WTIME.n_elem, t;
	int Psq=P*P;
	wZ.zeros();
	for(int wt=0; wt<wT; wt++){
		t = WTIME[wt];
		for(int l=1;l<=L;l++){
			if(sti[t-l]){
				for(int p=0;p<P;p++){
					uvec cond= (vxicat.col(p) >= l);
					wZ.slice(wt).col(p).subvec((l-1)*Psq+ p*P,(l-1)*Psq+(p+1)*P-1) = (cond)%trans(U.row(t-l));
				}
				wZ.slice(wt).rows((L+l-1)*Psq,(L+l)*Psq-1).zeros();
			}else{
				wZ.slice(wt).rows((l-1)*Psq,l*Psq-1).zeros();
				for(int p=0;p<P;p++){
					uvec cond= (vxicat.col(p) >= l);
					wZ.slice(wt).col(p).subvec((L+l-1)*Psq+p*P,(L+l-1)*Psq+(p+1)*P-1)=cond%trans(U.row(t-l));
				}
			}
		}
	}

}


void star(rowvec& vphis, const rowvec& vphi, const int j, const int cat, const int len, const int L){
	for(int l=1; l<=cat;l++){
		vphis[(l-1)*len+j] = vphi[(l-1)*len+j];
		vphis[(l-1+L)*len+j] = vphi[(l-1+L)*len+j];
	}
	for(int l=cat+1;l<=L;l++){
		vphis[(l-1)*len+j] = 0;
		vphis[(l-1+L)*len+j] = 0;		
	}
}
			
rowvec star(const rowvec& vphi, const umat& vxicat){
	rowvec vphis = vphi;
	umat::const_iterator iter = vxicat.begin();
	int len = vxicat.n_elem;
	int L = vphi.n_elem/2/len;
	for(int i=0; i< len; i++){
		for(int l=(*iter);l<L;l++){
			vphis[l*len+i] = 0;
			vphis[(l+L)*len+i] = 0;
		}
		iter++;
	}
	return vphis;
}


void white(const mat& y, const mat& y_center, const uvec& sti, const cube& matphi1, const cube& matphi2, mat& wwy, mat& wwy_center, const uvec& End, const uvec& wEnd){
	int L=matphi1.n_slices, m = End.n_elem-1; 
	mat multiply, wy, wy_center;
	if(L==0){wwy=y; wwy_center = y_center;}
	else{
		for(int run=0;run< m ;run++){
			wy = y.rows(End[run]+L,End[run+1]-1);
			wy_center = y_center.rows(End[run]+L,End[run+1]-1);
			for(unsigned int t=End[run]+L;t<End[run+1];t++){
				for(int l=0;l<L;l++){
					multiply = (sti[t-l-1])?(matphi1.slice(l)):(matphi2.slice(l));
					wy.row(t-End[run]-L) -= y.row(t-l-1)*multiply;					
					wy_center.row(t-End[run]-L) -= y_center.row(t-l-1)*multiply;
				}
			}
			wwy.rows(wEnd[run],wEnd[run+1]-1) = wy;
			wwy_center.rows(wEnd[run],wEnd[run+1]-1) = wy_center;			
		}
	}
}


void white(const mat& x1, const mat& x2, const uvec& sti, const cube& matphi1, const cube& matphi2, cube& wX, const uvec& End, const uvec & wEnd){
	int T = x1.n_rows, L = matphi1.n_slices, P = x1.n_cols, m = End.n_elem-1;
	mat multiply(P,P);
	cube X(3*P,P,T);
	for(int t=0;t<T;t++){
		X.slice(t)=join_cols(join_cols(diagmat(x1.row(t)),diagmat(x2.row(t))),diagmat(ones(1,P)));
	}
	for(int run=0;run<m;run++){
		wX.slices(wEnd[run],wEnd[run+1]-1) = X.slices(End[run]+L,End[run+1]-1);
		if(L>0){
			for(unsigned int t=End[run]+L;t<End[run+1];t++){
				for(int l=0;l<L;l++){
					multiply = sti[t-l-1]?(matphi1.slice(l)):(matphi2.slice(l));
					wX.slice(t-L*(run+1)) -= X.slice(t-l-1)*multiply;
				}
			}
		}
	}
}

void white(const mat& Hc, const rowvec vec_beta, const uvec& sti, const cube& matphi1, const cube& matphi2, cube& wHX, const uvec& End, const uvec & wEnd){
	int T = Hc.n_rows, J = Hc.n_cols/2, L = matphi1.n_slices, P = matphi1.n_rows, m = End.n_elem-1;
	mat multiply(P,P);
	cube X = zeros(J*P,P,T);
	for(int t=0;t<T;t++){
		for(int p=0; p<P;p++){
			X.slice(t).col(p).rows(J*p, J*(p+1)-1) = trans(Hc.submat(t, 0, t, J-1))*vec_beta[p] + trans(Hc.submat(t, J, t, 2*J-1))*vec_beta[P+p];
		}
	}
	for(int run=0;run<m;run++){
		wHX.slices(wEnd[run],wEnd[run+1]-1) = X.slices(End[run]+L,End[run+1]-1);
		if(L>0){
			for(unsigned int t=End[run]+L;t<End[run+1];t++){
				for(int l=0;l<L;l++){
					multiply = sti[t-l-1]?(matphi1.slice(l)):(matphi2.slice(l));
					wHX.slice(t-L*(run+1)) -= X.slice(t-l-1)*multiply;
				}
			}
		}
	}

}


void white(const mat& y, const mat& y_center, const uvec& sti, const cube& matphi1, const cube& matphi2, mat& wy, mat& wy_center, const uvec& WTIME){
	int L=matphi1.n_slices, wT = WTIME.n_elem, t=0; 
	mat multiply;
	if(L==0){wy=y; wy_center = y_center;}
	else{
		wy = y.rows(WTIME);
		wy_center = y_center.rows(WTIME);
		for(int wt=0; wt<wT ; wt++){
			t = WTIME[wt];
			for(int l=0; l<L; l++){
				multiply = (sti[t-l-1])?(matphi1.slice(l)):(matphi2.slice(l));
				wy.row(wt) -= y.row(t-l-1)*multiply;					
				wy_center.row(wt) -= y_center.row(t-l-1)*multiply;
			}
		}
	}
}

void white(const mat& x1, const mat& x2, const uvec& sti, const cube& matphi1, const cube& matphi2, cube& wX, const uvec& WTIME){
	int T = x1.n_rows, L = matphi1.n_slices, P = x1.n_cols, wT = WTIME.n_elem, t=0;
	mat multiply(P,P);
	cube X(3*P,P,T);
	for(int t=0;t<T;t++){
		X.slice(t)=join_cols(join_cols(diagmat(x1.row(t)),diagmat(x2.row(t))),diagmat(ones(1,P)));
	}
	wX = noncont_slices(X, WTIME);
	if(L>0){
		for(int wt=0; wt<wT ; wt++){
			t = WTIME[wt];
			for(int l=0; l<L; l++){
				multiply = sti[t-l-1]?(matphi1.slice(l)):(matphi2.slice(l));
				wX.slice(wt) -= X.slice(t-l-1)*multiply;
			}
		}
	}
}

void white(const mat& Hc, const rowvec vec_beta, const uvec& sti, const cube& matphi1, const cube& matphi2, cube& wHX, const uvec& WTIME){
	int T = Hc.n_rows, J = Hc.n_cols/2, L = matphi1.n_slices, P = matphi1.n_rows, wT = WTIME.n_elem, t=0;
	mat multiply(P,P);
	cube X = zeros(J*P,P,T);
	for(int t=0;t<T;t++){
		for(int p=0; p<P;p++){
			X.slice(t).col(p).rows(J*p, J*(p+1)-1) = trans(Hc.submat(t, 0, t, J-1))*vec_beta[p] + trans(Hc.submat(t, J, t, 2*J-1))*vec_beta[P+p];
		}
	}
	wHX = noncont_slices(X, WTIME);
	if(L>0){
		for(int wt=0; wt < wT; wt++){
			t = WTIME[wt];
			for(int l=0; l<L; l++){
				multiply = sti[t-l-1]?(matphi1.slice(l)):(matphi2.slice(l));
				wHX.slice(wt) -= X.slice(t-l-1)*multiply;
			}
		}
	}

}


void lmMean(const mat& x,const mat& z, const rowvec& vec_beta, mat& mean){
	int P= mean.n_cols;
	mat xt=x, zt=z;
	xt.each_row() %=vec_beta.subvec(0,P-1);
	zt.each_row() %=vec_beta.subvec(P,2*P-1);
	mean= xt + zt;
	mean.each_row() += vec_beta.subvec(P*2,P*3-1);
}

void lmMean(const cube& wX, const rowvec& vec_beta, mat& mean){
	int n = wX.n_slices, P=wX.n_cols;
	mean = zeros(n,P);
	for(int i=0;i<n;i++){
		mean.row(i) = vec_beta* wX.slice(i);
	}
}

void lm2(const mat& wy, const cube& wX, const mat& Omega, rowvec& numer, mat& denom){
	int P=wy.n_cols, n = wy.n_rows, k = wX.n_rows;
	denom.zeros(), numer.zeros();
	mat twX(P,k);
	for(int t=0; t<n; t++){
		twX = trans(wX.slice(t));
		denom=denom+wX.slice(t)*Omega*twX; 
		numer=numer+wy.row(t)*Omega*twX;
	}
}


uvec count(const umat& vxicat, const int L){
	uvec result(L+1);
	for(int j=0; j<=L; j++){
		result[j]= accu(vxicat==j);
	}
	return result;
}

vec Dir(vec alpha){
	int n = alpha.n_elem;
	vec temp(n);
	for(int i=0;i<n;i++){
		temp[i] = ::Rf_rgamma(alpha[i],1.0);
	}
	return temp/sum(temp);
}
