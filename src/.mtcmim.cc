//cmtmim.c: composite multiple trait multiple interval mapping

//#include "mtcmim.h"

/*-----------------------------------------
 Loglik: calculate loglikelihood
 y: n by p matrix
 P: n by np matrix, probability for a mixing component
 ma: n by p matrix (means from covariates a)
 mb: np by p matrix (means from G_j)
 invS: inverse of p by p residual covariance matrix sigma
 detS: det(sigma)
 -----------------------------------------*/
double Loglik(double** y,int n,int p,double** P,int np,
	double** ma,double** mb,double** invS,double detS){
	# define PI 3.141592653589793
	double* y0=new double[p];
	double den=pow(2*PI,p/2.0)*sqrt(detS);
	double tt, tmp;

	double lik=0.0;
	for(int in=0;in<n;in++){
		tt=0.0;
		if(np>1)for(int i=0;i<np;i++){
			if(P[in][i]>0){
				tmp=0.0;
				for(int j=0;j<p;j++) y0[j]=y[in][j]-ma[in][j];
				for(int ip=0;ip<p;ip++){
					for(int jp=0;jp<p;jp++){
						tmp += (y0[ip]-mb[i][ip])*invS[ip][jp]*(y0[jp]-mb[i][jp]);
					}
				}
				tmp = exp(-tmp/2.0);
				tt += P[in][i]*tmp/den;
			}
		}
		else if(np==1){
			tmp=0.0;
			for(int j=0;j<p;j++) y0[j]=y[in][j]-ma[in][j];
			for(int ip=0;ip<p;ip++){
				for(int jp=0;jp<p;jp++){
					tmp += y0[ip]*invS[ip][jp]*y0[jp];
				}
			}
			tmp = exp(-tmp/2.0);
			tt += tmp/den;
		}
		else {cout<<"np: wrong..."<<endl; exit(1);}

		if(tt>0) lik += log(tt);
	}

	delete[] y0;
	return(lik);
}

/*-----------------------------------------
 fma: calculate means associated with W and a
 W: n by m matrix: covariates including 1
 nws: nws_j columns of W for y_j, detailed by ws
 a: covariate effects
 ma: n by p matrix (means from a)
 -----------------------------------------*/
void fma(int n,int p,double** W,int* nws,int* ws,double* a,double** ma){
	int* ii=new int[p+1];
	ii[0]=0;
	for(int i=0;i<p;i++) {ii[i+1] = ii[i]+nws[i];}
	for(int i=0;i<n;i++){
		for(int j=0;j<p;j++){
			ma[i][j]=0.0;
			for(int k=0;k<nws[j];k++){
				ma[i][j] += W[i][ws[ii[j]+k]-1]*a[ii[j]+k];
			}
		}
	}
	delete[] ii;
}

/*-----------------------------------------
 fmb: calculate means associated with G and b
 G: np by m genetic matrices
 b: qtl effects
 mb: np by p matrix (means from G_j)
 -----------------------------------------*/
void fmb(int p,double** G,int np,int* ngs,int* gs,double* b,double** mb){
	if(np>1){
		int* ii=new int[p+1];
		ii[0]=0;
		for(int i=0;i<p;i++) {ii[i+1] = ii[i]+ngs[i];}
		for(int i=0;i<np;i++){
			for(int j=0;j<p;j++){
				mb[i][j]=0.0;
				for(int k=0;k<ngs[j];k++){
					mb[i][j] += G[i][gs[ii[j]+k]-1]*b[ii[j]+k];
				}
			}
		}
		delete[] ii;
	}
	else if(np==1)for(int j=0;j<p;j++) mb[0][j]=0.0;
	else {cout<<"np: wrong..."<<endl; exit(1);}
}
/*------------------------------------------------------------------
 update pi (i.e., z), a(covariate effects), b(i.e., beta) and sigma
 ------------------------------------------------------------------*/
void fz(double** y,int n,int p,double** P,int np,double** ma,double** mb,
	double** invS,double** z){
	double* y0=new double[p];
	double tt,tmp;

	if(np>1) for(int i=0;i<n;i++){
		tt=0.0;
		for(int j=0;j<np;j++){
			tmp=0.0;
			if(P[i][j]>0){
				for(int k=0;k<p;k++) y0[k]=y[i][k]-ma[i][k];
				for(int ip=0;ip<p;ip++){
					for(int jp=0;jp<p;jp++){
						tmp += (y0[ip]-mb[j][ip])*invS[ip][jp]*(y0[jp]-mb[j][jp]);
					}
				}
				tmp = exp(-tmp/2.0);
				tmp=P[i][j]*tmp;
			}
			z[i][j]=tmp;
			tt += tmp;
		}
		if(tt>0.0) for(int j=0;j<np;j++){
			z[i][j] /= tt;
		}
	}
	else if(np==1) for(int i=0;i<n;i++) z[i][0]=1.0;
	else {cout<<"np: wrong..."<<endl; exit(1);}
	delete[] y0;
}

//CM-step in ECM algorithm
void fa(double** y,int n,int p,double** W,int* nws,int* ws,
	double** z,int np,double** ma,double** mb,double** invS,double* a){
	int* kk=new int[p+1];
	double* y0=new double[p];
	int kmax=0;
	kk[0]=0;
	for(int i=0;i<p;i++){
		kk[i+1] = kk[i]+nws[i];
		if(nws[i]>kmax) kmax=nws[i];
	}
	double* x=new double[kmax];
	double* x0=new double[kmax];
	double** A=new double*[kmax]; for(int i=0;i<kmax;i++) A[i]=new double[kmax];

	int k;
	double tmp;
	for(int l=0;l<p;l++){
		k=nws[l]; if(k<1) continue;

		for(int ii=0;ii<k;ii++){
			x0[ii]=0.0;
			for(int i=0;i<n;i++){
				for(int t=0;t<p;t++){
					if(t==l) y0[t]=y[i][t];
					else y0[t]=y[i][t]-ma[i][t];
				}
				for(int j=0;j<np;j++){
					tmp=0.0;
					for(int t=0;t<p;t++) tmp += invS[l][t]*(y0[t]-mb[j][t]);
					x0[ii] += z[i][j]*W[i][ws[kk[l]+ii]-1]*tmp;
				}
			}
		}
		for(int i=0;i<k;i++){
			for(int j=0;j<k;j++){
				tmp=0.0; 
				for(int t=0;t<n;t++) tmp += W[t][ws[kk[l]+i]-1]*W[t][ws[kk[l]+j]-1];
				A[i][j]=invS[l][l]*tmp;
			}
		}
		ginv(A,k,k,A); arr_prod(A,x0,k,k,x);
//		cholsl(A,k,x0,x);
//		solve(A,k,x0,x);
		for(int i=0;i<k;i++) a[kk[l]+i]=x[i];
		for(int i=0;i<n;i++){
			tmp=0.0;
			for(int j=0;j<k;j++){
				tmp += W[i][ws[kk[l]+j]-1]*x[j];
			}
			ma[i][l]=tmp;
		}

	}

	for(int i=0;i<kmax;i++) delete[] A[i];
	delete[] x; delete[] x0; delete[] A;
	delete[] kk; delete[] y0;
}

//M-step in EM algorithm
void fb(double** y,int n,int p,double** G,int np,int* ngs,int* gs,
	double** z,double** ma,double** mb,double** invS,double* b){
	if(np<=1){cout<<"b: no exits..."<<endl; exit(1);}

	int* kk=new int[p+1];
	double* y0=new double[p];
	int kmax=0;
	kk[0]=0;
	for(int i=0;i<p;i++){
		kk[i+1] = kk[i]+ngs[i];
		if(ngs[i]>kmax) kmax=ngs[i];
	}
	double* x=new double[kmax];
	double* x0=new double[kmax];
	double** A=new double*[kmax]; for(int i=0;i<kmax;i++) A[i]=new double[kmax];

	int k;
	double tt1,tt2,tmp;
	for(int l=0;l<p;l++){
		k=ngs[l]; if(k<1) continue;

		for(int ii=0;ii<k;ii++){
			x0[ii]=0.0;
			for(int i=0;i<n;i++){
				for(int t=0;t<p;t++) y0[t]=y[i][t]-ma[i][t];
				for(int j=0;j<np;j++){
					tmp=0.0;
					for(int t=0;t<p;t++){
						if(t==l) tmp += invS[l][t]*y0[t];
						else tmp += invS[l][t]*(y0[t]-mb[j][t]);
					}
					x0[ii] += z[i][j]*G[j][gs[kk[l]+ii]-1]*tmp;
				}
			}
		}
		for(int i=0;i<k;i++){
			for(int j=0;j<k;j++){
				tt1=0.0; 
				for(int s=0;s<np;s++){
					tt2=0.0;
					for(int t=0;t<n;t++){
						tt2 += z[t][s];
					}
					tt1 += tt2*G[s][gs[kk[l]+i]-1]*G[s][gs[kk[l]+j]-1];
				}
				A[i][j]=invS[l][l]*tt1;
			}
		}
		ginv(A,k,k,A); arr_prod(A,x0,k,k,x);
//		cholsl(A,k,x0,x);
//		solve(A,k,x0,x);
		for(int i=0;i<k;i++) b[kk[l]+i]=x[i];
		for(int i=0;i<np;i++){
			tt1=0.0;
			for(int j=0;j<k;j++){
				tt1 += G[i][gs[kk[l]+j]-1]*x[j];
			}
			mb[i][l]=tt1;
		}

	}

	for(int i=0;i<kmax;i++) delete[] A[i];
	delete[] x; delete[] x0; delete[] A;
	delete[] kk; delete[] y0;
}

void fS(double** y,int n,int p,double** z,int np,double** ma,double** mb,double** sigma){
	double* y0=new double[p];
	double s;

	for(int is=0;is<p;is++){
		for(int js=0;js<p;js++){
			s=0.0;
			for(int i=0;i<n;i++){
				for(int j=0;j<p;j++) y0[j]=y[i][j]-ma[i][j];
				for(int j=0;j<np;j++){
					s += z[i][j]*(y0[is]-mb[j][is])*(y0[js]-mb[j][js]);
				}
			}
			sigma[is][js]=s/n;
		}
	}
	delete[] y0;
}

/*------------------------------------------------------------
 estimates: a(covariate effects), b(i.e., beta) and sigma
 y: n by p, traits
 P: n by np, mixing proportions
 G: np by nG, genetic matrix
 ngs: vector of length p, ngs_j QTL for y_j (detailed by gs)
 W: n by nW, covariates
 ------------------------------------------------------------*/
double mtcmimEst(double** y,int n,int p,double** P,int np,
	double** G,int* ngs,int* gs,double** W,int* nws,int* ws,
	double* a,double* b,double** sigma,int init,int iter,double tol){
	#define INF 1e+308
	double** z=new double*[n]; for(int i=0;i<n;i++) z[i]=new double[np];
	double** ma=new double*[n]; for(int i=0;i<n;i++) ma[i]=new double[p];
	double** mb=new double*[np]; for(int i=0;i<np;i++) mb[i]=new double[p];
	double** invS=new double*[p]; for(int i=0;i<p;i++) invS[i]=new double[p];
	double detS;

	int na=0, nb=0;
	for(int i=0;i<p;i++) {na += nws[i]; nb += ngs[i];}
	if(!init){
		for(int i=0;i<nb;i++) b[i] = 0.0;
		for(int i=0;i<na;i++) a[i]=0.0;
		for(int i=0;i<p;i++)
			for(int j=0;j<p;j++){
				if(i==j) sigma[i][j]=1.0; else sigma[i][j]=0.0;
			}
	}
	detS=inv_det(sigma,p,invS);
	fma(n,p,W,nws,ws,a,ma);
	fmb(p,G,np,ngs,gs,b,mb);

	double lik=0.0,lik1=-INF,lik2=-INF;
	double mx,la=INF,la1;

	do{
		fz(y,n,p,P,np,ma,mb,invS,z);
		fa(y,n,p,W,nws,ws,z,np,ma,mb,invS,a);
		if(np>1) fb(y,n,p,G,np,ngs,gs,z,ma,mb,invS,b);
		fS(y,n,p,z,np,ma,mb,sigma);
		detS=inv_det(sigma,p,invS);

		lik2 = lik1;
		lik1 = lik;
		lik = Loglik(y,n,p,P,np,ma,mb,invS,detS);

		if(lik==lik1) break;
		la1 = la;
		la = lik1+(lik-lik1)/(1-(lik-lik1)/(lik1-lik2));
		mx = abs(la-la1);
		iter--;
		if(iter<0){
			cout<<"mtcmim: convergence failed..."<<endl;
			break;
		}
	}while(mx>tol);

	for(int i=0;i<n;i++) delete[] z[i];
	for(int i=0;i<n;i++) delete[] ma[i];
	for(int i=0;i<np;i++) delete[] mb[i];
	for(int i=0;i<p;i++) delete[] invS[i];
	delete[] z; delete[] ma; delete[] mb; delete[] invS;

	return lik;
}

/*----------------
 calculate P(Q|M)
 pp: BC-1, RIL-selfing-2, RIL-brother-sister-mating-3
 ----------------*/
double getp(int m1,int m2,double d,int pp=1){
	double p;
	p=haldane_inv(d/100);
	if(pp==2) p=2*p/(1+2*p);
	else if(pp==3) p=4*p/(1+6*p);
	else if(pp!=1){
		cout<<"Only allowed: BC-1, RIL-selfing-2 or RIL-brother-sister-mating-3..."<<endl;
		exit(1);
	}
	if(m1==m2) p=1-p;

	return p;
}
void fP(int** A,int nP,int nQ,int** mdat,int n,int nm,
	double** mpos,int* dists_ch,int* dists_mid,double* dists_d,int* mid,int nmid,
	double**P,int pp){
	double* d=new double[nQ];
	double* dtmp=new double[nQ];
	int* x=new int[nQ+2];
	double d0,d1;
	int nd0,nd,m;

	for(int i=0;i<n;i++){
		for(int j=0;j<nP;j++){
			nd0=0;
			for(int k=0;k<nmid;k++){
				m=mid[k];
				nd=0;
				for(int l=0;l<nQ;l++) if(dists_mid[l]==m){dtmp[nd]=dists_d[l];nd++;}
				sort(dtmp,nd,d);
				d0=mpos[m-1][3];
				if(d[nd-1]>d0+0.001){
					d1=mpos[m][3];
					for(int t=0;t<nd;t++)x[t+1]=A[j][nd0+t];
					x[0]=mdat[i][m-1]; x[nd+1]=mdat[i][m];
					P[i][j] *= getp(x[0],x[1],d[0]-d0,pp);
					P[i][j] *= getp(x[nd],x[nd+1],d1-d[nd-1],pp);
					for(int t=1;t<nd;t++){
						P[i][j] *= getp(x[t],x[t+1],d[t]-d[t-1],pp);
					}
					P[i][j] /= getp(x[0],x[nd+1],d1-d0,pp);
				}else{
					for(int t=0;t<nd;t++)x[t+1]=A[j][nd0+t];
					x[0]=mdat[i][m-1]; x[nd+1]=mdat[i][m-1];
					for(int t=0;t<nd+1;t++){
						P[i][j] *= getp(x[t],x[t+1],0.0,pp);
					}
				}
				nd0 += nd;
			}
		}
	}

	delete[] d; delete[] dtmp; delete[] x;
}
