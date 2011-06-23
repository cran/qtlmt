//intervalMapping.c

//#include "mtcmim.h"

/*-----------------------------------------
 Loglik: calculate loglikelihood
 y: vector of length n
 P: n by m matrix (m: number of components)
 G: m by k matrix (k: number of parameters)
 W: n by l matrix (covariates)
 a: vector of length l
 b: vector of length k
 -----------------------------------------*/
double Loglik(double* y,double** P,double** G,double** W,int n,int m,int k,int l,
  double* a,double* b,double sigma){
  double* mu=new double[m];
  double* y0=new double[n];

  if(m>1) arr_prod(G,b,m,k,mu);
    else if(m==1) mu[0]=0.0;
      else {error(_("m: wrong info...\n"));}
  for(int i=0;i<n;i++){
    y0[i] = y[i];
    for(int j=0;j<l;j++) y0[i] -= W[i][j]*a[j];
  }
  double lik=0.0;
  double tmp;
  for(int i=0;i<n;i++){
    tmp=0.0;
    for(int j=0;j<m;j++){
      tmp += P[i][j]*dnorm(y0[i],mu[j],sigma);
    }
    if(tmp>0) lik += log(tmp);
  }
  delete[] mu; delete[] y0;
  return(lik);
}

/*--------------------------------------------
 update pi (i.e., z), b(i.e., beta) and sigma
 y: vector of length n
 P: n by m matrix
 G: m by k matrix
 W: n by l matrix (covariates)
 a: vector of length l
 b: vector of length k
 z: n by m matrix
 ---------------------------------------------*/
void fz(double* y,double** P,double** G,double** W, int n,int m,int k,int l,
  double* a,double* b,double sigma,double** z){
// update pi (i.e.,z)
  double *mu=new double[m];
  double* y0=new double[n];
  double tmp;

  for(int i=0;i<n;i++){
    y0[i] = y[i];
    for(int j=0;j<l;j++) y0[i] -= W[i][j]*a[j];
  }

  if(m>1) arr_prod(G,b,m,k,mu);
    else if(m==1) mu[0]=0.0;
      else {error(_("m: wrong info...\n"));}
  for(int i=0;i<n;i++){
    tmp=0.0;
    for(int j=0;j<m;j++){
      z[i][j]=P[i][j]*dnorm(y0[i],mu[j],sigma);
      tmp += z[i][j];
    }
    if(tmp>0.0) for(int j=0;j<m;j++){
      z[i][j] /= tmp;
    }
  }
  delete[] mu; delete[] y0;
}

//CM-step in ECM algorithm
void fa(double* y,double** z,double** G,double** W,int n,int m,int k,int l,
  double* b,double* a){
  double *mu=new double[m];
  double** A0=new double*[l]; for(int i=0;i<l;i++) A0[i]=new double[l];
  double* a0=new double[l];
  double tmp;

  if(m>1) arr_prod(G,b,m,k,mu);
    else if(m==1) mu[0]=0.0;
      else {error(_("m: wrong info...\n"));}
  for(int i=0;i<l;i++){
    a0[i]=0.0;
    for(int ir=0;ir<n;ir++){
      tmp=y[ir];
      for(int ic=0;ic<m;ic++){
        tmp -= z[ir][ic]*mu[ic];
      }
      a0[i] += W[ir][i]*tmp;
    }
  }
  for(int i=0;i<l;i++){
    for(int j=0;j<l;j++){
      A0[i][j]=0.0;
      for(int ir=0;ir<n;ir++){
        A0[i][j] += W[ir][i]*W[ir][j];
      }
    }
  }

  ginv(A0,l,l,A0); arr_prod(A0,a0,l,l,a);
//  solve(A0, l, a0, a);

  for(int i=0;i<l;i++) delete[] A0[i];
  delete[] mu; delete[] A0; delete[] a0;
}

//M-step in EM algorithm
void fb(double* y,double** z,double** G,double** W,int n,int m,int k,int l,
  double* a,double* b){
  double** A0=new double*[k]; for(int i=0;i<k;i++) A0[i]=new double[k];
  double* b0=new double[k];
  double* y0=new double[n];

  for(int i=0;i<n;i++){
    y0[i] = y[i];
    for(int j=0;j<l;j++) y0[i] -= W[i][j]*a[j];
  }
  for(int i=0;i<k;i++){
    b0[i]=0.0;
    for(int ir=0;ir<n;ir++){
      for(int ic=0;ic<m;ic++){
        b0[i] += y0[ir]*z[ir][ic]*G[ic][i];
      }
    }
  }
  for(int i=0;i<k;i++){
    for(int j=0;j<k;j++){
      A0[i][j]=0.0;
      for(int ir=0;ir<n;ir++){
        for(int ic=0;ic<m;ic++){
          A0[i][j] += z[ir][ic]*G[ic][i]*G[ic][j];
        }
      }
    }
  }

  ginv(A0,k,k,A0); arr_prod(A0,b0,k,k,b);
//  solve(A0, k, b0, b);

  for(int i=0;i<k;i++) delete[] A0[i];
  delete[] A0; delete[] b0; delete[] y0;
}

void fsigma(double* y,double** z,double** G,double** W,int n,int m,int k,int l,
  double* a,double* b,double& sigma){
  double* mu=new double[m];
  double* y0=new double[n];
  double s=0.0;

  if(m>1) arr_prod(G,b,m,k,mu);
    else if(m==1) mu[0]=0.0;
      else {error(_("m: wrong info...\n"));}
  for(int i=0;i<n;i++){
    y0[i] = y[i];
    for(int j=0;j<l;j++) y0[i] -= W[i][j]*a[j];
  }
  for(int ir=0;ir<n;ir++){
    for(int ic=0;ic<m;ic++){
      s += z[ir][ic]*(y0[ir]-mu[ic])*(y0[ir]-mu[ic]);
    }
  }
  sigma = sqrt(s/n); //cout<<sigma<<endl;
  delete[] mu; delete[] y0;
}

double mimEst(double* y,double** P,double** G,double** W,int n,int m,int k,int l,
  double* a,double* b,double& sigma,int init,int iter,double tol){
  #define INF 1e+308
  double** z=new double*[n]; for(int i=0;i<n;i++) z[i]=new double[m];
  double lik=0.0,lik1,lik2;
  double mx,la,la1;
  if(!init){
    if(k>0) for(int j=0;j<k;j++) b[j] = 0.0;
      else if(m>1){error(_("mimEst: either m or k wrong...\n")); }
  }
  if(sigma<0.0) sigma=sd(y,n);
  lik1 = -INF;
  lik2 = -INF;
  la = INF;

  do{
    fz(y,P,G,W,n,m,k,l,a,b,sigma,z);
    fa(y,z,G,W,n,m,k,l,b,a);
    if(k>0) fb(y,z,G,W,n,m,k,l,a,b);
    fsigma(y,z,G,W,n,m,k,l,a,b,sigma);

    lik2 = lik1;
    lik1 = lik;
    lik = Loglik(y,P,G,W,n,m,k,l,a,b,sigma);
    if(lik==lik1) break;
    la1 = la;
    la = lik1+(lik-lik1)/(1-(lik-lik1)/(lik1-lik2));
    mx = abs(la-la1);
    iter--;
    if(iter<0){
      cout<<"mim: convergence failed..."<<endl;
      break;
    }
  }while(mx>tol);

  for(int i=0;i<n;i++) delete[] z[i];
  delete[] z;
  return lik;
}

