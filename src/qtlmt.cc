
//mylib.cc

#include "mtcmim.h"
#include ".basic.cc"
#include ".sure.cc"
#include ".mappingBasic.cc"
#include ".intervalMapping.cc"
#include ".mtcmim.cc"

/********************************************************
 my general R functions
 ********************************************************/
extern "C"{
 	void runifc(double* x,int& n,long *seed){
		runif(x,n,seed);
	}
	void svdc(double *a, int& m, int& n, double* w, double *v){
		double** a0=new double*[m]; for(int i=0;i<m;i++) a0[i]=a+i*n;
		double** v0=new double*[n]; for(int i=0;i<n;i++) v0[i]=v+i*n;
		svd(a0, m, n, w, v0);
		delete[] a0; delete[] v0;
	}
       void cholc(double* A, int& n, double* p){
		double** a=new double*[n]; for(int i=0;i<n;i++) a[i]=A+i*n;
		chol(a,n,p);
		delete[] a;
       }
	void cholsolve(double* A,int& n,double* b,double* x){
		double** a=new double*[n]; for(int i=0;i<n;i++) a[i]=A+i*n;
		cholsl(a,n,b,x);
		delete[] a;
	}
	void lusolve(double* A,int& n,double* b,double* x){
		double** a=new double*[n]; for(int i=0;i<n;i++) a[i]=A+i*n;
		solve(a,n,b,x);
		delete[] a;
	}
}

/********************************************************
 my R functions for QTL mapping
 ********************************************************/
extern "C"{
//SURE models
	void sureEstc(double* y,int& n,int& p,double* x,int& m,int* nqs,int* qs,
		double* b,double* sigma,double& loglik,int& ini_sigma,int& iter,double& tol){
		double** y0=new double*[n];for(int i=0;i<n;i++) y0[i]=y+i*p;
		double** x0=new double*[n];for(int i=0;i<n;i++) x0[i]=x+i*m;
		double** sigma0=new double*[p];for(int i=0;i<p;i++) sigma0[i]=sigma+i*p;
		loglik = sureEst(y0,n,p,x0,m,nqs,qs,b,sigma0,ini_sigma,iter,tol);
		delete[] y0; delete[] x0; delete[] sigma0;
	}
/*	void sureAdd1c(double* y,int& n,int& p,double* x, int& m,int* nqs,int* qs,
		int* nupper,int* upper,double* sigma,int* which, bool& add,
		int& iter, double& tol){
		double** y0=new double*[n];for(int i=0;i<n;i++) y0[i]=y+i*p;
		double** x0=new double*[n];for(int i=0;i<n;i++) x0[i]=x+i*m;
		double** sigma0=new double*[p];for(int i=0;i<p;i++) sigma0[i]=sigma+i*p;
		add = sureAdd1(y0,n,p,x0,m,nqs,qs,nupper,upper,sigma0,which,iter,tol);
		delete[] y0; delete[] x0; delete[] sigma0;
	}
*/
	void sureStepc(double* y,int& n,int& p,double* x,int& m,
		int* nlower,int* lower,int* nupper,int* upper,double& k,int& direction,
		int* vin,double* rec,int& max_terms,int& steps,int& iter, double& tol){
		double** y0=new double*[n];for(int i=0;i<n;i++) y0[i]=y+i*p;
		double** x0=new double*[n];for(int i=0;i<n;i++) x0[i]=x+i*m;
		int** vin0=new int*[p];for(int i=0;i<p;i++) vin0[i]=vin+i*m;
		sureStep(y0,n,p,x0,m,nlower,lower,nupper,upper,k,direction,vin0,rec,
			max_terms,steps,iter,tol);
		delete[] y0; delete[] x0; delete[] vin0;
	}
//single-trait composite multiple-interval mapping
	void mimEstc(double* y,double* P,double* G,double* W,int& n,int& m,int& k,int& l,
		double*a,double* b,double& sigma,double& loglik,
		int& init,int& iter,double& tol){
		double** P0=new double*[n]; for(int i=0;i<n;i++) P0[i]=P+i*m;
		double** G0=new double*[m]; for(int i=0;i<m;i++) G0[i]=G+i*k;
		double** W0=new double*[n]; for(int i=0;i<n;i++) W0[i]=W+i*l;
		loglik=mimEst(y,P0,G0,W0,n,m,k,l,a,b,sigma,init,iter,tol);
		    delete[] P0; delete[] G0; delete[] W0;
	}
//multiple-trait composite multiple-interval mapping
	void mtcmimEstc(double* y,int& n,int& p,double* P,int& np,
		double* G,int& nG,int* ngs,int* gs,double* W,int& nW,int* nws,int* ws,
		double* a,double* b,double* sigma,double& loglik,
		int& init,int& iter,double& tol){
		double** y0=new double*[n]; for(int i=0;i<n;i++) y0[i]=y+i*p;
		double** P0=new double*[n]; for(int i=0;i<n;i++) P0[i]=P+i*np;
		double** G0=new double*[np]; for(int i=0;i<np;i++) G0[i]=G+i*nG;
		double** W0=new double*[n]; for(int i=0;i<n;i++) W0[i]=W+i*nW;
		double** sigma0=new double*[p]; for(int i=0;i<p;i++) sigma0[i]=sigma+i*p;
		loglik=mtcmimEst(y0,n,p,P0,np,G0,ngs,gs,W0,nws,ws,a,b,sigma0,
			init,iter,tol);
		delete[] y0; delete[] P0; delete[] G0; delete[] W0; delete[] sigma0;
	}
	void fPc(int* A,int& nP,int& nQ,int* mdat,int& n,int& nm,
		double* mpos,int* dists_ch,int* dists_mid,double* dists_d,int* mid,int& nmid,
		double* P,int& pp){
		int** A0=new int*[nP];for(int i=0;i<nP;i++) A0[i]=A+i*nQ;
		int** mdat0=new int*[n];for(int i=0;i<n;i++) mdat0[i]=mdat+i*nm;
		double** mpos0=new double*[nm];for(int i=0;i<nm;i++) mpos0[i]=mpos+i*4;
		double** P0=new double*[n];for(int i=0;i<n;i++) P0[i]=P+i*nP;
		fP(A0,nP,nQ,mdat0,n,nm,mpos0,dists_ch,dists_mid,dists_d,mid,nmid,P0,pp);
		delete[] A0; delete[] mdat0; delete[] mpos0; delete[] P0;
	}
}

/**************************************************************
# solve Ax =b by Cholesku decomposition  *
# A should be positive definite          *
#**************************************************************
mycholsl<- function(A,b){
  if(!is.matrix(A)) stop("mycholsl: A should be a matrix...")
  n<- dim(A)[1]
  if(n!=dim(A)[2]) stop("mycholsl: A should be a square matrix...") 
  x<- b 
  out<- .C("cholsolve",
          as.double(t(A)),
          as.integer(n),
          as.double(b),
          x=as.double(x))$x
  as.vector(out)
}

#****************************************************
# solve Ax =b by LU decomposition  *
# A should be positive definite    *
#****************************************************
mysolve<- function(A,b){
  if(!is.matrix(A)) stop("mycholsl: A should be a matrix...")
  n<- dim(A)[1]
  if(n!=dim(A)[2]) stop("mycholsl: A should be a square matrix...") 
  x<- b 
  out<- .C("lusolve",
          as.double(t(A)),
          as.integer(n),
          as.double(b),
          x=as.double(x))$x
  as.vector(out)
}
 
#******************************************************************
# generate n random variable from (0,1)     #
# seed: system time will be used if missing #
#******************************************************************
myrunif<- function(n,seed){
  if(missing(seed)) seed<- 0
  a<- rep(0.0,n)
  a<- .C("runifc",
    x=as.double(a),
    as.integer(n),
    as.integer(seed))$x
  a
}

#***********************************************************
# cholesky decomposition of n by n matrix a #
#***********************************************************
mychol<- function(a){
  n1<- nrow(a)
  n2<- ncol(a)
  if(n1!=n2) stop("cholesky decomposition: should be a square matrix...")
  p<- rep(0,n1)
  out<- .C("cholc",
    l=as.double(t(a)),
    as.integer(n1),
    p=as.double(p))
  L<- matrix(out$l,nrow=n1,byrow=T)
  p<- out$p
  for(i in 1:n1){
    for(j in 1:n2){
      if(j>i){
        L[i,j]<- 0
      }else if(j==i) L[i,j]<- p[i]
    }
  }
  L
}

#***********************************************************
# svd decomposition of n by m matrix a #
#***********************************************************
mysvd<- function(a){
  n1<- nrow(a)
  n2<- ncol(a)
  nr<- max(n1,n2)
  nc<- min(n1,n2)
  w<- rep(0,nc)
  if(n1>=n2) a<- t(a)
  v<- matrix(0,nrow=nc,ncol=nc)
  out<- .C("svdc",
    u=as.double(a),
    as.integer(nr),
    as.integer(nc),
    w=as.double(w),
    v=as.double(v))
  o<- list()
  o$d<- out$w
  if(n1>=n2){
    o$u<- matrix(out$u,nrow=nr,ncol=nc,byrow=T)
    o$v<- matrix(out$v,nrow=nc,ncol=nc,byrow=T)
  }else{
    o$u<- matrix(out$v,nrow=nc,ncol=nc,byrow=T)
    o$v<- matrix(out$u,nrow=nr,ncol=nc,byrow=T)
  }
  o
}

***********************************************************/
