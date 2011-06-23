//sure.c

//#include "mtcmim.h"

/*------------------------------------------------------------
 ML estimate b (of length kk+p) and sigma (p by p) by SURE
 y: n by p data matrix
 x: n by m data matrix (predictors)
 nqs (of length p): how many predictors to fit for y_j
 qs (of length kk): which predictors to fit for y_j
 ini_sigma: 1 if sigma is initial, 0 if not
 Return: log-likelihood
 ------------------------------------------------------------*/
double sureEst(double** y,int n,int p,double** x,int m,int* nqs,int* qs,
	double* b,double** sigma,int ini_sigma,int iter,double tol){
	#define INF 1e+308
	#define PI 3.141592653589793
	const double CONST = -n*p/2.0*(log(2.0*PI)+1.0);
	double lik, lik1, lik2;
	double la, la1;

	int kk = sum(nqs,p);
	int its = 2;
	int idx;
  	double mx;
  	double tmp;
	double* ybar=new double[p];
	double** e=new double*[n]; for(int i=0;i<n;i++) e[i]=new double[p];
	double* b0=new double[m+1];
	double* b1=new double[m+1];
	double* bb=new double[kk+p];
	double** sigma0=new double*[p]; for(int i=0;i<p;i++) sigma0[i]=new double[p];
	double** sinv=new double*[p]; for(int i=0;i<p;i++) sinv[i]=new double[p];
	double** den=new double*[m+1]; for(int i=0;i<m+1;i++) den[i]=new double[m+1];
	for(int j=0;j<p;j++){
		ybar[j] = 0.0;
		for(int i=0;i<n;i++){
			ybar[j] += y[i][j];
		}
		ybar[j] /= n;
	}
	for(int i=0;i<n;i++){
		for(int j=0;j<p;j++){
			e[i][j] = y[i][j]-ybar[j];
		}
	}
	if(ini_sigma==0){
		for(int j=0;j<p;j++){
			for(int k=j;k<p;k++){
				sigma[j][k] = 0.0;
				for(int i=0;i<n;i++){
					sigma[j][k] += e[i][j]*e[i][k];
				}
				sigma[j][k] /= n;
			}
			for(int k=0;k<j;k++){
				sigma[j][k] = sigma[k][j];
			}
		}
	}else if(ini_sigma==1){
		inv(sigma,p,sinv);
		while(its>0){
			idx = 0;
			for(int k=0;k<kk+p;k++){
				bb[k] = b[k];
			}
			for(int k=0;k<p;k++){
				den[0][0] = n*sinv[k][k];
				for(int s=1;s<=nqs[k];s++){
					den[0][s] = 0.0;
					for(int i=0;i<n;i++){
						den[0][s] += x[i][qs[idx+s-1]-1];
					}
					den[0][s] *= sinv[k][k];
				}
				for(int r=1;r<=nqs[k];r++){
					for(int s=r;s<=nqs[k];s++){
						den[r][s] = 0.0;
						for(int i=0;i<n;i++){
							den[r][s] += x[i][qs[idx+r-1]-1]*x[i][qs[idx+s-1]-1];
						}
						den[r][s] *= sinv[k][k];
					}
					for(int s=0;s<r;s++){
						den[r][s] = den[s][r];
					}
				}
				for(int r=0;r<nqs[k]+1;r++){
					b0[r] = 0.0;
					for(int i=0;i<n;i++){
						tmp = 0.0;
						for(int j=0;j<p;j++){
							if(j!=k) tmp += sinv[k][j]*e[i][j];
							else tmp += sinv[k][j]*y[i][j];
						}
						if(r==0) b0[r] += tmp;
						else b0[r] += tmp*x[i][qs[idx+r-1]-1];
					}
				}
				cholsl(den,nqs[k]+1,b0,b1);
				//solve(den,nqs[k]+1,nqs[k]+1,b0,b1);
				for(int r=0;r<nqs[k]+1;r++){
					b[idx+k+r] = b1[r];
				}
				for(int i=0;i<n;i++){
					e[i][k] = y[i][k]-b[idx+k];
					for(int r=0;r<nqs[k];r++){
						e[i][k] -= x[i][qs[idx+r]-1]*b[idx+k+r+1];
					}
				}
				idx += nqs[k];
			}
			mx = 0.0;
			for(int k=0;k<kk+p;k++){
				tmp = abs(b[k]-bb[k]);
				if(mx<tmp) mx = tmp;
			}
			its--;
		}
	}else{
      cout<<"ini_sigma="<<ini_sigma<<endl;
		error(_("sureEst: wrong ini_sigma ...\n"));
	}
	for(int j=0;j<p;j++){
		for(int k=0;k<p;k++){
			sigma0[j][k] = sigma[j][k];
		}
	}
	lik = CONST - n/2.0*log(det(sigma,p));
	lik1 = -INF;
	lik2 = -INF;
	la = INF;
	do{
		idx = 0;
		inv(sigma,p,sinv);
		for(int k=0;k<p;k++){
			den[0][0] = n*sinv[k][k];
			for(int s=1;s<=nqs[k];s++){
				den[0][s] = 0.0;
				for(int i=0;i<n;i++){
					den[0][s] += x[i][qs[idx+s-1]-1];
				}
				den[0][s] *= sinv[k][k];
			}
			for(int r=1;r<=nqs[k];r++){
				for(int s=r;s<=nqs[k];s++){
					den[r][s] = 0.0;
					for(int i=0;i<n;i++){
						den[r][s] += x[i][qs[idx+r-1]-1]*x[i][qs[idx+s-1]-1];
					}
					den[r][s] *= sinv[k][k];
				}
				for(int s=0;s<r;s++){
					den[r][s] = den[s][r];
				}
			}

			for(int r=0;r<nqs[k]+1;r++){
				b0[r] = 0.0;
				for(int i=0;i<n;i++){
					tmp = 0.0;
					for(int j=0;j<p;j++){
						if(j!=k) tmp += sinv[k][j]*e[i][j];
						else tmp += sinv[k][j]*y[i][j];
					}
					if(r==0) b0[r] += tmp;
					else b0[r] += tmp*x[i][qs[idx+r-1]-1];
				}
			}
			cholsl(den,nqs[k]+1,b0,b1);
			//solve(den,nqs[k]+1,nqs[k]+1,b0,b1);
			for(int r=0;r<nqs[k]+1;r++){
				b[idx+k+r] = b1[r];
			}
			for(int i=0;i<n;i++){
				e[i][k] = y[i][k]-b[idx+k];
				for(int r=0;r<nqs[k];r++){
					e[i][k] -= x[i][qs[idx+r]-1]*b[idx+k+r+1];
				}
			}
			idx += nqs[k];
		}
		for(int j=0;j<p;j++){
			for(int k=j;k<p;k++){
				sigma[j][k] = 0.0;
				for(int i=0;i<n;i++){
					sigma[j][k] += e[i][j]*e[i][k];
				}
				sigma[j][k] /= n;
			}
			for(int k=0;k<j;k++){
				sigma[j][k] = sigma[k][j];
			}
		}
		lik2 = lik1;
		lik1 = lik;
		lik = CONST - n/2.0*log(det(sigma,p));
		if(lik==lik1) break;
		la1 = la;
		la = lik1+(lik-lik1)/(1-(lik-lik1)/(lik1-lik2));
		mx = abs(la-la1);
		iter--;
		if(iter<0){
			cout<<"sureEst: convergence failed..."<<endl;
			break;
		}
	}while(mx>tol);
 	
	for(int i=0;i<n;i++) delete[] e[i];
	for(int i=0;i<p;i++) delete[] sigma0[i];
	for(int i=0;i<p;i++) delete[] sinv[i];
	for(int i=0;i<m+1;i++) delete[] den[i];
	delete[] ybar; delete[] e; delete[] b0; delete[] b1; delete[] bb;
	delete[] sigma0; delete[] sinv; delete[] den;

	return lik;
}

/*------------------------------------------------------------------
 Select the predictor  to add to the SURE model by ML
 y: n by p data matrix
 x: n by m data matrix (predictors)
 nqs (of length p): how many predictors to fit for y_j
 qs (of length kk): which predictors to fit for y_j
 nupper (of length p): how many predictors to choose to fit for y_j
 upper: which predictors to choose from to fit for y_j
 sigma: will change
 which[2]: which predictor to be selected
 	 which[0]--which y_i, which[1]--which x_j
 ------------------------------------------------------------------*/
bool sureAdd1(double** y,int n,int p,double** x, int m,int* nqs,int* qs,
	int* nupper,int* upper,double** sigma,int* which, int iter=100, double tol=1e-8){
	int idx1=0;
	int idx2=0;
	bool no=true;
	for(int i=0;i<p;i++){//check if out of upper
		for(int j=idx1;j<idx1+nqs[i];j++){
			for(int k=idx2;k<idx2+nupper[i];k++){
				if(qs[j]==upper[k]) {no = false; break;}
			}
			if(no){
				error(_("sure add1: out of upper...\n"));
			}
		}
		idx1 += nqs[i];
		idx2 += nupper[i];
	}
	
	#define PI 3.141592653589793
	#define INF 1e+308
//	const double CONST = -n*p/2.0*(log(2.0*PI)+1.0);
	bool add = false;
	int kk=sum(nqs,p); //total number of predictors in model
	
	int idx;
	int nl;
	int* nqs1=new int[p];
	int* qs1=new int[kk+1];
	double* b1=new double[kk+p+1];
	double** sigma1=new double*[p]; for(int i=0;i<p;i++) sigma1[i]=new double[p];
	int* vl = new int[m]; 
	for(int r=0;r<p;r++){
		for(int c=0;c<p;c++){
			sigma1[r][c] = sigma[r][c];
		}
	}
	double lik;
	double likold = -INF;
	idx1=0; idx2=0;
	for(int i=0;i<p;i++){
		nl = nupper[i]-nqs[i]; //number of predictors left
		if(nl>0){
			add = true;
			for(int j=0;j<p;j++){
				if(j==i) nqs1[j] = nqs[j]+1;
				else nqs1[j] = nqs[j];
			}
			idx = 0;
			for(int k=0;k<nupper[i];k++){
				no = true;
				for(int j=0;j<nqs[i];j++){
					if(upper[idx2+k]==qs[idx1+j]){no = false; break;}
				}
				if(no){vl[idx] = upper[idx2+k]; idx++;}
			}
			for(int j=0;j<idx1+nqs[i];j++){
				qs1[j] = qs[j];
			}
			for(int j=idx1+nqs[i]+1;j<kk+1;j++){
				qs1[j] = qs[j-1];
			}
			for(int k=0;k<nl;k++){
				qs1[idx1+nqs[i]] = vl[k];
				lik = sureEst(y,n,p,x,m,nqs1,qs1,b1,sigma1,1,iter,tol);
//				lik = CONST - n/2.0*log(det(sigma1,p));
				if(lik>likold+1e-8){
					which[0] = i;
					which[1] = vl[k];
					likold = lik;
					arr_copy(sigma1,p,p,sigma);
				}
			}
		}
		idx1 += nqs[i];
		idx2 += nupper[i];
	}

	for(int i=0;i<p;i++) delete[] sigma1[i];
    delete[] nqs1; delete[] qs1; delete[] b1; delete[] sigma1;
	delete[] vl;
	return add;
}
	
/*------------------------------------------------------------------
 Select the predictor  to drop from the SURE model by ML
 y: n by p data matrix
 x: n by m data matrix (predictors)
 nqs (of length p): how many predictors to fit for y_j
 qs (of length kk): which predictors to fit for y_j
 sigma: will change
 which[2]: which predictor to be selected
 	 which[0]--which y_i, which[1]--which x_j
 ------------------------------------------------------------------*/
bool sureDrop1(double** y,int n,int p,double** x, int m,int* nqs,int* qs,
	int* nlower,int* lower,double** sigma,int* which,int iter=100, double tol=1e-8){
	int idx1=0;
	int idx2=0;
	bool no=true;
	for(int i=0;i<p;i++){//check if include lower
		for(int j=idx1;j<idx1+nlower[i];j++){
			for(int k=idx2;k<idx2+nqs[i];k++){
				if(lower[j]==qs[k]) {no = false; break;}
			}
			if(no){
				error(_("sure drop1: not incude lower...\n"));
			}
		}
		idx1 += nlower[i];
		idx2 += nqs[i];
	}
	
	#define PI 3.141592653589793
	#define INF 1e+308
	bool drop = false;
//	const double CONST = -n*p/2.0*(log(2.0*PI)+1.0);
	int kk=sum(nqs,p); //total number of predictors in model
	if(kk<1) return drop;
	
	int* nqs1=new int[p];
	int* qs1=new int[kk-1];
	double* b1=new double[kk+p-1];
	double** sigma1=new double*[p]; for(int i=0;i<p;i++) sigma1[i]=new double[p];
	for(int r=0;r<p;r++){
		for(int c=0;c<p;c++){
			sigma1[r][c] = sigma[r][c];
		}
	}
	double lik;
	double likold = -INF;
	int idx;
	int nr;
	int* vr = new int[m]; 
	idx1=0; idx2=0;
	for(int i=0;i<p;i++){
		nr = nqs[i]-nlower[i];
		if(nr>0){
			drop = true;
			for(int j=0;j<p;j++){
				if(j==i) nqs1[j] = nqs[j]-1;
				else nqs1[j] = nqs[j];
			}
			idx = 0;
			for(int k=0;k<nqs[i];k++){
				no = true;
				for(int j=0;j<nlower[i];j++){
					if(qs[idx2+k]==lower[idx1+j]){no = false; break;}
				}
				if(no){vr[idx] = k; idx++;}
			}
			for(int j=0;j<nr;j++){
				for(int k=0;k<idx2+vr[j];k++){
					qs1[k] = qs[k];
				}
				for(int k=idx2+vr[j]+1;k<kk;k++){
					qs1[k-1] = qs[k];
				}
				lik = sureEst(y,n,p,x,m,nqs1,qs1,b1,sigma1,1,iter,tol);
//				lik = CONST - n/2.0*log(det(sigma1,p));
				if(lik>likold+1e-8){
					which[0] = i;
					which[1] = qs[idx2+vr[j]];
					likold = lik;
					arr_copy(sigma1,p,p,sigma);
				}
			}
		}
		idx1 += nlower[i];
		idx2 += nqs[i];
	}
	
	for(int i=0;i<p;i++) delete[] sigma1[i];
    delete[] nqs1; delete[] qs1; delete[] b1; delete[] sigma1; delete[] vr;
	return drop;
}
	
/*------------------------------------------------------------------
 Model selection for the SURE model by ML
 y: n by p data matrix
 x: n by m data matrix (m predictors)
 nupper (of length p): how many predictors to fit for y_j
 upper (of length kk): which predictors to fit for y_j
 k: penalty for each parameter
 direction: 0 if backward, 1 if forward, 2 if both
 vin: p by m matrix (0/1); 1 if selected, 0 if not selected
 	input initial model and return the final model
 record: return y_j for which to add (+) or drop (-)
 	x_j to add or drop, and loglikelihood
 max_terms: forward possible only if number of terms < max_terms
 steps: max number of steps
 iter: max number of iterations in estimation
 tol: accuracy of parameter estimates
 ------------------------------------------------------------------*/
void myf1(int** vin,int p,int m,int* nxs){
	for(int i=0;i<p;i++){
		nxs[i] = 0;
		for(int j=0;j<m;j++){
			if(vin[i][j]==1) nxs[i]++;
			else if(vin[i][j]!=0){
				error(_("vin in sureStep: wrong info...\n"));
			}
		}
	}
}
void myf2(int** vin,int* nxs,int p,int m,int* xs){
	int idx=0;
	for(int i=0;i<p;i++){
		for(int j=0;j<m;j++){
			if(vin[i][j]==1) {xs[idx] = j+1; idx++;}
		}
	}
	if(idx!=sum(nxs,p)){
		error(_("Number of predictors: something was wrong...\n"));
	}
}

void sureStep(double** y,int n,int p,double** x, int m,int* nlower,int* lower,
	int* nupper,int* upper,double k,int direction,int** vin,double* record,
	int max_terms,int steps,int iter, double tol){
	bool backwd;
	bool forwd;
	if(direction == 0){backwd = true; forwd = false;}
	else if(direction == 1){backwd = false; forwd = true;}
	else if(direction == 2){backwd = true; forwd = true;}
	else {
		error(_("sureStep: wrong direction...\n"));
	}

	int* nxs=new int[p]; myf1(vin,p,m,nxs);
	int kk = sum(nxs,p);
	int* xs=new int[sum(nupper,p)+1]; myf2(vin,nxs,p,m,xs);
	
	const double CONST = -n*p/2.0*(log(2.0*PI)+1.0);
	double* b=new double[kk+p];
	double** sigma=new double*[p]; for(int i=0;i<p;i++) sigma[i]=new double[p];
	double** sigma1=new double*[p]; for(int i=0;i<p;i++) sigma1[i]=new double[p];
	sureEst(y,n,p,x,m,nxs,xs,b,sigma,0,iter,tol);

	double lik;
	double lik1;
	double aic;
	double aic1;
	int idx = 0;
	lik = CONST - n/2.0*log(det(sigma,p));
	aic = -2*lik+k*(kk+p+p*(p+1)/2);
	record[idx] = 0; idx++;
	record[idx] = 0; idx++;
	record[idx] = lik; idx++;
	bool go;
	bool yes;
	int* which=new int[2];
	while(steps>0){
		go = false;
		if(backwd){
			yes = false;
			arr_copy(sigma,p,p,sigma1);
			yes = sureDrop1(y,n,p,x,m,nxs,xs,nlower,lower,sigma1,which,iter,tol);
			if(yes){
				lik1 = CONST - n/2.0*log(det(sigma1,p));
				aic1 = -2*lik1+k*(kk+p-1+p*(p+1)/2);
				if(aic1<aic-1e-8){
					vin[which[0]][which[1]-1] = 0;
					nxs[which[0]]--;
					myf2(vin,nxs,p,m,xs); 
					arr_copy(sigma1,p,p,sigma);
					lik = lik1;
					aic = aic1;
					kk--;
					go = true;
					record[idx] = -(which[0]+1); idx++; //"-" if drop
					record[idx] = which[1]; idx++;
					record[idx] = lik; idx++;
				}
			}
		}
		if(forwd && kk<max_terms){
			yes = false;
			arr_copy(sigma,p,p,sigma1);
			yes = sureAdd1(y,n,p,x,m,nxs,xs,nupper,upper,sigma1,which,iter,tol);
			if(yes){
				lik1 = CONST - n/2.0*log(det(sigma1,p));
				aic1 = -2*lik1+k*(kk+1+p+p*(p+1)/2);
				if(aic1<aic-1e-8){
					vin[which[0]][which[1]-1] = 1;
					nxs[which[0]]++;
					myf2(vin,nxs,p,m,xs);
					arr_copy(sigma1,p,p,sigma);
					lik = lik1; 
					aic = aic1;
					kk++;
					go = true;
					record[idx] = +(which[0]+1); idx++; //"+" if add
					record[idx] = which[1]; idx++;
					record[idx] = lik; idx++;
				}
			}
		}
		if(!go) break;
		steps--;
	}
	record[idx] = 9999; idx++; //end with 9999
	record[idx] = 9999; idx++;
	record[idx] = 9999; idx++;

    for(int i=0;i<p;i++){delete[] sigma[i]; delete[] sigma1[i];}
	delete[] nxs; delete[] xs; delete[] which;
	delete[] b; delete[] sigma; delete[] sigma1;
}

