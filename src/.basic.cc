//basic.c

//#include "mtcmim.h"

/*----------------------------------
 transpose of m by n array arr
 ----------------------------------*/
void arr_t(double** arr,int m, int n,double** result){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			result[j][i]=arr[i][j];
		}
	}
}

/*---------------------------------------------------
 Cholesky decomposition A = L*LT (positive definite)
 A: input (only the upper triangle is needed,
 		which is not modified), L is returned by the 
 		lower triangle of A and p
 p: vector of length n, contains the diagonal of L      
 --------------------------------------------------*/
void chol(double** A, int n, double* p){
	int i,j,k;
	double sum;
	for (i=0;i<n;i++) {
		for (j=i;j<n;j++) {
			for (sum=A[i][j],k=i-1;k>=0;k--) sum -= A[i][k]*A[j][k];
			if (i == j) {
				if (sum <= 0.0){
					error(_("Cholesky decomposition failed...\n"));
				}
				p[i]=sqrt(sum);
			} else A[j][i]=sum/p[i];
		}
	}
}

/*--------------------------------------------
 Solve A *x = b, where A is positive definite    
 ---------------------------------------------*/
void cholsl(double **A, int n, double b[], double x[]){
	double* p = new double[n];
	chol(A,n,p);
	
	int i,k;
	double sum;
	for (i=0;i<n;i++) {
		for (sum=b[i],k=i-1;k>=0;k--) sum -= A[i][k]*x[k];
		x[i]=sum/p[i];
	}
	for (i=n-1;i>=0;i--) {
		for (sum=x[i],k=i+1;k<n;k++) sum -= A[k][i]*x[k];
		x[i]=sum/p[i];
	}
	delete[] p;
}

/*-------------------------------------
 dnorm: normal density function
 defualt: mean=0 and sd=1
 --------------------------------------*/
double dnorm(double x, double mean, double sd){
	# define PI 3.141592653589793
	double pdf;
	double y;

	y = ( x - mean ) / sd;
	pdf = exp ( - 0.5*y*y )	/ ( sd * sqrt ( 2.0 * PI ) );

	return pdf;
	# undef PI
}

/*------------------------
 create a file name 
 ------------------------*/
void fchar(char *path,char *file,int num,char *type,char *buff){
	if(num==0){
		strcpy(buff,path);
		strcat(buff,file);
		strcat(buff,type);
	}
	else{
		char* p=new char[100];
		strcpy(buff,path);
		strcat(buff,file);
		itoa(num,p);
		strcat(buff,p);
		strcat(buff,type);
		delete[] p;
	}
}

/*--------------------------------------------------
 Determinant of n by n matrix A by LU decomposition
 A: n by n matrix, not change
 --------------------------------------------------*/
double det(double** A, int n){
	double** a=new double*[n]; for(int i=0;i<n;i++) a[i]=new double[n];
	int* indx=new int[n];
	double d;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			a[i][j] = A[i][j];
		}
	}
	lud(a,n,indx,&d);
	for(int j=0;j<n;j++) d *= a[j][j];
	
	for(int i=0;i<n;i++) delete[] a[i];
	delete[]a; delete[] indx;
	
	return d;
}

/*--------------------------------
 Generalized inverse of A by SVD
 A: m by n matrix
 ---------------------------------*/
void ginv(double** A, int m, int n, double** ginvA){
	double s;
	double** u=new double*[m]; for(int i=0;i<m;i++) u[i]=new double[n];
	double* w=new double[n];
	double** v=new double*[n]; for(int i=0;i<n;i++) v[i]=new double[n];
	double eps = numeric_limits<double>::epsilon();

	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			u[i][j] = A[i][j];
		}
	}
	svd(u, m, n, w, v);	
	for (int j=0;j<n;j++) {
		if(abs(w[j])>sqrt(eps)) s = 1/w[j];
			else s = 0.0;
		for(int i=0;i<n;i++){
			v[i][j] *= s;
		}
	}
	for(int i=0;i<n;i++){
		for (int j=0;j<m;j++) {
			ginvA[i][j] = 0.0;
			for(int k=0;k<n;k++){
				ginvA[i][j] += v[i][k]*u[j][k];
			}
		}
	}

	for(int i=0;i<m;i++) delete[] u[i];
	for(int i=0;i<n;i++) delete[] v[i];
	delete[] u; delete[] w; delete[] v;
}

/*-----------------------------------------------
 Inverse of n by n matrix A by LU decomposition
 A: n by n matrix, not change
 -----------------------------------------------*/
void inv(double** A, int n, double** invA){
	double** a=new double*[n]; for(int i=0;i<n;i++) a[i]=new double[n];
	int* indx=new int[n];
	double* col=new double[n];
	double d;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			a[i][j] = A[i][j];
		}
	}
	lud(a,n,indx,&d);
	for(int j=0;j<n;j++) {
		for(int i=0;i<n;i++) col[i]=0.0;
		col[j]=1.0;
		lubksb(a,n,indx,col);
		for(int i=0;i<n;i++) invA[i][j]=col[i];
	}

	for(int i=0;i<n;i++) delete[] a[i];
	delete[]a; delete[] indx; delete[] col;
}

/*--------------------------------------------------------------
 Inverse and determinant of n by n matrix A by LU decomposition
 A: n by n matrix, not change
 --------------------------------------------------------------*/
double inv_det(double** A, int n, double** invA){
	double** a=new double*[n]; for(int i=0;i<n;i++) a[i]=new double[n];
	int* indx=new int[n];
	double* col=new double[n];
	double d;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			a[i][j] = A[i][j];
		}
	}
	lud(a,n,indx,&d);
	for(int j=0;j<n;j++) {
		d *= a[j][j];
		for(int i=0;i<n;i++) col[i]=0.0;
		col[j]=1.0;
		lubksb(a,n,indx,col);
		for(int i=0;i<n;i++) invA[i][j]=col[i];
	}

	for(int i=0;i<n;i++) delete[] a[i];
	delete[]a; delete[] indx; delete[] col;

	return d;
}

/*----------------------------------------------
 itoa: convert an integer into a string
 default: base=10
 -----------------------------------------------*/
void itoa(int i,char buff[],int base){
	int len=100;
	char* tmp=new char[len];
	int k,l=i;
	int indx=0;
	while(l){
		k=l%base;
		tmp[indx]=48+k;
		l=(l-k)/base;
		indx++;
	}
	for(int n=0;n<indx;n++){
		buff[n]=tmp[indx-n-1];
	}
	buff[indx]='\0';
	if(indx>len){
		error(_("Error: array provided is too short!\n"));
	}
	delete[] tmp;
}

/*-----------------------------------------------------
 LU decomposition: A = L*U (non-singluar)
 a: n by n matrix, will store L and U
 indx: vector of length n, records the row permutation
 	effected by the partial pivoting (0...n-1)
 d: output as ±1 depending on whether the number of row
   interchanges was even or odd, respectively
 -----------------------------------------------------*/
void lud(double **a, int n, int *indx, double *d){
	#define TINY 1.0e-38
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=new double[n];
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++){
			temp=abs(a[i][j]);
			if (temp > big) big=temp;
		}
		if (big == 0.0){
			error(_("Singular matrix in routine ludcmp...\n"));
		}
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		imax=j;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*abs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	delete[] vv;
}

/*-------------------------------
 Generate n uniform number (0,1)
 returned by x
 ------------------------------*/
double runif0(long *seed){
	#define IM1 2147483563
	#define IM2 2147483399
	#define AM (1.0/IM1)
	#define IMM1 (IM1-1)
	#define IA1 40014
	#define IA2 40692
	#define IQ1 53668
	#define IQ2 52774
	#define IR1 12211
	#define IR2 3791
	#define NTAB 32
	#define NDIV (1+IMM1/NTAB)
	#define EPS 1.2e-7
	#define RNMX (1.0-EPS)

	int j;
	long k;
	static long seed2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	if (*seed <= 0) {
		if (-(*seed) < 1) *seed=1;
			else *seed = -(*seed);
	}
	seed2=(*seed);
	for (j=NTAB+7;j>=0;j--) {
		k=(*seed)/IQ1;
		*seed=IA1*(*seed-k*IQ1)-k*IR1;
		if (*seed < 0) *seed += IM1;
		if (j < NTAB) iv[j] = *seed;
	}
	iy=iv[0];
	k=(*seed)/IQ1;
	*seed=IA1*(*seed-k*IQ1)-k*IR1;
	if (*seed < 0) *seed += IM1;
	k=seed2/IQ2;
	seed2=IA2*(seed2-k*IQ2)-k*IR2;
	if (seed2 < 0) seed2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-seed2;
	iv[j] = *seed;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
		else return temp;
}

void runif(double* x,int n,long *seed){
	if(*seed==0){
		*seed = (long int)time(0);
	}
	for(int i=0;i<n;i++){
		x[i] = runif0(seed);
	}
}

/*----------------------------------------------------------
 Solve Ax=b for x returned by b
 A: n by n matrix from LU decomposition, not change
 indx: vector of length n from LU decomposition, not change
 ----------------------------------------------------------*/
void lubksb(double **a, int n, int *indx, double b[]){
	int i,ii=0,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii-1;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void solve(double** A, int n, double b[], double x[]){
	double** a=new double*[n];for(int i=0;i<n;i++) a[i]=new double[n];
	int* indx=new int[n];
	double d;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			a[i][j] = A[i][j];
		}
		x[i] = b[i];
	}
	lud(a,n,indx,&d);
	lubksb(a,n,indx,x);

	for(int i=0;i<n;i++) delete[] a[i];
	delete[] a; delete[] indx;
}

/*--------------------------------
 Solve Ax=b for x returned by SVD
 A: m by n matrix (m >= n)
 Note: much slower than LU
 ---------------------------------*/
void solve(double** A, int m, int n, double b[], double x[]){
	int jj,j,i;
	double s,*tmp=new double[n];
	double** u=new double*[m];for(int i=0;i<m;i++) u[i]=new double[n];
	double* w=new double[n];
	double** v=new double*[n];for(int i=0;i<n;i++) v[i]=new double[n];
	double eps = numeric_limits<double>::epsilon();

	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			u[i][j] = A[i][j];
		}
	}
	svd(u, m, n, w, v);
	for (j=0;j<n;j++) {
		s=0.0;
		if (abs(w[j])>sqrt(eps)) {
			for (i=0;i<m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=0;j<n;j++) {
		s=0.0;
		for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}

	for(int i=0;i<m;i++) delete[] u[i];
	for(int i=0;i<n;i++) delete[] v[i];
	delete[] tmp; delete[] u; delete[] w; delete[] v;
}

/*--------------------------------------
 SVD decomposition: A = U*W*VT
 a: m by n matrix
 U: replace a on output
 W: vector of length n; singular values
 V: n by n matrix
 Note: not accurate if m<n
 --------------------------------------*/
double pythag(double a, double b){ 
	#define SQR(a) (a == 0.0 ? 0.0 : a*a)
	double absa,absb; 
	absa=abs(a); 
	absb=abs(b); 
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa)); 
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))); 
} 

void svd(double **a, int m, int n, double w[], double **v){ 
	#define SQR(a) (a == 0.0 ? 0.0 : a*a)
	#define FMAX(a,b) (a > b ? a : b)
	#define FMIN(a,b) (a < b ? a : b)
	#define SIGN(a,b) (b >= 0.0 ? abs(a) : -abs(a))
	#define eps numeric_limits<double>::epsilon()

	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1; 

	rv1=new double[n];
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i< m) {
			for (k=i;k<m;k++) scale += abs(a[k][i]);
			if (scale) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i< m && i!= n-1) {
			for (k=l;k<n;k++) scale += abs(a[i][k]);
			if (scale) {
				for (k=l;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k]; 
					for (k=l;k<n;k++) a[j][k] += s*rv1[k]; 
				} 
				for (k=l;k<n;k++) a[i][k] *= scale; 
			}
		}
		anorm=FMAX(anorm,(abs(w[i])+abs(rv1[i])));
	} 
	for (i=n-1;i>=0;i--) {
		if (i< n-1) { 
			if(g) { 
				for (j=l;j<n;j++)
				v[j][i]=(a[i][j]/a[i][l])/g; 
				for (j=l;j<n;j++) { 
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j]; 
					for (k=l;k<n;k++) v[k][j] += s*v[k][i]; 
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	} 
	for (i=FMIN(m,n)-1;i>=0;i--) {
		l=i+1; 
		g=w[i]; 
		for (j=l;j<n;j++) a[i][j]=0.0; 
		if (g) { 
			g=1.0/g; 
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			} 
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	} 
	for (k=n-1;k>=0;k--) {
		for (its=1;its<=250;its++) {
			flag=1; 
			for (l=k;l>=0;l--) {
				nm=l-1;
				if (abs(rv1[l]) < eps) { 
					flag=0; 
					break; 
				} 
				if (abs(w[nm]) < eps) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if (abs(f) < eps) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					} 
				} 
			} 
			z=w[k]; 
			if(l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				} 
				break;
			}
			if (its == 250){
				Rf_warning("svd: convergence might have failed in 250 iterations"); 
//				cout<<"\a	svd: convergence might have failed in 250 iterations"<<endl; 
//				exit(1);
			}
			x=w[l];
			nm=k-1; 
			y=w[nm]; 
			g=rv1[nm]; 
			h=rv1[k]; 
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y); 
			g=pythag(f,1.0); 
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x; 
			c=s=1.0;
			for (j=l;j<=nm;j++) { 
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y*=c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}

	delete[] rv1;
} 

/*--------------------------------------------
 print vetor arr into m rows and n columns
 default: width=12
 --------------------------------------------*/
template <class T>
void arr_print(T* arr,int m,int n,int width){
	cout<<endl;
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			cout<<setw(width)<<arr[i*n+j];
		}
		cout<<'\n';
	}
	cout<<endl;
}

/*--------------------------------------------
 print m by n arr into n columns
 default: width=12
 --------------------------------------------*/
template <class T>
void arr_print(T** arr,int m,int n,int width){
	cout<<endl;
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			cout<<setw(width)<<arr[i][j];
		}
		cout<<'\n';
	}
	cout<<endl;
}

/*--------------------------------------------------
 copy array1 into array (of length n)
 --------------------------------------------------*/
template <class T>
void arr_copy(T* arr1,int n,T* arr){
	for(int i=0;i<n;i++){
		arr[i]=arr1[i];
	}
}

/*--------------------------------------------------
 copy m by n array1 into m by n array
 --------------------------------------------------*/
template <class T>
void arr_copy(T** arr1,int m,int n,T** arr){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			arr[i][j]=arr1[i][j];
		}
	}
}

/*--------------------------------------------
 kronecker product of two arrays
 arr1: nr1 by nc1; arr2: nr2 by nc2
 arr: nr1*nr2 by nc1*nc2 to be returned
 --------------------------------------------*/
template <class T>
void arr_kronecker(T** arr1,int nr1,int nc1,T** arr2,int nr2,int nc2,T** arr){
	int i; int j;
	for(int i1=0;i1<nr1;i1++){
		for(int j1=0;j1<nc1;j1++){
			i=i1*nr2;
			for(int i2=0;i2<nr2;i2++){
				j=j1*nc2;
				for(int j2=0;j2<nc2;j2++){
					arr[i][j] = arr1[i1][j1]*arr2[i2][j2];
					j++;
				}
				i++;
			}
		}
	}
}

/*----------------------------------
 arr1: m by n; arr2: n by 1
 arr: m by 1 to be returned
 ----------------------------------*/
template <class T>
void arr_prod(T** arr1,T* arr2,int m,int n,T* arr){
	for(int i=0;i<m;i++){
		arr[i]=0;
		for(int l=0;l<n;l++)
			arr[i] += arr1[i][l]*arr2[l];
	}
}

/*----------------------------------
 arr1: m by k; arr2: k by n
 arr: m by n to be returned
 ----------------------------------*/
template <class T>
void arr_prod(T** arr1,T** arr2, int m,int k,int n,T** arr){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			arr[i][j]=0;
			for(int l=0;l<k;l++){
				arr[i][j] += arr1[i][l]*arr2[l][j];
			}
		}
	}
}

/*----------------------------------
 arr1: n by n; arr2: n by n
 arr: n by n to be returned
 ----------------------------------*/
template <class T>
void arr_prod(T** arr1,T** arr2, int n,T** arr){
	arr_prod(arr1, arr2, n, n, n, arr);
}

/*------------------------------------------------------------
 cumsum: cumulative sumation of an array of length n
 ------------------------------------------------------------*/
template <class T>
void cumsum(T* x,int n,T* arr){
	T s=0;
	for(int i=0;i<n;i++){
		s+=x[i];
		arr[i]=s;
	}
}

/*-----------------------------------------------
 largest: obtain largest k values
 from an array x of length n
 returned by arr (k<=n)
 choice: indices of elements being chosen
 default choice: NULL
 -----------------------------------------------*/
template <class T>
void largest(T *x,int n,T *arr,int k,int *which){
	if(k>n){
		error(_("largest: cann't get more numbers than given...\n"));
	}
	int* p=new int[n];
	order(x,n,p,false);
	if(which==NULL){
		for(int i=0;i<k;i++){
			arr[i]=x[p[i]];
		}
	}
	else{
		for(int i=0;i<k;i++){
			arr[i]=x[p[i]]; which[i]=p[i];
		}
	}
	delete[] p;
}

/*---------------------------------------------------------
 max: obtain the maximum of an array x of length n
 ---------------------------------------------------------*/
template <class T>
T max(T* x,int n,int* which){
	T m=x[0];
	if(which==NULL){
		for(int i=0;i<n;i++){
			if(x[i]>m)m=x[i];
		}
	}
	else{
		*which=0;
		for(int i=0;i<n;i++){
			if(x[i]>m){m=x[i]; *which=i;}
		}	
	}
	return(m);
}

/*---------------------------------------------------------
 meam: obtain the mean of an array x of length n
 ---------------------------------------------------------*/
template <class T>
double mean(T* x,int n){
	double m;
	m=static_cast<double>(sum(x,n))/n;
	return m;
}

/*---------------------------------------------------------
 min: obtain the minimum of an array x of length n
 ---------------------------------------------------------*/
template <class T>
T min(T* x,int n,int* which){
	T m=x[0];
	if(which==NULL){
		for(int i=0;i<n;i++){
			if(x[i]<m)m=x[i];
		}
	}
	else{
		*which=0;
		for(int i=0;i<n;i++){
			if(x[i]<m){m=x[i]; *which=i;}
		}
	}
	return(m);
}

/*----------------------------------------------------------
 order: obtain the order of an array x of length n
 returned by integeral arr of length n
 default: increasing=true
 ----------------------------------------------------------*/
template <class T>
void order(T* x, int n, int* arr, bool increasing){
	T tmp;
	int k;
	T* p=new T[n];

	for(int i=0;i<n;i++){
		p[i]=x[i];
		arr[i]=i;
	}
	if(increasing==true){
		for(int i=0;i<n-1;i++){
			for(int j=i+1;j<n;j++){
				if(p[j]<p[i]){
					tmp=p[i];
					p[i]=p[j];
					p[j]=tmp;
					k=arr[i];
					arr[i]=arr[j];
					arr[j]=k;
				}
			}
		}
	}
	else if(increasing==false){
		for(int i=0;i<n-1;i++){
			for(int j=i+1;j<n;j++){
				if(p[j]>p[i]){
					tmp=p[i];
					p[i]=p[j];
					p[j]=tmp;
					k=arr[i];
					arr[i]=arr[j];
					arr[j]=k;
				}
			}
		}		
	}
	delete[] p;
}

/*--------------------------------------------------------------
 quantile: obtain p-quantile of an array x of length n
 --------------------------------------------------------------*/
template <class T>
double quantile(T* x, int n, double p){
	if(p<0||p>1){
		error(_("misspecified p!\n"));
	}

	T* ptr=new T[n];
	sort(x,n,ptr);
	T q;
	int indx;
	double y1,y2;

	if(p<=1.0/n) q=ptr[0];
	else{
		for(int i=1;i<n;i++){
			if(p>(i+0.0)/n&&p<=(i+1.0)/n){
				indx=i; 
				break;
			}
		}
		y1=static_cast<double>(ptr[indx-1]); 
		y2=static_cast<double>(ptr[indx]);
		q=y1+(y2-y1)/((indx+1.0)/n-(indx+0.0)/n)*(p-(indx+0.0)/n);
	}

	delete[] ptr;
	return(q);
}

/*----------------------------------------------------
 read a distk file -- infile into an array x
 ----------------------------------------------------*/
template <class T>
void read_T(T* x,char infile[]){
	T* p=x;
	ifstream inf(infile);
	while(inf){
		inf>>*p;
		p++;
	}
	inf.close();
}

/*--------------------------------------------------------------
 rep: genarate an array of length n with same element x
 returned by arr of length n
 ---------------------------------------------------------------*/
template <class T>
void rep(T x,int n,T* arr){
	for(int i=0;i<n;i++){
		arr[i]=x;
	}
}

/*-------------------------------------------
 round: round a number x to n-th place
 default: n=0
 -------------------------------------------*/
template <class T>
double round(T x,int n){
	double y;
	double rd;
	rd=static_cast<double>(n);
	y=static_cast<double>(x);
	y=floor(y*pow(10,rd)+0.5)/pow(10,rd);
	return y;
}

/*-----------------------------------------------------------
 Sample: sample k observations from array x of n
 to array arr of k
 default: replace=false, seed=0 (time-dependent)
 Note: a permutation when k=n and replace=false
 ------------------------------------------------------------*/
template <class T>
void sample(T* x,int n,T* arr,int k,long int *seed,bool replace){
	if(replace==false){
		if(n<k){
			error(_("Cannot sample from a smaller group without replacement...\n"));
		}
		double *rs=new double[n];
		runif(rs,n,seed);

		int* ord=new int[n];
		order(rs,n,ord);
		for(int i=0;i<k;i++){
			arr[i]=x[ord[i]];
		}
		delete[] rs;
		delete[] ord;
	}
	else if(replace==true){
		int irs;
		double tmp;
		for(int i=0;i<k;i++){
			runif(&tmp,1,seed);
			irs=static_cast<int>(floor(tmp));
			arr[i]=x[irs];
		}
	}
}

/*---------------------------------------------------
 standard deviation of an array of length n
 ---------------------------------------------------*/
template <class T>
double sd(T* x,int n){
	double sdv;

	sdv=sqrt(var(x,n));
	return sdv;
}

/*-----------------------------------------------
 smallest: obtain smallest k values
 from an array x of length n
 returned by arr (k<=n)
 choice: indices of elements being chosen
 default choice: NULL
 -----------------------------------------------*/
template <class T>
void smallest(T *x,int n,int k,T *arr,int *choice){
	if(k>n){
		error(_("largest: cann't get more numbers than given...\n"));
	}
	int* p=new int[n];
	order(x,n,p,true);
	if(choice==NULL){
		for(int i=0;i<k;i++){
			arr[i]=x[p[i]];
		}
	}
	else{
		for(int i=0;i<k;i++){
			arr[i]=x[p[i]]; choice[i]=p[i];
		}
	}
	delete[] p;
}

/*----------------------------------------
 sort: sort an array x of length n
 returned by arr of length n
 default: increasing=true
 ----------------------------------------*/
template <class T>
void sort(T* x,int n,T* arr,bool increasing){
	T tmp;

	for(int i=0;i<n;i++){
		arr[i]=x[i];
	}
	if(increasing==true){
		for(int i=0;i<n-1;i++){
			for(int j=i+1;j<n;j++){
				if(arr[j]<arr[i]){
					tmp=arr[i];
					arr[i]=arr[j];
					arr[j]=tmp;
				}
			}
		}
	}
	else if(increasing==false){
		for(int i=0;i<n-1;i++){
			for(int j=i+1;j<n;j++){
				if(arr[j]>arr[i]){
					tmp=arr[i];
					arr[i]=arr[j];
					arr[j]=tmp;
				}
			}
		}		
	}
}

/*-------------------------------------------
 subset: select a subset of k elements
 by choice from an array x of length n
 returned by arr of length k
 -------------------------------------------*/
template <class T>
void subset(T* x,bool* choice,int n,T* arr){
	int indx=0;
	for(int i=0;i<n;i++){
		if(choice[i]==true){
			arr[indx]=x[i];
			indx++;
		}
	}
}

/*----------------------------------
	sum of an array of length n
 -----------------------------------*/
template <class T>
T sum(T* x,int n){
	T s=0;
	for(int i=0;i<n;i++){
		s+=x[i];
	}
	return s;
}

int sum(bool* x,int n){
	int s=0;
	for(int i=0;i<n;i++){
		if(x[i]==true) s++;
	}
	return(s);
}

/*-----------------------------------------
 variance of an array of length n
 -----------------------------------------*/
template <class T>
double var(T* x,int n){
	if(n<2){
		error(_("var not exits!\n"));
	}
	double m, v=0;

	m=mean(x,n);
	for(int i=0;i<n;i++){
		v += (x[i]-m)*(x[i]-m)/(n-1);
	}
	return v;
}

/*---------------------------------------------------
 write an array of length n 
 into a distk file (outfile) of ncol columns 
 ---------------------------------------------------*/
template <class T>
void write_T(T* x,int n,int ncol,char outfile[]){
	T* p=x;
	ofstream outf(outfile);
	int indx=1;
	for(int i=0;i<n;i++,indx++,p++){
		outf<<*p;
		if(indx%ncol==0){
			outf<<'\n';
		}
		else outf<<' ';
	}
	outf<<'\n';
	outf.close();
}

