
//mappingBasic.c

//#include "mtcmim.h"

/////////////////////////////////////////////////////////
/*--------------------------------------------
 haldane: convert recombination rate c
 to distance d (M)
 -------------------------------------------*/
double haldane(double c){
	double d;
	d=-0.5*log(1-2*c);

	return d;
}

/*------------------------------------------
 haldane_inv: convert distance d (M)
 to recombination rate c
 ------------------------------------------*/
double haldane_inv(double d){
	double c;
	c=(1-exp(-2*d))/2;

	return c;
}

/*---------------------------------------------------------
  compute recombination rate ("a-b-c")
  no interference assumed
-----------------------------------------------------------*/
double rate_ac(double r_ab, double r_bc){
  double r_ac;
  r_ac=r_ab+r_bc-2*r_ab*r_bc;

  return r_ac;
}

double rate_bc(double r_ab,double r_ac){
  double r_bc;
  r_bc=(r_ac-r_ab)/(1-2*r_ab);

  return r_bc;
}

/*------------------------------------------------
  cal condl prob of QTL(QQ)
  m1 and m2: flanking markers that take 0 or 1
  r:m1-m2 and r1:m1-qtl
-------------------------------------------------*/
double getprob(int m1,int m2,double r1,double r){
  double p;
  double r2;
  r2=rate_bc(r1,r);
  if(m1==1&&m2==1){
    p=(1-r1)*(1-r2)/(1-r);
  }
  else if(m1==1&&m2==0){
    p=(1-r1)*r2/r;
  }
  else if(m1==0&&m2==1){
    p=r1*(1-r2)/r;
  }
  else if(m1==0&&m2==0){
    p=r1*r2/(1-r);
  }
  else{
    error(_("wrong marker information!\n"));
  }
  return p;
}

/*--------------------------------------------------------------
 cal condl prob of QTL(QQ)
 mdat: nrow by ncol matrix (0/1)
 mid: which column of mdata (or marker id)
 d1: distance to mid
 d: length of marker interval with left flanking marker id mid
 p: array of length nrow
---------------------------------------------------------------*/
double getp_1(int mL,int mR,double d1,double d){
  double r1,r;
  double p;
  r1=haldane_inv(d1/100);
  r=haldane_inv(d/100);
  p=getprob(mL,mR,r1,r);

  return p;
}

void getp(int** mdat,int nrow,int ncol,int mid,double d1,double d,double* p){
  for(int n=0;n<nrow;n++){
    p[n]=getp_1(mdat[n][mid-1],mdat[n][mid],d1,d);
  }
}

/*--------------------------------------------------------
 i:vector of marker interval index
 nmark:vector of marker numbers on chromosomes
 returned by array m of length i_len
---------------------------------------------------------*/
void itom(int *i,int i_len,int *nmark,int nm_len,int* m){
	int* ninterval=new int[nm_len];
	for(int j=0;j<nm_len;j++){
		ninterval[j]=nmark[j]-1;
	}
	int* cumsuminterval=new int[nm_len];
	cumsum(ninterval,nm_len,cumsuminterval);
	for(int j=0;j<i_len;j++){
		int indx=0;
		for(int k=0;k<nm_len;k++){
			if(i[j]<=cumsuminterval[k])break;
			indx++;
		}
		m[j]=i[j]+indx;
	}
	delete[] ninterval;
	delete[] cumsuminterval;
}

/*----------------------------------------------------------
 m:vector of marker index
 nmark:vector of marker numbers on chromosomes
 returned by array i of length m_len
-----------------------------------------------------------*/
void mtoi(int* m,int m_len,int* nmark,int nm_len,int* i){
	int* cumsumnmark=new int[nm_len];
	cumsum(nmark,nm_len,cumsumnmark);
	for(int j=0;j<m_len;j++){
		int indx=0;
		for(int k=0;k<nm_len;k++){
			if(m[j]<=cumsumnmark[k])break;
			indx++;
	 }
		i[j]=m[j]-indx;
	}
	delete[] cumsumnmark;
}

/*------------------------------------------------------------------
 read a distk file -- infile into a CLSmpos or CLSqtlpos 
 ------------------------------------------------------------------*/
template <class T> //T: mpos or qtlpos
void read_pos(T& pos,char infile[]){
	int n=pos.getn();
	ifstream inf(infile);
	for(int i=0;i<n;i++){
		inf>>pos.getid()[i]
		>>pos.getch()[i]
		>>pos.getm()[i]
		>>pos.getdist()[i];
	}
	inf.close();
}

/*------------------------------------------------------------------
 write a CLSmpos or CLSqtlpos into a distk file -- outfile 
 -------------------------------------------------------------------*/
template <class T> //T: mpos or qtlpos
void write_pos(T& pos,char outfile[]){
	int n=pos.getn();
	ofstream outf(outfile);
	for(int i=0;i<n;i++){
		outf<<pos.getid()[i]<<' '
		<<pos.getch()[i]<<' '
		<<pos.getm()[i]<<' '
		<<pos.getdist()[i]<<'\n';
	}
	outf<<'\n';
	outf.close();
}


