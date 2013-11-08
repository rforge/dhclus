#include "def.h"

SEXP icav_c_silhouette(SEXP DMAT, SEXP LABELS, SEXP ULABELS, SEXP PAR){
// 	PAR[0] length of labels,PAR[1] = 0 AVG, ,PAR[1] = 1 silhouette vector;
	SEXP ret;
	double *dmat;
	double a,*b,*sil,bmin,den,*S,*S_cl,*INFO;
	
	int *labels,*ulabels,*par,ulen,*cl,cl_a,cl_b;
	
	int i,j,k,ii;
	char bool;
	PROTECT(DMAT = coerceVector(DMAT, REALSXP));
	PROTECT(LABELS = coerceVector(LABELS, INTSXP));
	PROTECT(ULABELS = coerceVector(ULABELS, INTSXP));
	PROTECT(PAR = coerceVector(PAR, INTSXP));
	
	dmat = REAL(DMAT); labels = INTEGER(LABELS); ulabels = INTEGER(ULABELS); ulen = length(ULABELS); par = INTEGER(PAR);
	
	PROTECT(ret = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ret, 0,  allocMatrix(REALSXP, par[0],3));
	SET_VECTOR_ELT(ret, 1,  allocVector(REALSXP, ulen));
	SET_VECTOR_ELT(ret, 2,  allocVector(REALSXP, 1));

	S = REAL(VECTOR_ELT(ret, 2));
	S_cl = REAL(VECTOR_ELT(ret, 1));
	INFO = REAL(VECTOR_ELT(ret, 0));
	sil = &(INFO[2*par[0]]);
	
	b = (double*) R_alloc(ulen,sizeof(double));
	cl = (int*) R_alloc(ulen,sizeof(int));
	
	for (j=0;j<ulen;j++) {cl[j]=0;S_cl[j]=0;}

	for (i=0;i<par[0];i++) for (j=0;j<ulen;j++) cl[j] += (ulabels[j]==labels[i]);
	for (i=0;i<par[0];i++){
		a = 0.;sil[i]=0.;
		ii = i*par[0];
		for (j=0;j<ulen;j++) b[j]=0.;
		INFO[i] = labels[i];
		for (k=0;k<par[0];k++){
			bool = (labels[i]==labels[k]);
			a += (bool) ? dmat[ii+k]:0.;
			if (bool==0) {for (j=0;j<ulen;j++) b[j] += (ulabels[j]==labels[k]) ?dmat[ii+k]:0.;}
		}
		for (j=0;j<ulen;j++) b[j] /= (cl[j]>1) ? cl[j]:1;
		bmin = MAXDOUBLE;
		for (j=0;j<ulen;j++) {
			bool = ((ulabels[j]!=labels[i])&&(bmin>b[j]));
			cl_a = (ulabels[j]==labels[i]) ? cl[j]:cl_a;
			INFO[i+par[0]] = bool ? ulabels[j]:INFO[i+par[0]];
			bmin = bool ? b[j]:bmin;
		}
		a /= (cl_a-1)>0 ? (cl_a-1):1;
		den = (bmin>a)? bmin:a ;den = (den==0)? 1:den ;
		sil[i] = (bmin-a)/den;
	}
	
	S[0] = 0.;
	for (i=0;i<par[0];i++) {
		S[0] += sil[i];
		for (j=0;j<ulen;j++) {S_cl[j] += (ulabels[j]==labels[i]) ? sil[i]:0;}
	}
	S[0] /= par[0];
	for (j=0;j<ulen;j++) S_cl[j] /= (cl[j]>1) ? cl[j]:1;
	
	UNPROTECT(5);
	return(ret);
}

SEXP icav_euc_metric(SEXP Data, SEXP par){//par[0] rows (dimension of Data), par[1] cols (amount of Data)

	double *d;
	int *p;

	int i,k,j,matlen,idx,kk,jj,kkk;
	double *m,dif;
	SEXP mimat;
	
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(par = coerceVector(par, INTSXP));
	
	d = REAL(Data);p = INTEGER(par);

	PROTECT(mimat = allocMatrix(REALSXP, p[1],p[1]));
	m = REAL(mimat);
	for (k=0;k<p[1];k++){
		m[k*p[1]+k] = 0;
		kk = k*p[0];kkk = k*p[1];
		for (j=k+1;j<p[1];j++){
			jj = j*p[0];
			dif = 0.;m[k*p[1]+j] = 0;
			for (i = 0;i<p[0];i++){
				dif = d[kk+i]- d[jj+i];
				m[kkk+j] += dif*dif;
			}
			m[j*p[1]+k] = m[kkk+j] = sqrt(m[kkk+j]);
		}
	}
	UNPROTECT(3);
	return(mimat);
}
