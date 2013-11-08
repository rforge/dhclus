#include "def.h"

SEXP icav_sq_euc_metric(SEXP Data, SEXP par){//par[0] rows (dimension of Data), par[1] cols (amount of Data)

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
			m[j*p[1]+k] = m[kkk+j];
		}
	}
	UNPROTECT(3);
	return(mimat);
}


SEXP icav_gap_withinSum(SEXP DMAT, SEXP LABELS, SEXP ULABELS, SEXP PAR){
// 	PAR[0] length of labels, DMAT is PAR[0]xPAR[0]
	SEXP ret;
	double *dmat,*_ret, d_zero=0.;
	int *labels,*ulabels,*par,ulen,*cl;
	int i,j,current_cluster,COLUMN,bogus=-1;
	char bool;
	PROTECT(DMAT = coerceVector(DMAT, REALSXP));
	PROTECT(LABELS = coerceVector(LABELS, INTSXP));
	PROTECT(ULABELS = coerceVector(ULABELS, INTSXP));
	PROTECT(PAR = coerceVector(PAR, INTSXP));
	
	dmat = REAL(DMAT); labels = INTEGER(LABELS); ulabels = INTEGER(ULABELS); ulen = length(ULABELS); par = INTEGER(PAR);
	
	PROTECT(ret = allocVector(REALSXP,  ulen));
	
	_ret = REAL(ret);
	cl = (int*) R_alloc(ulen,sizeof(int));
	
	for (j=0;j<ulen;j++) {_ret[j]=0.;cl[j]=0;}

	for (i=0;i<par[0];i++) for (j=0;j<ulen;j++) cl[j] += (ulabels[j]==labels[i]);
	for (i=0;i<par[0];i++){
		COLUMN = i*par[0];
		current_cluster= bogus;
		for (j=0;j<ulen;j++) current_cluster = (ulabels[j]==labels[i])? j : current_cluster;
		
		for (j=i+1;j<par[0];j++){
			_ret[current_cluster] += (ulabels[current_cluster]==labels[j]) ? dmat[j+COLUMN]*dmat[j+COLUMN] : d_zero;
		}
	}
	for (j=0;j<ulen;j++) _ret[j] /= (cl[j]>0) ? cl[j] : 1;
	
	
	UNPROTECT(5);
	return(ret);
}

// Same as above (withinSum) different implementation
SEXP icav_gap_withinSum_mean(SEXP Data, SEXP LABELS, SEXP ULABELS, SEXP PAR){
// 	PAR[0] length of labels, Data is PAR[0]xPAR[1]
	SEXP ret_list,mean,ret;
	double *dmat,*_ret, d_zero=0.,*_mean;
	int *labels,*ulabels,*par,ulen,*cl;
	int i,j,current_cluster;
	char bool;
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(LABELS = coerceVector(LABELS, INTSXP));
	PROTECT(ULABELS = coerceVector(ULABELS, INTSXP));
	PROTECT(PAR = coerceVector(PAR, INTSXP));
	
	dmat = REAL(Data); labels = INTEGER(LABELS); ulabels = INTEGER(ULABELS); ulen = length(ULABELS); par = INTEGER(PAR);
	
	PROTECT(ret_list = allocVector(VECSXP,  2));
	PROTECT(mean = allocMatrix(REALSXP, ulen,par[1] ));
	PROTECT(ret = allocVector(REALSXP, ulen));
	
	SET_VECTOR_ELT(ret_list, 0,  mean);
	SET_VECTOR_ELT(ret_list, 1,  ret);

_ret = REAL(ret);_mean = REAL(mean);

	cl = (int*) R_alloc(ulen,sizeof(int));
	
	for (j=0;j<ulen;j++) {_ret[j]=0.;cl[j]=0;}
	for (j=0;j<ulen*par[1];j++) _mean[j] = 0.; 

	// calculate number of elements in each cluster
	for (i=0;i<par[0];i++) for (j=0;j<ulen;j++) cl[j] += (ulabels[j]==labels[i]);

	// calculate mean
	for (i=0;i<par[0];i++){
		for (j=0;j<ulen;j++) current_cluster = (ulabels[j]==labels[i])? j : current_cluster;
		for (j=0;j<par[1];j++){
		  _mean[current_cluster+j*ulen] += (ulabels[current_cluster]==labels[i]) ? dmat[i+j*par[0]] : d_zero;
		}
	}
	for (i=0;i<ulen;i++) for (j=0;j<par[1];j++) _mean[i+j*ulen] /= (cl[i]>0) ? cl[i] : 1;

	// calculate Mean Square Dist. to the center of the cluster
	for (i=0;i<par[0];i++){
		for (j=0;j<ulen;j++) current_cluster = (ulabels[j]==labels[i])? j : current_cluster;
		for (j=0;j<par[1];j++){
		   _ret[current_cluster] +=  (dmat[i+j*par[0]]-_mean[current_cluster+j*ulen])*(dmat[i+j*par[0]]-_mean[current_cluster+j*ulen]);
		}
	}
	UNPROTECT(7);
	return(ret_list);
}
