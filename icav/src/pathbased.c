#include "def.h"

SEXP icav_mstPath(SEXP X, SEXP par){
	
	int i,j,k,OFFSET;
	int *p;
	char *visited,logic;
	int *merge,*count;
	double *diss,*dd,*D;
	SEXP ret,node,w;
	
	double dmin,*ws;
	int nmin,*nodo,cnt;
	
	PROTECT(par = coerceVector(par, INTSXP));
	PROTECT(X = coerceVector(X, REALSXP));
	PROTECT(ret = allocVector(VECSXP,2));
	p = INTEGER(par);	D = REAL(X);
	
	SET_VECTOR_ELT(ret, 0,  allocVector(VECSXP, p[0]));
	SET_VECTOR_ELT(ret, 1,  allocVector(VECSXP, p[0]));
	
	visited = (char*) R_alloc(p[0],sizeof(char));
	merge = (int*) R_alloc(p[0],sizeof(int));
	count = (int*) R_alloc(p[0],sizeof(int));
	diss = (double*) R_alloc(p[0],sizeof(double));
	dd = (double*) R_alloc(p[0],sizeof(double));
	
	memset (visited, 0, p[0]);
	visited[0] = 1;
	for (i=0;i<p[0];i++) {dd[i] = -1.;diss[i]=D[i],merge[i]=0;count[i]=1;}
	count[0]=0;
	for (j=1;j<p[0];j++){
		dmin = DBL_MAX;
		for (i=0;i<p[0];i++){
			logic = dmin>diss[i] && visited[i]==0;
			dmin = logic? diss[i]:dmin;
			nmin = logic? i:nmin;
		}
		OFFSET = nmin*p[0];
		visited[nmin] = 1;
		dd[nmin] = dmin;
		for (i=0;i<p[0];i++){
			logic = diss[i]>D[i+OFFSET] && visited[i]==0;
			diss[i] = logic? D[i+OFFSET]:diss[i];
			merge[i] = logic? nmin:merge[i];
		}
	}
	
	cnt = 0;
	for (i=1;i<p[0];i++) count[merge[i]]++;
	
	node = VECTOR_ELT(ret, 1);
	w = VECTOR_ELT(ret, 0);
	for (i=0;i<p[0];i++){
		SET_VECTOR_ELT(node, i,  allocVector(INTSXP, count[i]));
		SET_VECTOR_ELT(w, i,  allocVector(REALSXP, count[i]));
		nodo = INTEGER(VECTOR_ELT(node, i));
		ws = REAL(VECTOR_ELT(w, i));
		k = 0;
		for (j=1;j<p[0];j++) {
			logic = (merge[j]==i);
			nodo[k] =  logic ? j:nodo[k];
			ws[k] =   logic ? dd[j]:ws[k];
			k += logic ? 1:0;
		}
	}
	for (i=1;i<p[0];i++) {
		INTEGER(VECTOR_ELT(node, i))[count[i]-1] = merge[i];
		REAL(VECTOR_ELT(w, i))[count[i]-1] = dd[i];
	}
	UNPROTECT(3);
	return(ret);
}


void icav_pathBased(SEXP list, char* visited,int startpoint,double *diss,int *p, int *P, double *Pd, int *index){
	
	SEXP nodes_vec,nodes_list;
	SEXP diss_vec,diss_list;
	
	int len,i,idx,j;
	int *nodes;
	double *d;
	double val;
	
	nodes_list = VECTOR_ELT(list, 1);
	nodes_vec = VECTOR_ELT(nodes_list, startpoint);
	
	diss_list = VECTOR_ELT(list, 0);
	diss_vec = VECTOR_ELT(diss_list, startpoint);
	
	len = length(nodes_vec);
	nodes = INTEGER(nodes_vec);
	d = REAL (diss_vec);
	visited[startpoint] = 1;
	for (i=0;i<len;i++) {
		idx = (int)nodes[i];
		if (idx>-1 && visited[idx]==0){
			P[index[0]]=idx; Pd[index[0]] = d[i];index[0]++;
			val = 0;
			for(j = index[0]-2;j>-1;j--) {
				val = val < Pd[j+1] ? Pd[j+1]:val;
				diss[ P[index[0]-1]*p[0] + P[j] ] = diss[ P[j]*p[0] + P[index[0]-1] ] = val;
			}
			icav_pathBased(list,visited,idx,diss,p,P,Pd,index);
			index[0]--;
		}
	}
	return;
}

SEXP icav_pathDiss(SEXP list,SEXP par){
	
	SEXP DISS,nodes_list;
	double *diss;
	int i,j,*p;
	char *visited;
	int *P,index,*startPoints;
	double *Pd;
	
	PROTECT(par = coerceVector(par, INTSXP));
	p = INTEGER(par);
	
	nodes_list = VECTOR_ELT(list, 1);
	

	
	PROTECT(DISS = allocMatrix(REALSXP,p[0],p[0]));	

	visited = (char*) R_alloc(p[0],sizeof(char));
	P = (int*) R_alloc(p[0],sizeof(int));
	Pd = (double*) R_alloc(p[0],sizeof(double));
	diss = REAL(DISS);
	
	
	for (i=0;i<p[0];i++) {
		P[0] = i; Pd[0] = 0.;index = 1;
		memset (visited, 0, p[0]);
		icav_pathBased(list,visited,i,diss,p,P,Pd,&index);
	}
	UNPROTECT(2);
	return(DISS);
}

SEXP icav_pathDiss_1(SEXP list,SEXP par){
	
	SEXP DISS,nodes_list;
	double *diss;
	int i,j,*p;
	char *visited,logic;
	int *P,index,*startPoints;
	double *Pd;
	
	PROTECT(par = coerceVector(par, INTSXP));
	p = INTEGER(par);
	
	nodes_list = VECTOR_ELT(list, 1);
	

	
	PROTECT(DISS = allocMatrix(REALSXP,p[0],p[0]));	

	visited = (char*) R_alloc(p[0],sizeof(char));
	P = (int*) R_alloc(p[0],sizeof(int));
	Pd = (double*) R_alloc(p[0],sizeof(double));

	startPoints = (int*) R_alloc(p[0],sizeof(int));
	diss = REAL(DISS);
	
	j=0;for (i=0;i<p[0];i++) {
	  logic = length(VECTOR_ELT(nodes_list, i))==1;
	  startPoints[j] = logic ? i:startPoints[j];
	  j += logic;
	}
	for (i=0;i<p[0]*p[0];i++) diss[i] = -1.;
	for (i=0;i<(j-1);i++) {
		P[0] = startPoints[i]; Pd[0] = 0.;index = 1;
		memset (visited, 0, p[0]);
		icav_pathBased(list,visited,startPoints[i],diss,p,P,Pd,&index);
	}
	for (i=0;i<p[0];i++) diss[i+i*p[0]] = 0.;
	UNPROTECT(2);
	return(DISS);
}


void icav_reduceIndex(SEXP list, SEXP _len,SEXP _n){
  
  SEXP nodes_vec;
  int *nodes,*n,*len;
  int i,j;
  
  
  
  PROTECT(_len = coerceVector(_len, INTSXP));
  PROTECT(_n = coerceVector(_n, INTSXP));
  n = INTEGER(_n);len = INTEGER(_len);
	for (i=0;i<*len;i++){

	  nodes_vec = VECTOR_ELT(list, i);
	  nodes = INTEGER(nodes_vec);
	  for (j=0;j<length(nodes_vec);j++) nodes[j] -= *n; 
	}
	UNPROTECT(2);
}
