#include "def.h"

int icav_compare_doubles (const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;
     
	return (*da > *db) - (*da < *db);
}


SEXP icav_Matrix_chunk(SEXP _M, SEXP Col_index, SEXP Row_index,SEXP Param)
{
	//copies values in (Col_index and Row_index) of _M to a new matrix mimat. size of mimat < size of _M
	
	// param[0] rows of _M
	// param[1] cols of _M
	// param[2] chunk size
	SEXP mimat;
	int i,j,_COL;
	int len_row,len_col;
	
	double *M,*mm;
	int * col_index, * row_index , *param, *copy_col, *copy_row;

	PROTECT(_M = coerceVector(_M, REALSXP));
	PROTECT(Col_index = coerceVector(Col_index, INTSXP));
	PROTECT(Row_index = coerceVector(Row_index, INTSXP));
	PROTECT(Param = coerceVector(Param, INTSXP));
	
	len_row = length(Row_index);
	len_col = length(Col_index);
	PROTECT(mimat = allocMatrix(REALSXP,len_row,len_col));

	M = REAL(_M);col_index = INTEGER(Col_index); row_index = INTEGER(Row_index); param = INTEGER(Param);
	mm = REAL(mimat);

	copy_col = (double*) R_alloc(len_col,sizeof(double));
	copy_row = (double*) R_alloc(len_row,sizeof(double));
	
	//we order the cols
	for(i=0;i<len_col;i++) copy_col[i]=col_index[i]-1;
	for(j=0;j<len_row;j++) copy_row[j]=row_index[j]-1;
	for(i=0;i<len_col;i++){
	  _COL = copy_col[i]*param[0];
	  for(j=0;j<len_row;j++){
	    mm[j+len_row*i] = M[copy_row[j]+_COL];
	  }
	}
	
	UNPROTECT(5);
	return(mimat);
}

SEXP icav_Matrix_sort(SEXP _M,SEXP Param)
{
	// sorts a cloned version
	//copies values in (Col_index and Row_index) of _M to a new matrix mimat. size of mimat < size of _M
	
	// param[0] rows of _M
	// param[1] cols of _M
	// param[2] chunk size
	SEXP mimat;
	int i,j,_COL;
	int MAX_T,Tnum;
	int one = 1;
	
	double *M,*mm;
	int *param;

	PROTECT(_M = coerceVector(_M, REALSXP));
	PROTECT(Param = coerceVector(Param, INTSXP));

	M = REAL(_M); param = INTEGER(Param);
	PROTECT(mimat = allocMatrix(REALSXP,param[0],param[1]));
	mm=REAL(mimat);
  #ifdef USE_OPENMP    // md EEV 10-11-2013  
	MAX_T = omp_get_max_threads();
#endif
	//we order the cols
	for(i=0;i<param[1];i++){
	  _COL = i*param[0];
	  for(j=0;j<param[0];j++){
	    mm[j+_COL]=M[j+_COL];
	  }
	}
	for(i=0;i<param[1];i++){
	  _COL = i*param[0];
	  qsort (&mm[_COL], param[0], sizeof(double), icav_compare_doubles);
	}
	
	UNPROTECT(3);
	return(mimat);
}


// used in gap.inter.R
SEXP icav_wss_matrix(SEXP _M, SEXP Param)
{
	//copies values in (Col_index and Row_index) of _M to a new matrix mimat. size of mimat < size of _M
	
	// param[0] rows of _M
	// param[1] cols of _M
//param[0] == param[1] this matrix must be symmetric
	SEXP mimat;
	int i,j,_COL;
	int MAX_T,Tnum;
	
	double *M,*mm,tmp_max;
	double tmp_acc;
	int *param;

	PROTECT(_M = coerceVector(_M, REALSXP));
	PROTECT(Param = coerceVector(Param, INTSXP));

	M = REAL(_M); param = INTEGER(Param);

	PROTECT(mimat = allocVector(REALSXP, 1));

	mm = REAL(mimat);
		tmp_acc = 0.;
	//we order the cols
	for(i=0;i<param[1];i++){
	  _COL = i*param[0];
	  for(j=0;j<param[0];j++){
		  tmp_acc += M[_COL+j]*M[_COL+j];
	  }
	}

	*mm = tmp_acc/(2*param[0]);
	UNPROTECT(3);
	return(mimat);
}
