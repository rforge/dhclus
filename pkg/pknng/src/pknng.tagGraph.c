#include "def.h"

SEXP connect_kneigbours (SEXP list,SEXP par){
	/* par[0] len of data, par[1] dim of data ,p[2]conn type,p[3] 0 no pen, 1 exponencial, 2 lineal, 3 power, 4 exponencial (exp(x-1)) ,5 power (x^n)-1
	1 pen,p[4] metric,p[5] lineal constant*/
	
	
	//min_val changed for euclidean
	
	int *p,*var_i;
	lnode *root,*tmp,*min_root,*n,*m,*addroot;
	SEXP nodes,ret,grps,Dist;
	char *visited;
	int cnt,startpoint,grpups,e;
	int i,k,j;
	int from,to,min_from,min_to,bins,len_i,len_d;
	double dif,val,min_val,d_one=1.;
	char bool;

	PROTECT(par = coerceVector(par, INTSXP));
	
	p = INTEGER(par);
	e = p[5];
	nodes = VECTOR_ELT(list, 1);
	visited = (char*) R_alloc(p[0],sizeof(char));
	memset (visited, 0, p[0]);
	tmp = root = (lnode*) R_alloc(1,sizeof(lnode));
	addroot = (lnode*) R_alloc(1,sizeof(lnode));
	PROTECT(ret = allocVector(VECSXP,1));
	cnt = 0;startpoint = 0;
	init_root_node(tmp,startpoint);

	while(1){
		get_groups(nodes,visited,startpoint,tmp);
		cnt += tmp->cnt+1;// we have to count the root (+1)
		startpoint = -1;
		for (i = 0;i<p[0];i++){ if (visited[i]==0){ startpoint = i;break;}}
		if (startpoint == -1) break;
		tmp = add_root_node(root,startpoint);
	}

	grpups = root->cnt_root+1;
	SET_VECTOR_ELT(ret, 0,  allocVector(VECSXP, grpups));
	grps = VECTOR_ELT(ret, 0);
	tmp = root;
	j=0;
	while (tmp != NULL){
		n = tmp;
		SET_VECTOR_ELT(grps, j,  allocVector(INTSXP, n->cnt+1));
		var_i = INTEGER(VECTOR_ELT(grps, j));
		k = 0;
		while (n != NULL){
			var_i[k] = n->itag+1;
			n = n->next;
			k++;
		}
		tmp = tmp->down;
		j++;
	}
	UNPROTECT(2);
	return(ret);
}

SEXP NMI_intersect_l(SEXP list_set,SEXP param){
	// param[0] total of points, param[1] repetitions
	// set matrix param[0]xparam[1] 1 indicates point i of repetition j (set[i][j]) present
	// labels matrix param[0]xparam[1] indicates labels of point i of repetition j (labels[i][j])
	//
	int *p,*idx_set,idx_set_len;
	double *cmatrix;
	int i,j,k,z,ii,ii2,jj;
	
	double *probs_1,*probs_2,*probs_12;
	double H_1,H_2,H_12,lg,lg12;
	int tMAX_label,MAX_1,MAX_2,A;
	char *class_translator,*class_translator_2;
	int *class_roseta,*class_roseta_2;
	SEXP cons_vector;
	
	PROTECT(param = coerceVector(param, INTSXP));

	
	p = INTEGER(param);
	PROTECT(cons_vector = allocVector(REALSXP, p[1]));
	cmatrix = REAL(cons_vector);
	
	class_translator = (char*) R_alloc(p[0],sizeof(char));
	class_translator_2 = (char*) R_alloc(p[0],sizeof(char));
	
	class_roseta = (int*) R_alloc(p[0],sizeof(int));
	class_roseta_2 = (int*) R_alloc(p[0],sizeof(int));
	
	tMAX_label=0;
	
	for(i=0;i<p[1];i++){
			// obtener tMAX_label en esta etapa
		memset (class_translator, 0, p[0]);	memset (class_translator_2, 0, p[0]);
		idx_set_len = length(VECTOR_ELT(list_set, i))/2;
		idx_set = INTEGER(VECTOR_ELT(list_set, i));
		
// 		Rprintf("idx_set_len: %i\n",idx_set_len);
		ii = 0;
		ii2 = idx_set_len;
		j = 0; jj = 0;
		for(k=0;k<idx_set_len;k++){
			if (class_translator[idx_set[k+ii]] == 0){
				class_translator[idx_set[k+ii]] = 1;
				class_roseta[j] = idx_set[k+ii];
				j++;
			}
			if (class_translator_2[idx_set[k+ii2]] == 0){
				class_translator_2[idx_set[k+ii2]] = 1;
				class_roseta_2[jj] = idx_set[k+ii2];
				jj++;
			}
		}
		j++;
		tMAX_label = (tMAX_label<j)?j:tMAX_label;
		jj++;
		tMAX_label = (tMAX_label<jj)?jj:tMAX_label;
		for(k=0;k<idx_set_len;k++){
			for(z=0;z<j;z++) idx_set[k+ii] = (idx_set[k+ii]==class_roseta[z])?(z+1):idx_set[k+ii];
		}
		for(k=0;k<idx_set_len;k++){
			for(z=0;z<jj;z++) idx_set[k+ii2] = (idx_set[k+ii2]==class_roseta_2[z])?(z+1):idx_set[k+ii2];
		}
	}
// 	Rprintf("tMAX_label: %i\n",tMAX_label);
// 	Rprintf("*\n");
	tMAX_label++;
	probs_1 = (double*) R_alloc(tMAX_label,sizeof(double));
	probs_2 = (double*) R_alloc(tMAX_label,sizeof(double));  
	probs_12 = (double*) R_alloc((tMAX_label*tMAX_label),sizeof(double));
	tMAX_label--;
// 	Rprintf("**\n");
	for(i=0;i<p[1];i++){
		idx_set_len = length(VECTOR_ELT(list_set, i))/2;
		idx_set = INTEGER(VECTOR_ELT(list_set, i));
		
		ii = 0;
		ii2 = idx_set_len;
// 		Rprintf("\n");
// 		for(k=0;k<idx_set_len;k++) if (idx_set[k+ii] > tMAX_label || idx_set[k+ii2]>tMAX_label || idx_set[k+ii2]<1 || idx_set[k+ii]<1) Rprintf("ERROR idx_set[k+ii2]: %i %i MAX: %i on %i\n",idx_set[k+ii],idx_set[k+ii2],tMAX_label,k);
		
		A = tMAX_label*tMAX_label;
		for(k=0;k<=A;k++) probs_12[k] = 0.;
		for(k=0;k<=tMAX_label;k++) probs_2[k] = probs_1[k] = 0.;
		MAX_2 = 0;MAX_1 = 0;
		for(k=0;k<idx_set_len;k++) {
			probs_1[idx_set[k+ii]]++;MAX_1 = (MAX_1<idx_set[k+ii]) ? idx_set[k+ii]:MAX_1;
			probs_2[idx_set[k+ii2]]++;MAX_2 = (MAX_2<idx_set[k+ii2]) ? idx_set[k+ii2]:MAX_2;
		}
		for(k=0;k<idx_set_len;k++) probs_12[((idx_set[k+ii2]-1)*MAX_1) + idx_set[k+ii]]++;
		H_1=0;for(k=1;k<=MAX_1;k++) {probs_1[k] /= p[0];H_1 -= (probs_1[k]>0)?probs_1[k]*log2(probs_1[k]):0.;}
		H_2=0;for(k=1;k<=MAX_2;k++) {probs_2[k] /= p[0];H_2 -= (probs_2[k]>0)?probs_2[k]*log2(probs_2[k]):0.;}
// 		Rprintf("\n");
		A = MAX_1* MAX_2;
		H_12 = 0;for(k=1;k<=A;k++) {probs_12[k] /= p[0];/*H_12 -= (probs_12[k]>0)?probs_12[k]*log2(probs_12[k]):0;*/}
		
		cmatrix[i] = 0.;
		for(k=1;k<=MAX_2;k++) {
			lg = (probs_2[k]>0)?log2(probs_2[k]):0;
			for(z=1;z<=MAX_1;z++){
				lg12 = (probs_12[z+(k-1)*MAX_1]>0) ? log2(probs_12[z+(k-1)*MAX_1]):0;
				cmatrix[i] += (probs_1[z]>0)?probs_12[z+(k-1)*MAX_1]* (lg12 -lg - log2(probs_1[z])):probs_12[z+(k-1)*MAX_1]*(lg12 -lg);
			}
		}
// 		Rprintf("%lf  %lf  %lf  %lf -- %i\n",H_1, H_2, H_12, cmatrix[i],i);
// 		cmatrix[i] = H_1 + H_2 + 2*H_12 - 4*cmatrix[i];
// 		cmatrix[i] = H_1 + H_2 - 2*cmatrix[i];
		cmatrix[i] =  2*cmatrix[i]/(H_1 + H_2);
	}
	UNPROTECT(2);
	return(cons_vector);
}

SEXP NMI_intersect_l_OMP(SEXP list_set,SEXP param){
	// param[0] total of points, param[1] repetitions
	// set matrix param[0]xparam[1] 1 indicates point i of repetition j (set[i][j]) present
	// labels matrix param[0]xparam[1] indicates labels of point i of repetition j (labels[i][j])
	//
	int *p,*idx_set,idx_set_len;
	double *cmatrix;
	int i,j,k,z,ii,ii2,jj;
	
	double **probs_1,**probs_2,**probs_12;
	
	double H_1,H_2,H_12,lg,lg12;
	int tMAX_label,MAX_1,MAX_2,A;
	char **class_translator,**class_translator_2;
	int **class_roseta,**class_roseta_2;
	int MAX_T,Tnum,*tMAX_label_tmp;
	SEXP cons_vector;
	
	PROTECT(param = coerceVector(param, INTSXP));
	p = INTEGER(param);
	PROTECT(cons_vector = allocVector(REALSXP, p[1]));
	cmatrix = REAL(cons_vector);
	
	MAX_T = omp_get_max_threads();
	
	class_translator = (char**) R_alloc(MAX_T,sizeof(char*));
	for(i=0;i<MAX_T;i++) class_translator[i] = (char*) R_alloc(p[0],sizeof(char));
	
	class_translator_2 = (char**) R_alloc(MAX_T,sizeof(char*));
	for(i=0;i<MAX_T;i++) class_translator_2[i] = (char*) R_alloc(p[0],sizeof(char));
	
	class_roseta = (int**) R_alloc(MAX_T,sizeof(int*));
	for(i=0;i<MAX_T;i++) class_roseta[i] = (int*) R_alloc(p[0],sizeof(int));
	
	class_roseta_2 = (int**) R_alloc(MAX_T,sizeof(int*));
	for(i=0;i<MAX_T;i++) class_roseta_2[i] = (int*) R_alloc(p[0],sizeof(int));
	
	tMAX_label_tmp = (int*) R_alloc(MAX_T,sizeof(int));
	
	tMAX_label=0;
# pragma omp parallel private(Tnum,i,idx_set_len,idx_set,ii,ii2,j,jj,k,z)
{
  Tnum = omp_get_thread_num();
  tMAX_label_tmp[Tnum] = 0;
# pragma omp for schedule(static)
	for(i=0;i<p[1];i++){
			// obtener tMAX_label en esta etapa
		memset (class_translator[Tnum], 0, p[0]);memset (class_translator_2[Tnum], 0, p[0]);
		idx_set_len = length(VECTOR_ELT(list_set, i))/2;
		idx_set = INTEGER(VECTOR_ELT(list_set, i));
		
// 		Rprintf("idx_set_len: %i\n",idx_set_len);
		ii = 0;
		ii2 = idx_set_len;
		j = 0; jj = 0;
		for(k=0;k<idx_set_len;k++){
			if (class_translator[Tnum][idx_set[k+ii]] == 0){
				class_translator[Tnum][idx_set[k+ii]] = 1;
				class_roseta[Tnum][j] = idx_set[k+ii];
				j++;
			}
			if (class_translator_2[Tnum][idx_set[k+ii2]] == 0){
				class_translator_2[Tnum][idx_set[k+ii2]] = 1;
				class_roseta_2[Tnum][jj] = idx_set[k+ii2];
				jj++;
			}
		}
		j++;
		tMAX_label_tmp[Tnum] = (tMAX_label_tmp[Tnum]<j)?j:tMAX_label_tmp[Tnum];
		jj++;
		tMAX_label_tmp[Tnum] = (tMAX_label_tmp[Tnum]<jj)?jj:tMAX_label_tmp[Tnum];
		for(k=0;k<idx_set_len;k++){
			for(z=0;z<j;z++) idx_set[k+ii] = (idx_set[k+ii]==class_roseta[Tnum][z])?(z+1):idx_set[k+ii];
		}
		for(k=0;k<idx_set_len;k++){
			for(z=0;z<jj;z++) idx_set[k+ii2] = (idx_set[k+ii2]==class_roseta_2[Tnum][z])?(z+1):idx_set[k+ii2];
		}
	}
}
	for(i=0;i<MAX_T;i++) tMAX_label = (tMAX_label<tMAX_label_tmp[i])? tMAX_label_tmp[i] : tMAX_label;
	Rprintf("tMAX_label: %i\n",tMAX_label);
	Rprintf("*\n");
	tMAX_label++;
	probs_1 = (double**) R_alloc(MAX_T,sizeof(double*));
	for(i=0;i<MAX_T;i++) probs_1[i] = (double*) R_alloc(tMAX_label,sizeof(double));
	
	probs_2 = (double**) R_alloc(MAX_T,sizeof(double*));
	for(i=0;i<MAX_T;i++) probs_2[i] = (double*) R_alloc(tMAX_label,sizeof(double));
	
	probs_12 = (double**) R_alloc(MAX_T,sizeof(double*));
	for(i=0;i<MAX_T;i++) probs_12[i] = (double*) R_alloc( (tMAX_label*tMAX_label) ,sizeof(double));
	
	tMAX_label--;
// 	Rprintf("**\n");
# pragma omp parallel private(Tnum,i,idx_set_len,idx_set,ii,ii2,A,k,MAX_1,MAX_2,H_1,H_2,H_12,z,lg,lg12)
{
  Tnum = omp_get_thread_num();
# pragma omp for schedule(static)
	for(i=0;i<p[1];i++){
		idx_set_len = length(VECTOR_ELT(list_set, i))/2;
		idx_set = INTEGER(VECTOR_ELT(list_set, i));
		
		ii = 0;
		ii2 = idx_set_len;
// 		Rprintf("\n");
// 		for(k=0;k<idx_set_len;k++) if (idx_set[k+ii] > tMAX_label || idx_set[k+ii2]>tMAX_label || idx_set[k+ii2]<1 || idx_set[k+ii]<1) Rprintf("ERROR idx_set[k+ii2]: %i %i MAX: %i on %i\n",idx_set[k+ii],idx_set[k+ii2],tMAX_label,k);
		
		A = tMAX_label*tMAX_label;
		for(k=0;k<=A;k++) probs_12[Tnum][k] = 0.;
		for(k=0;k<=tMAX_label;k++) probs_2[Tnum][k] = probs_1[Tnum][k] = 0.;
		MAX_2 = 0;MAX_1 = 0;
		for(k=0;k<idx_set_len;k++) {
			probs_1[Tnum][idx_set[k+ii]]++;MAX_1 = (MAX_1<idx_set[k+ii]) ? idx_set[k+ii]:MAX_1;
			probs_2[Tnum][idx_set[k+ii2]]++;MAX_2 = (MAX_2<idx_set[k+ii2]) ? idx_set[k+ii2]:MAX_2;
		}
		for(k=0;k<idx_set_len;k++) probs_12[Tnum][((idx_set[k+ii2]-1)*MAX_1) + idx_set[k+ii]]++;
		H_1=0;for(k=1;k<=MAX_1;k++) {probs_1[Tnum][k] /= p[0];H_1 -= (probs_1[Tnum][k]>0)?probs_1[Tnum][k]*log2(probs_1[Tnum][k]):0.;}
		H_2=0;for(k=1;k<=MAX_2;k++) {probs_2[Tnum][k] /= p[0];H_2 -= (probs_2[Tnum][k]>0)?probs_2[Tnum][k]*log2(probs_2[Tnum][k]):0.;}
// 		Rprintf("\n");
		A = MAX_1* MAX_2;
		H_12 = 0;for(k=1;k<=A;k++) {probs_12[Tnum][k] /= p[0];/*H_12 -= (probs_12[k]>0)?probs_12[k]*log2(probs_12[k]):0;*/}
		
		cmatrix[i] = 0.;
		for(k=1;k<=MAX_2;k++) {
			lg = (probs_2[Tnum][k]>0)?log2(probs_2[Tnum][k]):0;
			for(z=1;z<=MAX_1;z++){
				lg12 = (probs_12[Tnum][z+(k-1)*MAX_1]>0) ? log2(probs_12[Tnum][z+(k-1)*MAX_1]):0;
				cmatrix[i] += (probs_1[Tnum][z]>0)?probs_12[Tnum][z+(k-1)*MAX_1]* (lg12 -lg - log2(probs_1[Tnum][z])):probs_12[Tnum][z+(k-1)*MAX_1]*(lg12 -lg);
			}
		}
// 		Rprintf("%lf  %lf  %lf  %lf -- %i\n",H_1, H_2, H_12, cmatrix[i],i);
// 		cmatrix[i] = H_1 + H_2 + 2*H_12 - 4*cmatrix[i];
// 		cmatrix[i] = H_1 + H_2 - 2*cmatrix[i];
		cmatrix[i] =  2*cmatrix[i]/(H_1 + H_2);
	}
}

	UNPROTECT(2);
	return(cons_vector);
}

SEXP tagger (SEXP Data, SEXP param){
	// param[0] data length
	register int i,j,k;
	int *r,*n_from,*_m,len,len_from;
	int class_val;
	SEXP ret,nodes,nodes_i_from;
	
	PROTECT(param = coerceVector(param, INTSXP));
	_m = INTEGER(param); 
	PROTECT(ret = allocVector(INTSXP, _m[0]));
	r = INTEGER(ret);

	nodes = VECTOR_ELT(Data, 0);
	len = length(nodes);
// # pragma omp parallel private(Tnum)
// 	{	
// 	Tnum = omp_get_thread_num();
// 	# pragma omp for private(i,j,to,val,flag,k) schedule(dynamic,8)
	for(i = 0; i < _m[0]; i++) r[i] = 0;
	class_val = 1;
	for(i = 0; i < len; i++){
		nodes_i_from = VECTOR_ELT(nodes, i);
		n_from = INTEGER(nodes_i_from);
		len_from = length(VECTOR_ELT(nodes, i));
		for(j = 0; j < len_from; j++){
			if (n_from[j] <= 0 && r[n_from[j]-1]!=0) Rprintf("Warning: malformation on list position %i point %i\n",j,(n_from[j]-1));
			r[n_from[j]-1] = class_val;
		}
		class_val++;
	}
	class_val--;
// 	Rprintf("\n\n");
//	}	
	UNPROTECT(2);
	return(ret);
}

SEXP outlayer_merger (SEXP list, SEXP Data, SEXP par){
	/* par[0] len of data, par[1] dim of data ,p[2]conn type,p[3] 0 no pen, 1 exponencial, 2 lineal, 3 power, 4 exponencial (exp(x-1)) ,5 power (x^n)-1
	1 pen,p[4] metric,p[5] lineal constant*/
	
	
	//min_val changed for euclidean
	
	lnode *root,*tmp,*min_root,*n,*m,*addroot,*tmp_roota,*tmp_rootb;
	SEXP nodes,ret,grps;
	
	int cnt,startpoint,grpups,e;
	int *p,*var_i;
	int i,k,j;
	int from,to,min_from,min_to,bins,len_i,len_d,*_ns;
	int p1,p2;
// 	int I;
	double *data,dif,val,min_val,*var_ret,d_one=1.;
	
	char bool,bool1;
	char *visited;
	
	PROTECT(par = coerceVector(par, INTSXP));
	PROTECT(Data = coerceVector(Data, REALSXP));

	data = REAL(Data);
	p = INTEGER(par);
	e = p[5];
// 	Rprintf("%i \n",p[5]);
	
	nodes = VECTOR_ELT(list, 1);
	visited = (char*) R_alloc(p[0],sizeof(char));
	memset (visited, 0, p[0]);
	tmp = root = (lnode*) R_alloc(1,sizeof(lnode));
	addroot = (lnode*) R_alloc(1,sizeof(lnode));
	PROTECT(ret = allocVector(VECSXP,2));
// 	Rprintf("%i %i %i %i %i %lf\n",p[0],p[1],p[2],p[3],p[4],_mean[0]);

	cnt = 0;startpoint = 0;
	init_root_node(tmp,startpoint);
	while(1){
		get_groups(nodes,visited,startpoint,tmp);
		cnt += tmp->cnt+1;// we have to count the root (+1)
		startpoint = -1;
		for (i = 0;i<p[0];i++){ if (visited[i]==0){ startpoint = i;break;}}
		if (startpoint == -1) break;
		tmp = add_root_node(root,startpoint);
// 			Rprintf("cnt: %i\n",cnt);
	}
// 	Rprintf("*\n");
	grpups = root->cnt_root+1;
	SET_VECTOR_ELT(ret, 0,  allocVector(VECSXP, grpups));
	grps = VECTOR_ELT(ret, 0);
	tmp = root;
	j=0;
	while (tmp != NULL){
		n = tmp;
		SET_VECTOR_ELT(grps, j,  allocVector(INTSXP, n->cnt+1));
		var_i = INTEGER(VECTOR_ELT(grps, j));
		k = 0;
		while (n != NULL){
//  			Rprintf("%i ",n->itag);
			var_i[k] = n->itag+1;
			n = n->next;
			k++;
		}
//  		Rprintf("\n\n");
		tmp = tmp->down;
		j++;
	}
	/* one to the rest */
	if (grpups > 1){
		init_root_node(addroot,0);
// 		I = 0;
		do{
// 			Rprintf("I: %i\n",I);
// 			I++;
			tmp_roota = root;
			tmp_rootb = root;
		
			while (tmp_roota->down != NULL && !(tmp_roota->cnt <p[6])) {tmp_roota = tmp_roota->down;}
			if (tmp_roota->down == NULL && !(tmp_roota->cnt <p[6])) {tmp_roota = NULL;tmp_rootb = NULL;}
			else tmp_roota = (tmp_rootb==tmp_roota)?tmp_roota->down:tmp_roota;
			m = tmp_roota;
			
			min_val = DBL_MAX;
			bool1=0;
			while (tmp_rootb != NULL){
			bool1=1;
				while (m != NULL){
				n = tmp_rootb;
				from = (m->itag);
// 				Rprintf("from: %i\n",from);
					while (n != NULL){
						to = (n->itag);
// 							Rprintf("	to: %i\n",to);
						val = 0;
						val =  data[from+to*p[0]];
						bool = (min_val > val);
						min_root = (bool)?tmp_rootb:min_root;
						min_from = (bool)?from:min_from;
						min_to = (bool)?to:min_to;
						min_val = (bool)?val:min_val;
						n = n->next;
					}
				m = m->next;
				}
				tmp_rootb = tmp_rootb->down;
				m = tmp_roota;
				tmp_rootb = (tmp_rootb==tmp_roota)?tmp_roota->down:tmp_rootb;
			}
			if(bool1){
				add_ij_node(addroot,min_from,min_to,min_val);
				merge(tmp_roota,min_root);
// 				Rprintf("addroot->cnt %i\n",addroot->cnt);
			}
		}while(bool1);
	}
	
//  	tmp_roota = root;
//  	while (tmp_roota->down != NULL) {Rprintf("tmp_roota->cnt: %i\n",tmp_roota->cnt);tmp_roota = tmp_roota->down;}
//  	Rprintf("tmp_roota->cnt: %i\n",tmp_roota->cnt);
// 	Rprintf("grpups: %i	addroot->cnt %i\n",grpups,addroot->cnt);

	if (grpups > 1 && addroot->cnt > 0){
		SET_VECTOR_ELT(ret, 1,  allocMatrix(REALSXP,addroot->cnt,3));
		var_ret = REAL(VECTOR_ELT(ret, 1));
// 		Rprintf("2\n");
		j=0;
		tmp = addroot->next;
		while (tmp != NULL){
			var_ret[j] = tmp->itag+1;
// 			Rprintf("tmp->itag+1:  %i\n",(tmp->itag+1));
			var_ret[j+addroot->cnt] = tmp->jtag+1;
// 			Rprintf("tmp->jtag+1:  %i\n",(tmp->jtag+1));
			var_ret[j+2*addroot->cnt] = tmp->dist;
// 			Rprintf("tmp->dist:  %lf\n",(tmp->dist));
			tmp = tmp->next;
			j++;
		}
// 		Rprintf("3 j: %i\n",j);
	}
	else{
		SET_VECTOR_ELT(ret, 1,  allocMatrix(REALSXP,1,1));
		var_ret = REAL(VECTOR_ELT(ret, 1));
		var_ret[0] = 0;
	}
// 	Rprintf("4\n");
	UNPROTECT(3);
	return(ret);
}

void matrix_edge (SEXP list, SEXP DistMat, SEXP par){
  /* par[0] len of list , and size of DistMat, DistMat is a square matrix
   * 
   * SEXP list type list()
   * 
  */

	SEXP nodes,*nodes_vec;

	double *data;
	int *p,*n;
	int len,i,j,Tnum,MAX_T;
	int one = 1;
	PROTECT(par = coerceVector(par, INTSXP));
	PROTECT(DistMat = coerceVector(DistMat, REALSXP));

	data = REAL(DistMat);
	p = INTEGER(par);

	nodes = VECTOR_ELT(list, 1);
	MAX_T = omp_get_max_threads();
	nodes_vec = (SEXP*) R_alloc(MAX_T,sizeof(SEXP));
# pragma omp parallel private(Tnum,i,j,n,len)
{
	Tnum = omp_get_thread_num();
# pragma omp for schedule(dynamic,32)
	for(i = 0; i < p[0]; i++){
	  nodes_vec[Tnum] = VECTOR_ELT(nodes, i);
	  n = INTEGER(nodes_vec[Tnum]);
	  len = length(nodes_vec[Tnum]);
	  for(j = 0; j < len; j++){
	    /* n stores the nodes using R addresing convention start at 1 and end at n*/
	    data[i*p[0]+ (n[j]-1)] += one;
	  }
	}
	#pragma omp barrier
# pragma omp for schedule(dynamic,32)
	for(i = 0; i < p[0]; i++){
	  nodes_vec[Tnum] = VECTOR_ELT(nodes, i);
	  n = INTEGER(nodes_vec[Tnum]);
	  len = length(nodes_vec[Tnum]);
	  for(j = 0; j < len; j++){
	    /* n stores the nodes using R addresing convention start at 1 and end at n*/
	    data[i+ (n[j]-1)*p[0]] += one;
	  }
	}

}
	UNPROTECT(2);

}

void matrix_div (SEXP DistMat, SEXP D_val, SEXP par){
  /* par[0] len of list , and size of DistMat, DistMat is a square matrix
   * 
   * SEXP list type list()
   * 
  */

	double *data,*d;
	int *p;
	int i;

	PROTECT(par = coerceVector(par, INTSXP));
	PROTECT(DistMat = coerceVector(DistMat, REALSXP));
	PROTECT(D_val = coerceVector(D_val, REALSXP));
	
	data = REAL(DistMat);
	d = REAL(D_val);
	p = INTEGER(par);

# pragma omp parallel private(Tnum,i,j,n,len)
{
# pragma omp for schedule(dynamic,128)
	for(i = 0; i < p[0]*p[0]; i++) data[i] /= d[0];

}
	UNPROTECT(3);
}