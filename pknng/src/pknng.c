#include "def.h"


SEXP all_dijkstra (SEXP ,SEXP );
inline void dijkstra(SEXP , int , int , int* , double *,int *,struct fibheap * ,struct fibheap_el ** );
inline int which_is_min(double *, int );
inline int which_is_max(double *, int );
inline double intpow(double,int);

SEXP get_k_neighbors (SEXP Data, SEXP param);
SEXP make_symmetric (SEXP Data, SEXP param, SEXP thr);
void get_groups(SEXP list,char* visited,int startpoint,lnode* root);
SEXP connect_groups_oneToTheRest (SEXP list, SEXP Data, SEXP par,SEXP mean);
SEXP connect_groups_all (SEXP list, SEXP Data, SEXP par,SEXP mean);




SEXP all_dijkstra (SEXP list,SEXP par){
	
	SEXP ret;
	int *p,*idx,i,*previous;
	double *dist;
	struct fibheap_el **x;
	struct fibheap *l;
	
	PROTECT(par = coerceVector(par, INTSXP));
	p = INTEGER(par);
	
	PROTECT(ret = allocMatrix(REALSXP, p[0],p[0]));
	dist = REAL(ret);
	previous = (int*) R_alloc(p[0],sizeof(int));
	idx = (int*) R_alloc(p[0],sizeof(int));
	x = (struct fibheap_el **) R_alloc(p[0],sizeof (struct fibheap_el*));
	l = (struct fibheap *) R_alloc(1,sizeof (struct fibheap));
	
	l = fh_makekeyheap();
	for(i=0;i<p[0];i++){
		x[i] = (struct fibheap_el *) R_alloc(1,sizeof (struct fibheap_el));
		idx[i]=i;
	}
	for (i = 0;i<p[0];i++) {
		dijkstra(list,p[0],i,previous,&dist[i*p[0]],idx,l,x);
		fh_clearyheap(l);
	}
	
	UNPROTECT(2);
	return(ret);
}

inline void dijkstra(SEXP list, int par, int source, int* prev, double *d,int *data,struct fibheap * l,struct fibheap_el ** x){
	
	int i,len_nodes,len,n,cnt;
	double *adj_w,alt;
	int *adj;
	SEXP nodes,weights;
	struct fibheap_el *node;
 
	nodes = VECTOR_ELT(list, 1);
	weights = VECTOR_ELT(list, 0);
	len_nodes = par;
	for (i=0;i < par;i++){
		d[i] = DBL_MAX;prev[i] = -1;
		fh_insertkey_el(l,x[i], d[i],(void *)&data[i]);
	}
	d[source] = 0;
	fh_replacekeydata(l,x[source],d[source], (void *)&data[source]);
	
	node = l->fh_min;
	n = ((int*)node->fhe_data)[0];
	cnt = 0;
	while (1){
		fh_extractmin(l);
		len = length(VECTOR_ELT(nodes, n));
		adj = INTEGER(VECTOR_ELT(nodes, n));
		adj_w = REAL(VECTOR_ELT(weights, n));
		for (i=0;i < len;i++){
			alt = d[n]+adj_w[i];
			if (alt < d[adj[i]-1]){
				d[adj[i]-1] = alt;
				prev[(adj[i])-1] = n;
				fh_replacekeydata(l,x[adj[i]-1],d[adj[i]-1], (void *)&data[adj[i]-1]);
			}
		}
		node = l->fh_min;
		if (l->fh_min==NULL) break;
		else{
			n = ((int*)node->fhe_data)[0];
			cnt++;
		}
	}
}

inline int which_is_min(double a[], int l){
  
	int min_index = 0;
	double min_a = a[0];  
	char bool;
	int i;
  
	for(i = 0; i < l ; i++){
		bool = (a[i] < min_a);
		min_index = (bool)?i:min_index;
		min_a = (bool)?a[i]:min_a;
    
	}
	return(min_index);  
}


inline int which_is_max(double a[], int l){
  
	int max_index = 0;
	double max_a = a[0];  
	char bool;
	int i;
  
	for(i = 0; i < l ; i++){
		bool = (a[i] > max_a);
		max_index = (bool)?i:max_index;
		max_a = (bool)?a[i]:max_a;
	}
	return(max_index);  
}


inline double intpow(double b,int e){
	int i;
	double r;
	r = (e==0)?1.:b;
	for (i=1;i<e;i++) r *=b;
	return r;
}

SEXP get_k_neighbors (SEXP Data, SEXP param) 
// Data is a mxm matrix of distancies or similarities
//param = m cant de elementos de Data, dim dimensiones de Data (not used left for compatibility), k vecinos
{
	register int i, j, k, nb, ii;
	int *xb,histoparam[4];
	double *xa,*xab;
	
	int *index,*idx,*mchunk_i;
	int max_index,min_index;
  
	double *values,*values_2;
	double max_value;
	int len;
	int Tnum,MAX_T,dummy;
	char bool_var,bool_equal,bool_equal_o;
	int *rep_neighbours,*len_rep_neigh;
	
	SEXP ab,idxab,_ab,_idxab;
	SEXP ret;
	
	
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(param = coerceVector(param, INTSXP));
	
	nb = length(param); 
	xa = REAL(Data); xb = INTEGER(param);
	
	PROTECT(ret = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret, 0,  allocVector(VECSXP,xb[0]));
	SET_VECTOR_ELT(ret, 1,  allocVector(VECSXP,xb[0]));
	
	ab = VECTOR_ELT(ret, 0);
	idxab = VECTOR_ELT(ret, 1);
	
	values = (double*) R_alloc((xb[2]),sizeof(double));
	values_2 = (double*) R_alloc((xb[2]),sizeof(double));
	index = (int*) R_alloc((xb[2]),sizeof(int));
	len_rep_neigh = (int*) R_alloc((xb[2]),sizeof(int));
	rep_neighbours = (int*) R_alloc((xb[2]*xb[0]),sizeof(int));
	
	for(j = 0; j < xb[2]; j++) len_rep_neigh[j] = 0;
	for (i = 0;i<xb[0];i++){
		for(j = 0; j < xb[2]; j++){
			len_rep_neigh[j] = 0;
			index[j] = -1;
			values[j] = DBL_MAX; 
		}
		ii = i*xb[0];
		len = 0;
		max_index = 0;
		for(j = 0 ; j < xb[0]; j++){
			/* Replace Maximal Value in List*/
			bool_equal_o = 0;
			for(k = 0; k < xb[2]; k++) {
				bool_equal = (values[k] == xa[j+ii]) && j!=i;
				bool_equal_o = bool_equal_o || bool_equal;
				rep_neighbours[ len_rep_neigh[k] + k*xb[0]] = (bool_equal) ? j:rep_neighbours[ len_rep_neigh[k] + k*xb[0]];
				len_rep_neigh[k] += (bool_equal) ? 1 : 0; 
				k = (bool_equal)? xb[2]:k;
			}
			bool_equal = bool_equal^0x01;
			bool_var = bool_equal && (xa[j+ii] < values[max_index]) && j!=i;
			values[max_index] = (bool_var) ? xa[j+ii]:values[max_index];
			index [max_index] = (bool_var) ? j : index [max_index];
			len_rep_neigh[max_index] = (bool_var) ? 0 : len_rep_neigh[max_index]; 
			/* Set new Maximum */
			
			max_index = (bool_var) ? which_is_max(values,xb[2]):max_index;
		}
		for(j = 0; j < xb[2]+1; j++) values_2[j] = values[j]; 
		/* Mark all "non-neighbors" with -1 */
		while(len<xb[2]) {
			min_index = which_is_min(values_2, xb[2]);
			len += len_rep_neigh[min_index] + 1; // + 1 por el elemento en values[] ... primer minimo encontrado
			values_2[min_index] = DBL_MAX;
		}
		SET_VECTOR_ELT(ab, i,  allocVector(REALSXP,len));
		SET_VECTOR_ELT(idxab, i,  allocVector(INTSXP,len));
		_ab = VECTOR_ELT(ab, i);
		_idxab = VECTOR_ELT(idxab, i);
		xab = REAL(_ab);idx = INTEGER(_idxab);
		k=0;
		while(k<len){
			min_index = which_is_min(values, xb[2]);
			xab[k] = values[min_index];
			idx[k] = index[min_index];
			k++;
			for(j=0;j<len_rep_neigh[min_index];j++){
				xab[k] = values[min_index];
				idx[k] = rep_neighbours[ j + min_index*xb[0]];
				k++;
			}
			values[min_index] = DBL_MAX;
		}
	}/* Next row */
	UNPROTECT(3);
	return(ret);
}

SEXP make_symmetric (SEXP Data, SEXP param, SEXP thr){
	// param[0] data length, param[1] k vecinos (no meaning in here cos we re-writed the thing but we kept it for compatibility)
	register int i,j,k;
	double tmp;
	double *d_to,*d_from;
	int *_m,*n_to,*n_from,*m,len_to,len_from;
	double *sg;
	int to,flag,len,*array_off;
	double val;
	SEXP dist,nodes,ret,ret_dist,ret_nodes,array_dist,array_nodes,dist_i_from,nodes_i_from,dist_i_to,nodes_i_to;
	double *var_d;
	int *var_i;
	double *var_int,media;
	int MAX_T,Tnum;
	lnode *addon,*node;
	unsigned long int TOT_conns;
	
	PROTECT(thr = coerceVector(thr, REALSXP));
	PROTECT(param = coerceVector(param, INTSXP));
	
	
	_m = INTEGER(param);sg = REAL(thr);
	PROTECT(ret = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ret, 0,  allocVector(VECSXP, _m[0]));
	SET_VECTOR_ELT(ret, 1,  allocVector(VECSXP, _m[0]));
	ret_dist = VECTOR_ELT(ret, 0);
	ret_nodes = VECTOR_ELT(ret, 1);
	
	array_off = (int*) R_alloc(_m[0],sizeof(int));
	addon = (lnode*) R_alloc(_m[0],sizeof(lnode));
	dist = VECTOR_ELT(Data, 0);
	nodes = VECTOR_ELT(Data, 1);
	
	for(i = 0; i < _m[0]; i++) {init_root_node(&addon[i], i);array_off[i] = 0;}
	for(i = 0; i < _m[0]; i++){
		dist_i_from = VECTOR_ELT(dist, i); nodes_i_from = VECTOR_ELT(nodes, i);
		d_from = REAL(dist_i_from); n_from = INTEGER(nodes_i_from);
		len_from = length(dist_i_from);
		for(j = 0; j < len_from; j++){
			to = n_from[j];
			val = d_from[j];
			dist_i_to = VECTOR_ELT(dist, to); nodes_i_to = VECTOR_ELT(nodes, to);
			d_to = REAL(dist_i_to); n_to = INTEGER(nodes_i_to);
			len_to = length(dist_i_to);
			flag = 0;
			for(k = 0; k < len_to; k++) flag = (n_to[k] == i) ? 1:flag;
			if (flag==0){ 
				if (val <= sg[0]) 
					add_node(&addon[to], i,val); 					//**************
				else{ n_from[j] = -1;d_from[j] = -1.;array_off[i]--;} 
			}
		}
	}
	media = 0.;/* ver */
	for(i = 0; i < _m[0]; i++){
		dist_i_from = VECTOR_ELT(dist, i); nodes_i_from = VECTOR_ELT(nodes, i);
		d_from = REAL(dist_i_from); n_from = INTEGER(nodes_i_from);
		_m[1] = length(dist_i_from);
		len = addon[i].cnt+_m[1]+array_off[i];
		TOT_conns += len;
		SET_VECTOR_ELT(ret_dist,i,  allocVector(REALSXP,len));
		SET_VECTOR_ELT(ret_nodes,i ,allocVector(INTSXP,len));
		array_dist = VECTOR_ELT(ret_dist, i);
		array_nodes = VECTOR_ELT(ret_nodes, i);
		var_d = REAL(array_dist);/* ver */
		var_i = INTEGER(array_nodes);/* ver */
		for (j = 0, k = 0;j<_m[1];j++){
			if (n_from[j] > -1){
				var_d[k] = d_from[j];
				var_i[k] = n_from[j]+1;
				media += var_d[k];/* ver */
				k++;
			}
		}
		if (addon[i].cnt>0) {
			j = k;
			node = &addon[i];/* ver */
			while(node->next != NULL){
				node=node->next;
				var_d[j] = node->dist; var_i[j] = (node->itag)+1;
				media += var_d[j];
				j++;
			}
		}
	}
	SET_VECTOR_ELT(ret, 2,  allocVector(REALSXP, 2));
	var_int = REAL(VECTOR_ELT(ret, 2));
	var_int[0] = TOT_conns;
	var_int[1] = media/TOT_conns;
	UNPROTECT(3);
	return(ret);
}

void get_groups(SEXP list,char* visited,int startpoint,lnode* root){
	
	SEXP nodes_vec;
	int len,i,idx;
	int *nodes;
// 	Rprintf("--- \n");
	nodes_vec = VECTOR_ELT(list, startpoint);
	len = length(nodes_vec);
// 	Rprintf("--- \n");
	nodes = INTEGER(nodes_vec);
	visited[startpoint] = 1;
// 	Rprintf("--- \n");
	for (i=0;i<len;i++) {
		idx = (int)nodes[i]-1;
		if (visited[idx]==0){
			//Rprintf("idx  %i \n",idx);
			add_node(root,idx,-1.);
			get_groups(list,visited,idx,root);
		}
	}
	return;
}

SEXP connect_groups (SEXP list, SEXP Data, SEXP par,SEXP mean){
	/* par[0] len of data, par[1] dim of data ,p[2]conn type,p[3] 0 no pen, 1 exponencial, 2 lineal, 3 power, 4 exponencial (exp(x-1)) ,5 power (x^n)-1
	1 pen,p[4] metric,p[5] lineal constant*/
	
	
	//min_val changed for euclidean
	
	int *p,*var_i;
	lnode *root,*tmp,*min_root,*n,*m,*addroot,*tmp_roota,*tmp_rootb;
	SEXP nodes,ret,grps;
	char *visited;
	int cnt,startpoint,grpups,e;
	int i,k,j;
	int from,to,min_from,min_to,bins,len_i,len_d;
	double *data,dif,val,min_val,*_mean,*var_ret,d_one=1.;
	char bool;
	PROTECT(par = coerceVector(par, INTSXP));
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(mean = coerceVector(mean, REALSXP));
	data = REAL(Data);
	p = INTEGER(par);
	_mean = REAL(mean);
	e = p[5];

	nodes = VECTOR_ELT(list, 1);
	visited = (char*) R_alloc(p[0],sizeof(char));
	memset (visited, 0, p[0]);
	tmp = root = (lnode*) R_alloc(1,sizeof(lnode));
	addroot = (lnode*) R_alloc(1,sizeof(lnode));
	PROTECT(ret = allocVector(VECSXP,2));
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
	/* one to the rest */
	if (grpups > 1){
		init_root_node(addroot,0);
		while (root->down != NULL){
			tmp_rootb = root->down; 
			m = root;
			min_val = DBL_MAX;
			while (tmp_rootb != NULL){
				while (m != NULL){
					n = tmp_rootb;
					from = (m->itag);
					while (n != NULL){
						to = (n->itag);
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
				m = root;
			}
			min_val = (p[3]==1)? min_val*exp(min_val/_mean[0]):min_val;
			min_val = (p[3]==2)? min_val*e:min_val;
			min_val = (p[3]==3)? min_val*intpow((min_val/_mean[0]),e):min_val;
			min_val = (p[3]==4)? min_val*exp((min_val/_mean[0])-d_one):min_val;
			min_val = (p[3]==5)? min_val*(intpow((min_val/_mean[0]),e)+d_one):min_val;
			add_ij_node(addroot,min_from,min_to,min_val);
			merge(root,min_root);
		}
	}
	if (grpups > 1){
		SET_VECTOR_ELT(ret, 1,  allocMatrix(REALSXP,addroot->cnt,3));
		var_ret = REAL(VECTOR_ELT(ret, 1));
		j=0;
		tmp = addroot->next;
		while (tmp != NULL){
			var_ret[j] = tmp->itag+1;
			var_ret[j+addroot->cnt] = tmp->jtag+1;
			var_ret[j+2*addroot->cnt] = tmp->dist;
			tmp = tmp->next;
			j++;
		}
	}
	else{
		SET_VECTOR_ELT(ret, 1,  allocMatrix(REALSXP,1,1));
		var_ret = REAL(VECTOR_ELT(ret, 1));
		var_ret[0] = 0;
	}
	UNPROTECT(4);
	return(ret);
}

SEXP connect_groups_oneToTheRest (SEXP list, SEXP Data, SEXP par,SEXP mean){
	/* par[0] len of data, par[1] dim of data ,p[2]conn type,p[3] 0 no pen, 1 exponencial, 2 lineal, 3 power,  4 exponencial (exp(x-1)) ,5 power (x^n)-1
	1 pen,p[4] metric,p[5] lineal constant*/
	
	
	//min_val changed for euclidean
	
	int *p,*var_i;
	lnode *root,*tmp,*min_root,*n,*m,*mm,*addroot,*tmp_roota,*tmp_rootb;
	SEXP nodes,ret,grps;
	char *visited;
	int cnt,startpoint,grpups,e;
	int i,k,j;
	int from,to,min_from,min_to,bins,len_i,len_d;
	double *data,dif,val,min_val,*_mean,*var_ret,d_one=1.;
	char bool;
	PROTECT(par = coerceVector(par, INTSXP));
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(mean = coerceVector(mean, REALSXP));
	data = REAL(Data);
	p = INTEGER(par);
	_mean = REAL(mean);
	e = p[5];
	
	nodes = VECTOR_ELT(list, 1);
	visited = (char*) R_alloc(p[0],sizeof(char));
	memset (visited, 0, p[0]);
	tmp = root = (lnode*) R_alloc(1,sizeof(lnode));
	addroot = (lnode*) R_alloc(1,sizeof(lnode));
	PROTECT(ret = allocVector(VECSXP,2));

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
	/* one to the rest */
	if (grpups > 1){
		init_root_node(addroot,0);
		mm=root;
		while (mm->down != NULL){
			tmp_rootb = mm->down; 
			m = mm;
			min_val = DBL_MAX;
			while (tmp_rootb != NULL){
				while (m != NULL){
					n = tmp_rootb;
					from = (m->itag);
					while (n != NULL){
						to = (n->itag);
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
				min_val = (p[3]==1)? min_val*exp(min_val/_mean[0]):min_val;
				min_val = (p[3]==2)? min_val*e:min_val;
				min_val = (p[3]==3)? min_val*intpow((min_val/_mean[0]),e):min_val;
				min_val = (p[3]==4)? min_val*exp((min_val/_mean[0])-d_one):min_val;
				min_val = (p[3]==5)? min_val*(intpow((min_val/_mean[0]),e)+d_one):min_val;
				add_ij_node(addroot,min_from,min_to,min_val);
				tmp_rootb = tmp_rootb->down;
				m = mm;
			}
			mm=mm->down;
		}
	}
	if (grpups > 1){
		SET_VECTOR_ELT(ret, 1,  allocMatrix(REALSXP,addroot->cnt,3));
		var_ret = REAL(VECTOR_ELT(ret, 1));
		j=0;
		tmp = addroot->next;
		while (tmp != NULL){
			var_ret[j] = tmp->itag+1;
			var_ret[j+addroot->cnt] = tmp->jtag+1;
			var_ret[j+2*addroot->cnt] = tmp->dist;
			tmp = tmp->next;
			j++;
		}
	}
	else{
		SET_VECTOR_ELT(ret, 1,  allocMatrix(REALSXP,1,1));
		var_ret = REAL(VECTOR_ELT(ret, 1));
		var_ret[0] = 0;
	}
	UNPROTECT(4);
	return(ret);
}

SEXP connect_groups_all (SEXP list, SEXP Data, SEXP par,SEXP mean){
	/* par[0] len of data, par[1] dim of data ,p[2]conn type,p[3] 0 no pen, 1 exponencial, 2 lineal, 3 power, 1 pen,p[4] metric,p[5] lineal constant*/
	
	
	//min_val changed for euclidean
	
	int *p;
	SEXP nodes,dist;
	char *visited;
	int e;
	int i,ii,k,j,LEN,*nodes_tmp;
	double *data,val,min_val,*_mean,*dist_tmp;
	char bool,c_one=1;
	int *index,i_one=1;
	double *distance;
	PROTECT(par = coerceVector(par, INTSXP));
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(mean = coerceVector(mean, REALSXP));
	data = REAL(Data);
	p = INTEGER(par);
	_mean = REAL(mean);
	e = p[5];

	visited = (char*) R_alloc(p[0],sizeof(char));
	memset (visited, 0, p[0]);
	
	dist_tmp = (double*) R_alloc(p[0],sizeof(double));
	dist = VECTOR_ELT(list, 0);
	nodes = VECTOR_ELT(list, 1);
	if(p[3]==1){
		for (i=0;i<p[0];i++){
			ii = p[0]*i;
			dist_tmp = (double*) R_alloc(p[0],sizeof(double));
			nodes_tmp = (int*) R_alloc(p[0],sizeof(int));
			memset (visited, 0, p[0]); // is this thread safe ?
			index = VECTOR_ELT(nodes, i);
			distance = VECTOR_ELT(dist, i);
			LEN = length(VECTOR_ELT(dist, i));
			for (j=0;j<LEN;j++) visited[index[j]] = c_one;
			k=0;
			for (j=0;j<p[0];j++){
				nodes_tmp[j]=(j+i_one);
				dist_tmp[j] = (visited[j])? distance[k]:data[j+ii]*exp(data[j+ii]/_mean[0]);
				k += (visited[j])? 1:0;
			}
			SET_VECTOR_ELT(index, i,  nodes_tmp);
			SET_VECTOR_ELT(distance, i,  dist_tmp);
		}
	}
	else if(p[3]==2){
		for (i=0;i<p[0];i++){
			ii = p[0]*i;
			dist_tmp = (double*) R_alloc(p[0],sizeof(double));
			nodes_tmp = (int*) R_alloc(p[0],sizeof(int));
			memset (visited, 0, p[0]); // is this thread safe ?
			index = VECTOR_ELT(nodes, i);
			distance = VECTOR_ELT(dist, i);
			LEN = length(VECTOR_ELT(dist, i));
			for (j=0;j<LEN;j++) visited[index[j]] = c_one;
			k=0;
			for (j=0;j<p[0];j++){
				nodes_tmp[j]=(j+i_one);
				dist_tmp[j] = (visited[j])? distance[k]:data[j+ii]*e;
				k += (visited[j])? 1:0;
			}
			SET_VECTOR_ELT(index, i,  nodes_tmp);
			SET_VECTOR_ELT(distance, i,  dist_tmp);
		}
	}
	else if(p[3]==3){
		for (i=0;i<p[0];i++){
			ii = p[0]*i;
			dist_tmp = (double*) R_alloc(p[0],sizeof(double));
			nodes_tmp = (int*) R_alloc(p[0],sizeof(int));
			memset (visited, 0, p[0]); // is this thread safe ?
			index = VECTOR_ELT(nodes, i);
			distance = VECTOR_ELT(dist, i);
			LEN = length(VECTOR_ELT(dist, i));
			for (j=0;j<LEN;j++) visited[index[j]] = c_one;
			k=0;
			for (j=0;j<p[0];j++){
				nodes_tmp[j]=(j+i_one);
				dist_tmp[j] = (visited[j])? distance[k]:data[j+ii]*intpow((data[j+ii]/_mean[0]),e);
				k += (visited[j])? 1:0;
			}
			SET_VECTOR_ELT(index, i,  nodes_tmp);
			SET_VECTOR_ELT(distance, i,  dist_tmp);
		}
	}
	else{
		Rprintf("only values {1,2,3} apply for (ALL) penalization");
	}
	UNPROTECT(3);
	return(0);
}
