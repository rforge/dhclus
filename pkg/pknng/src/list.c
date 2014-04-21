#include "def.h"

SEXP
printl (SEXP list)
{
	SEXP elmt,s,elmtcp;
	int listlen,elmtlen;
	double *var,*var2;
	int i,j;
	
	PROTECT(s = allocVector(VECSXP, 2));
	Rprintf("%i\n",length(list));
	listlen = length(list);
	elmt = VECTOR_ELT(list, 0);
	elmtlen = length(elmt);
	
	var = REAL(elmt);
/*	Rprintf("%i\n",elmtlen);
	for (i = 0;i<elmtlen;i++) Rprintf("%lf ",var[i]);
	Rprintf("\n");*/
	SET_VECTOR_ELT(s, 0,  allocMatrix(REALSXP, elmtlen,2));
	SET_VECTOR_ELT(s, 1,  allocVector(REALSXP, elmtlen));
	elmtcp = VECTOR_ELT(s, 0);
	var2 = REAL(elmtcp);
	
	for (j = 0;j<2;j++) for (i = 0;i<elmtlen;i++) var2[i+j*elmtlen] = var[i];
	
	elmtcp = VECTOR_ELT(s, 1);
	var2 = REAL(elmtcp);
	for (i = 0;i<elmtlen;i++) var2[elmtlen-i-1] = var[i];

	UNPROTECT(1);
	return (s);
} 

inline int add_node(lnode *root, unsigned long int node, double d){
	lnode *tmp,*n;
	
/*	Rprintf("root: %i\n",root->tag);
	Rprintf("root cnt: %i\n",root->cnt);
	Rprintf("to: %i\n",node);*/
	
	tmp = (lnode*) R_alloc(1,sizeof(lnode));
	tmp->itag = node;
	tmp->jtag = -1;
	tmp->cnt = -1;
	tmp->mark = 0;
	tmp->dist = d;
	tmp->next = NULL;
	tmp->down = NULL;
	if (root->next == NULL) root->next = tmp;
	else{
		n = root->next;
		while(n->next != NULL) n = n->next;
		n->next = tmp;
	}
	root->cnt++;
	return(root->cnt);
}

inline void add_ij_node(lnode *root, unsigned long int inode,unsigned long int jnode, double d){
	lnode *tmp,*n;
	
/*	Rprintf("root: %i\n",root->tag);
	Rprintf("root cnt: %i\n",root->cnt);
	Rprintf("to: %i\n",node);*/
	
	tmp = (lnode*) R_alloc(1,sizeof(lnode));
	tmp->itag = inode;
	tmp->jtag = jnode;
	tmp->cnt = -1;
	tmp->mark = 0;
	tmp->dist = d;
	tmp->next = NULL;
	tmp->down = NULL;
	tmp->up = NULL;
	if (root->next == NULL) root->next = tmp;
	else{
		n = root->next;
		while(n->next != NULL) n = n->next;
		n->next = tmp;
	}
	root->cnt++;
}


inline void init_root_node(lnode *root, unsigned long int tag){
	root->itag = tag;
	root->jtag = -1;
	root->cnt = 0;
	root->mark = 0;
	root->cnt_root = 0;
	root->next = NULL;
	root->down = NULL;
	root->up = NULL;
	root->dist = -1.;
}

inline lnode * add_root_node(lnode *root, unsigned long int tag){
	lnode *tmp,*n;
	
/*	Rprintf("root: %i\n",root->tag);
	Rprintf("root cnt: %i\n",root->cnt);
	Rprintf("to: %i\n",node);*/
	
	tmp = (lnode*) R_alloc(1,sizeof(lnode));
	tmp->itag = tag;
	tmp->jtag = tag;
	tmp->cnt = 0;
	tmp->mark = 0;
	tmp->cnt_root = 0;
	tmp->dist = -1.;
	tmp->next = NULL;
	tmp->down = NULL;
	tmp->up = NULL;
	if (root->down == NULL){
		 root->down = tmp;
		 tmp->up = root;
	}
	else{
		n = root->down;
		while(n->down != NULL) n = n->down;
		n->down = tmp;
		tmp->up = n;
	}
	root->cnt_root++;
	return(tmp);
}

// inline void merge(lnode *root,lnode *min_root){
// 	
// 	lnode *n;
// 	
// 	if (min_root->down != NULL) min_root->down->up = min_root->up;
// 	min_root->up->down = min_root->down;
// 	
// 	min_root->down = NULL;
// 	min_root->up = NULL;
// 	n = min_root;
// 	while(n->next != NULL) n = n->next;
// 	n->next = root->next;
// 	root->next = min_root;
// 	root->cnt += min_root->cnt;
// 	min_root->jtag = -1;
// 	min_root->cnt = 0;
// 	min_root->mark = 0;
// 	min_root->cnt_root = 0;
// 	
// }

inline void merge(lnode *root,lnode *min_root){
	
	lnode *n;
	
	if (min_root->up == NULL) {n = min_root;min_root = root;root = n;}
	
	if (min_root->down != NULL) min_root->down->up = min_root->up;
	if (min_root->up != NULL) min_root->up->down = min_root->down;
	

	min_root->down = NULL;
	min_root->up = NULL;
	n = min_root;
	
	while(n->next != NULL) n = n->next;
	
	n->next = root->next;
	root->next = min_root;
	root->cnt += (min_root->cnt) ? min_root->cnt:1;
	
	min_root->jtag = -1;
	min_root->cnt = 0;
	min_root->mark = 0;
	min_root->cnt_root = 0;
	
}
