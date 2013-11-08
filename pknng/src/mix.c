#include "def.h"

typedef struct{
	int idx;
	double val;
} Sample;

int compare_doubles (const void *, const void*);
inline int breaks(double *,double *,int *, double *, double *, int *,int *);

void mi_metric_cfun (double *, int *, double *);
SEXP correlation_metric(SEXP , SEXP );
SEXP euc_metric(SEXP , SEXP );
SEXP histo_fun(SEXP , SEXP );
SEXP mi_metric(SEXP , SEXP );
SEXP mi_estimator(SEXP , SEXP );

int  compare (const void *a, const void *b)
{
	const Sample *da = (const Sample *) a;
	const Sample *db = (const Sample *) b;
     
	return (da->val < db->val) - (da->val > db->val);
}

int compare_doubles (const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;
     
	return (*da > *db) - (*da < *db);
}


SEXP correlation_metric(SEXP Data, SEXP par){//par[0] rows (dimension of Data), par[1] cols (amount of Data)

	double *d;
	double *E,*E2;
	int *p;

	int i,k,kk,kkk,j,jjj;
	double Exy,ro,den,*m,inv_p0;
	SEXP mimat;
	
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(par = coerceVector(par, INTSXP));
	
	d = REAL(Data);p = INTEGER(par);

	PROTECT(mimat = allocMatrix(REALSXP, p[1],p[1]));
	m = REAL(mimat);
	inv_p0 = 1./p[0];
	E = (double*) R_alloc(p[1],sizeof(double));
	E2 = (double*) R_alloc(p[1],sizeof(double));
	for (k=0;k<p[1];k++){
		E[k] = 0;E2[k] = 0;
		kkk = k*p[0];
		for (i = 0;i<p[0];i++){
		E[k] += d[kkk+i]; E2[k] += d[kkk+i] * d[kkk+i];
		}
		E[k] *= inv_p0; E2[k] *= inv_p0;
	}
	for (k=0;k<p[1];k++){
		m[k*p[1]+k] = 0;
		kk = k*p[1];kkk = k*p[0];
		for (j=k+1;j<p[1];j++){
			m[kk+j] = 1.;
			den = (E2[k]-E[k]*E[k])*(E2[j]-E[j]*E[j]);
			Exy = 0;
			jjj = j*p[0];
			for (i = 0;i<p[0];i++)	Exy += d[kkk+i]*d[jjj+i];
			Exy *= inv_p0;
			ro = Exy - E[k]*E[j];
			ro /= (den > 0) ? sqrt(den):DBL_MIN;
			m[kk+j] -=  ro;
			m[j*p[1]+k] = m[kk+j];
		}
	}
	UNPROTECT(3);
	return(mimat);
}

SEXP euc_metric(SEXP Data, SEXP par){//par[0] rows (dimension of Data), par[1] cols (amount of Data)

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

SEXP histo_fun(SEXP Data, SEXP par){//Data vector, par[0] Data length

	double *d,*dclon,*duniques,*bs;
	double inter,MPB,ftmp;
	int *p;
	int *fduniques;
	int j,i,uniques,bins,ac,*fbins;

	SEXP breaks;
	
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(par = coerceVector(par, INTSXP));
	
	d = REAL(Data);p = INTEGER(par);
	
	bins = (int)(1.5 +  log2 (p[0]));
	
	dclon = (double*) R_alloc(p[0],sizeof(double));
	duniques = (double*) R_alloc(p[0],sizeof(double));
	fduniques = (int*) R_alloc(p[0],sizeof(int));
		
	for (j=0;j<p[0];j++) {dclon[j]=d[j];fduniques[j]=0;}
	qsort (dclon, p[0], sizeof(double), compare_doubles); //sort data 
	
	duniques[0] = dclon[0];
	uniques = 0;
	fduniques[0] = 1;
	for (j=1;j<p[0];j++) {
		if (dclon[j] != duniques[uniques]){
			uniques++;
			duniques[uniques] = dclon[j];
			fduniques[uniques] = 1;
		}
		else fduniques[uniques]++;
	}

	uniques++;
	
	if (uniques<=bins) {
		bins = uniques;
		PROTECT(breaks = allocVector(REALSXP, (bins+1)));
		bs = REAL(breaks);
		inter = duniques[uniques-1] - duniques[0];
		bs[0] = duniques[0];
		bs[bins] = dclon[p[0]-1];
		ac = 0;
		for (j=0;j<uniques-1;j++){
			ac += fduniques[j];
			bs[j+1] = duniques[0]+ac*(inter/p[0]);
			if (bs[j+1] > duniques[j+1]) bs[j+1] = duniques[j+1]-((duniques[j+1] - duniques[j])/10000);
		}
	}
	else{
		PROTECT(breaks = allocVector(REALSXP, (bins+1)));
		fbins = (int*) R_alloc(bins,sizeof(int));
		bs = REAL(breaks);
		for (j=0;j<bins;j++) {fbins[j] = 0;bs[j]=0.;}
		bs[bins]=0.;
		inter = dclon[p[0]-1] - dclon[0];
		MPB = (((double)p[0])/bins);
		bs[0] = duniques[0];
		bs[bins] = dclon[p[0]-1];
		i = 0;
		j = 0;
		ac = 0;
		for(i=0;i<uniques;){
			if (ac < (MPB+(MPB*j))){
				fbins[j] += fduniques[i];
				ac += fduniques[i];
				i++;
			}
			else {j++;}
			if (j==uniques) break;
		}
		ac = 0;
		for (j=0;j<bins-1;j++){
			ac += fbins[j];
			bs[j+1] = duniques[0]+ac*(inter/p[0]);
			if (bs[j+1] > dclon[ac+1]) bs[j+1] = dclon[ac+1]-((dclon[ac+1] - dclon[ac])/10000);
		}
	}
	
	
	UNPROTECT(3);
	return(breaks);
}

inline int breaks(double *d, double *breaks , int *p, double *dclon, double *duniques, int *fduniques, int *fbins){
// d is a pointer to double (data), breaks is a pointer to double (breas for histogram), 
// p is a pointer to int (parameters) p[0] d length, p[1] = bins

	double *bs;
	double inter,MPB,ftmp;
	int j,i,uniques,bins,ac;

	bins = p[1];
	for (j=0;j<p[0];j++) {
		dclon[j]=d[j];
		fduniques[j]=0;
	}
	for (j=0;j<bins+1;j++) {breaks[j] = 0.;/*Rprintf("%i / %i - %i ",&dclon[j],&breaks[j],j);*/}
	qsort (dclon, p[0], sizeof(double), compare_doubles); //sort data 
	duniques[0] = dclon[0];
	uniques = 0;
	fduniques[0] = 1;
	for (j=1;j<p[0];j++) {
		if (dclon[j] != duniques[uniques]){
			uniques++;
			duniques[uniques] = dclon[j];
			fduniques[uniques] = 1;
		}
		else fduniques[uniques]++;
	}
	uniques++;
	if (uniques<=bins) {
		bins = uniques;
		bs = breaks;
		inter = duniques[uniques-1] - duniques[0];
		bs[0] = duniques[0];
		bs[bins] = dclon[p[0]-1];
		ac = 0;
		for (j=0;j<uniques-1;j++){
			ac += fduniques[j];
			bs[j+1] = duniques[0]+ac*(inter/p[0]);
			if (bs[j+1] > duniques[j+1]) bs[j+1] = duniques[j+1]-((duniques[j+1] - duniques[j])/10000);
		}
	}
	else{
		bs = breaks;
		for (j=0;j<bins;j++) {fbins[j] = 0;bs[j]=0.;}
		bs[bins]=0.;
		inter = dclon[p[0]-1] - dclon[0];
		MPB = (((double)p[0])/bins);
		bs[0] = duniques[0];
		bs[bins] = dclon[p[0]-1];
		i = 0;
		j = 0;
		ac = 0;
		for(i=0;i<uniques;){
			if (ac < (MPB+(MPB*j))){
				fbins[j] += fduniques[i];
				ac += fduniques[i];
				i++;
			}
			else j++;
			if (j==uniques) break;
		}
		ac = 0;
		for (j=0;j<bins-1;j++){
			ac += fbins[j];
			bs[j+1] = duniques[0]+ac*(inter/p[0]);
			if (bs[j+1] > dclon[ac+1]) bs[j+1] = dclon[ac+1]-((dclon[ac+1] - dclon[ac])/10000);
		}
	}
	return(bins);
}

SEXP mi_metric(SEXP Data, SEXP par){
//par[1] rows (dimension of Data), par[0] cols (amount of Data), par[2] bins, p[3] = 1 matrix, 0  vector

	
	double *xbreak, *ybreak;
	double *ret,*d;
	double *px,*py,*pxy;
	double Hx,Hy,H;
	double *dclon,*duniques;
	double cte;
	
	int *fduniques,*fbins;
	int *param;
	int xbin,ybin,*p,bins;
	int i,j,u,v;
	int xidx,yidx,*xindexs;
	int bbb;
	int may=0,men=0;
	SEXP mimat;
	
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(par = coerceVector(par, INTSXP));
	
	d = REAL(Data);p = INTEGER(par);
	
	if (p[2] == 0) bins = (int)(1.5 +  log2 (p[1]));
	else bins = p[2];
	
	//cte = log2(bins) ;
	if (p[3]) PROTECT(mimat = allocMatrix(REALSXP,p[0],p[0]));
	else PROTECT(mimat = allocVector(REALSXP,(p[0]*p[0])));
	ret = REAL(mimat);	
	param = (int*) R_alloc(2,sizeof(int));
	bbb = (bins*bins);
			
	param[0] = p[1];
	param[1] = bins;
	ybreak = (double*) R_alloc((bins+1),sizeof(double));
	xbreak = (double*) R_alloc((bins+1),sizeof(double));
	pxy = (double*) R_alloc(bbb,sizeof(double));
	px = (double*) R_alloc(bins,sizeof(double));
	py = (double*) R_alloc(bins,sizeof(double));

	dclon = (double*) R_alloc(p[1],sizeof(double));
	duniques = (double*) R_alloc(p[1],sizeof(double));
	fduniques = (int*) R_alloc(p[1],sizeof(int));
	xindexs = (int*) R_alloc(p[1],sizeof(int));
	fbins = (int*) R_alloc(bins,sizeof(int));

	for(i=0;i<p[0];i++){
		xbin = breaks(&d[i*p[1]], xbreak , param, dclon, duniques, fduniques, fbins);
		ret[i+i*p[0]] = 0;
		for(u=0;u<bins;u++) px[u] = 0;
		for(u=0;u<p[1];u++){
			xindexs[u] = 0; // stores x indexes so you don't have to recalculate for each y
			for(v=1;v<xbin;v++) if (d[i*p[1]+u] > xbreak[v]) xindexs[u]++;
			px[xindexs[u]]+=1.;
		}
		Hx = 0;
		for(v=0;v<bins;v++) {px[v]/=(double)p[1];Hx-=(px[v]>0)?px[v]*log2(px[v]):0;}
		for(j=i+1;j<p[0];j++){
			ret[j+j*p[0]] = 0;
			for(u=0;u<bbb;u++) pxy[u] = 0;
			for(u=0;u<bins;u++) py[u] = 0;
			ybin = breaks(&d[j*p[1]], ybreak , param, dclon, duniques, fduniques, fbins);
			for(u=0;u<p[1];u++){
				yidx = 0;
				for(v=1;v<ybin;v++) if (d[j*p[1]+u] > ybreak[v]) yidx++;
				pxy[xindexs[u]+yidx*bins]+=1.;py[yidx]+=1.;
			}
			Hy = 0;
			for(v=0;v<bins;v++) {py[v]/=(double)p[1];Hy-=(py[v]>0)?py[v]*log2(py[v]):0;}
			for(u=0;u<bbb;u++) pxy[u]/=(double)p[1];
			ret[i+j*p[0]] = 0;
			for(u=0;u<bins;u++){//u=x;v=y
				for(v=0;v<bins;v++){
			ret[i+j*p[0]] += (pxy[u+v*bins]>0 && py[v]>0 && px[u]>0 )? pxy[u+v*bins]*(log2(pxy[u+v*bins])-log2(py[v])-log2(px[u])) :0;
				}
			}
			H = Hx+Hy;
			ret[i+j*p[0]] /= (H != 0)?H:1;
			ret[i+j*p[0]] = ret[j+i*p[0]] =  1.-ret[i+j*p[0]] ;
		}
	}
	UNPROTECT(3);
	return(mimat);
}

SEXP mi_estimator(SEXP Data, SEXP par){
//	par[0] rows (dimension of Data) , par[1] col (amount of Data)

	
	double *d,**fxy,*ret;
	double dif,ddif,dddif,H,H_sq,H2d,H_sq2d,*mean,RBF_cte,RBF_cte2d,cte_pow2d,cte_pow,sfx,sfy,sfxy;
	double Hxy,upper,lower;
	double *f_1d,*sigma_2d,*sigma_,**S,**S_inv,determinant;
	double raiz_2pi,_2pi,half,_invp0_sq,_inv,MEAN_1,MEAN_2,d1,d2,H_sq2d_,H_sq_,ro,max_err;
	int *p;
	int i,j,k,u,jj,iii,ii,kk,jjj,*cnt,_2by2;
	int Tnum,MAX_T,T;
	SEXP mimat;
	
	_2by2 = 4;
	
	// From openmp version
	MAX_T = 1;
	
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(par = coerceVector(par, INTSXP));
	d = REAL(Data);p = INTEGER(par);
	PROTECT(mimat = allocMatrix(REALSXP,p[1],p[1]));  // covariance matrix p[1] x p[1]  ;  ret is using first as covariance matrix a then a MI matrix
	ret = REAL(mimat);

	fxy = (double**) R_alloc(MAX_T,sizeof(double*));
	for(i=0;i<MAX_T;i++) fxy[i] = (double*) R_alloc(p[0],sizeof(double));
	mean = (double*) R_alloc((p[1]),sizeof(double));
	cnt = (int*) R_alloc(MAX_T,sizeof(int));
	
	
	S = (double**) R_alloc((MAX_T),sizeof(double*)); // covariance matrix p[1] x p[1]
	S_inv = (double**) R_alloc((MAX_T),sizeof(double*));
	for(i=0;i<MAX_T;i++) {S[i] = (double*) R_alloc((_2by2),sizeof(double)); S_inv[i] = (double*) R_alloc((_2by2),sizeof(double));}
	f_1d = (double*) R_alloc(p[1]*p[0],sizeof(double));
	
	cte_pow = powl((4./(3*p[0])), 0.2); //hopt 1d
	cte_pow2d = powl( (1./(p[0])), (1/6.)); //hopt 2d <---------------- (1/6.) valor exacto segun Steuer 02 Bioinformatics
	raiz_2pi = sqrt(2*M_PI);
	_2pi = 2*M_PI;
	half = 1./2;
	_invp0_sq =  1./(p[0]*(p[0]+1)*0.5);
	determinant = 1.;
	
	H_sq2d_ = H_sq2d = 1/(2*cte_pow2d*cte_pow2d);
	H_sq_ = H_sq = 1/(2 * cte_pow* cte_pow);

	// From openmp version
	Tnum = 0;

	for(i=0;i<p[1];i++){
		ii = i*p[0];
		mean[i] = 0.;
		for(k=0;k<p[0];k++) for(u=k+1;u<p[0];u++) mean[i] += sqrt((d[k+ii]-d[u+ii])*(d[k+ii]-d[u+ii]));
		mean[i] *= _invp0_sq;
	}
	for(i=0;i<p[1];i++){
		iii = i*p[1];
		ii = i*p[0];
		MEAN_1= mean[i];
		for(j=i;j<p[1];j++){
			ret[iii+j]=0.;
			jj = j*p[0];
			jjj = j*p[1];
			MEAN_2= mean[j];
			for(k=0;k<p[0];k++){
				for(u=k;u<p[0];u++){
					dif = (sqrt((d[k+ii]-d[u+ii])*(d[k+ii]-d[u+ii]))-MEAN_1);
					ddif = (sqrt((d[k+jj]-d[u+jj])*(d[k+jj]-d[u+jj]))-MEAN_2);
					ret[iii+j] += dif*ddif;
				}
			}
			ret[iii+j] *= _invp0_sq;
			ret[jjj+i] = ret[iii+j];
		}
	}

	for(i=0;i<p[1];i++) {ii=i*p[0];for(k=0;k<p[0];k++) f_1d[ii+k] = 0.;}
	for(i=0;i<p[1];i++){
		ii = i*p[0];
		RBF_cte = 1/(p[0]*sqrt(ret[i*p[1]+i])*cte_pow*raiz_2pi); // reuse mean array (as standard deviation inverse) 1/sqrt(ret[i+i]) where ret[i+i] is variance for patern i
		H_sq = 1/(2 * cte_pow* cte_pow*ret[i*p[1]+i]);
		/*	F_1d	*/ 
		for(k=0;k<p[0];k++){
			for(u=0;u<p[0];u++){
				dif = (d[k+ii]-d[u+ii]);
				f_1d[ii+k] += exp((-dif*dif*H_sq));	//  ### x
			}
			f_1d[ii+k] *= RBF_cte;
		}
	}
	for(i=0;i<p[1];i++){
			// ******* careful ret has covariance *******
		iii = i*p[1];
		ii = i*p[0];
		for(j=i+1;j<p[1];j++){
			jjj = j*p[1];
			jj = j*p[0];
			S[Tnum][0] = ret[iii+i]; S[Tnum][1] = ret[iii+j]; 
			S[Tnum][2] = ret[jjj+i]; S[Tnum][3] = ret[jjj+j]; /* S[1] == S[2] = ret[iii+j] == ret[jjj+i]*/
			determinant = (S[Tnum][0]*S[Tnum][3])-(S[Tnum][1]*S[Tnum][2]);
			
			S_inv [Tnum][0] = S[Tnum][3]/determinant; S_inv [Tnum][1] = -S[Tnum][1]/determinant;
			S_inv [Tnum][2] = -S[Tnum][2]/determinant; S_inv [Tnum][3] = S[Tnum][0]/determinant;
			ret[iii+j]=0.;
			RBF_cte2d = 1/(p[0]* sqrt(determinant)* cte_pow2d* cte_pow2d *_2pi);
			RBF_cte2d = (determinant<0.1 && ro>0.85) ? 0.5*RBF_cte2d:RBF_cte2d;
			for(k=0;k<p[0];k++){
				fxy[Tnum][k]=0.;
				for(u=0;u<p[0];u++){
					d1 = d[k+ii]-d[u+ii]; d2 = d[k+jj]-d[u+jj];
					dif = (d2*d2*S_inv[Tnum][3]);
					dddif = (2*d1*d2*S_inv[Tnum][1]);
					ddif = (d1*d1*S_inv[Tnum][0]);
					ddif += dif+dddif; 
					ddif *= H_sq2d;
					fxy[Tnum][k] += exp(-ddif);	//  ### xy 
				}
				fxy[Tnum][k] *= RBF_cte2d;
			}
			Hxy = 0.;// Estimador de H(x,y)
			ret[iii+j] = 0.;
			for(k=0;k<p[0];k++) {
				ddif = f_1d[ii+k]*f_1d[jj+k];
				dif = (fxy[Tnum][k]<1)?log2(fxy[Tnum][k]):0;
				ddif = -log2(ddif);
				ret[iii+j]+= (dif+ddif);
				Hxy -=  dif;
			}
			ret[iii+j] /= (Hxy!=0)? Hxy:1; 
			ret[iii+j] = 1.-ret[iii+j]; // ret = I(x,y)/H(x,y) ### Informacion mutua / Entropia conjunta
			max_err = (ret[iii+j]>1 && ret[iii+j]>max_err)?ret[iii+j]:max_err;
			ret[iii+j] = (ret[iii+j]>1)? 1:ret[iii+j];
			ret[jjj+i] = ret[iii+j];
		}
	}

	for(i=0;i<p[1];i++) ret[i*p[1]+i]=0.;
	UNPROTECT(3);
	return(mimat);
}

SEXP mutualInfo (SEXP Data, SEXP par){
//par[1] rows (dimension of Data), par[0] cols (amount of Data), par[2] bins

	
	double *xbreak, *ybreak;
	double *ret,*d;
	double *px,*py,*pxy;
	double *dclon,*duniques;
	double cte;
	
	int *fduniques,*fbins;
	int *param;
	int xbin,ybin,*p,bins;
	int i,j,u,v;
	int xidx,yidx,*xindexs;
	int bbb;
	int may=0,men=0;
	SEXP mimat;
	
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(par = coerceVector(par, INTSXP));
	
	d = REAL(Data);p = INTEGER(par);
	
	if (p[2] == 0) bins = (int)(1.5 +  log2 (p[1]));
	else bins = p[2];
	
	cte = log2(bins) ;
	PROTECT(mimat = allocMatrix(REALSXP,p[0],p[0]));
	ret = REAL(mimat);	
	param = (int*) R_alloc(2,sizeof(int));
	bbb = (bins*bins);
			
	param[0] = p[1];
	param[1] = bins;
	ybreak = (double*) R_alloc((bins+1),sizeof(double));
	xbreak = (double*) R_alloc((bins+1),sizeof(double));
	pxy = (double*) R_alloc(bbb,sizeof(double));
	px = (double*) R_alloc(bins,sizeof(double));
	py = (double*) R_alloc(bins,sizeof(double));

	dclon = (double*) R_alloc(p[1],sizeof(double));
	xindexs = (int*) R_alloc(p[1],sizeof(int));
	duniques = (double*) R_alloc(p[1],sizeof(double));
	fduniques = (int*) R_alloc(p[1],sizeof(int));
	fbins = (int*) R_alloc(bins,sizeof(int));
	
	for(i=0;i<p[0];i++){
		xbin = breaks(&d[i*p[1]], xbreak , param, dclon, duniques, fduniques, fbins);
		ret[i+i*p[0]] = 0;
		for(u=0;u<bins;u++) px[u] = 0;
		for(u=0;u<p[1];u++){
			xindexs[u] = 0;
			for(v=1;v<xbin;v++) if (d[i*p[1]+u] > xbreak[v]) xindexs[u]++;
			px[xindexs[u]]+=1.;
		}
		for(v=0;v<bins;v++) px[v]/=(double)p[1];
		for(j=i+1;j<p[0];j++){
			ret[j+j*p[0]] = 0;
			for(u=0;u<bbb;u++) pxy[u] = 0;
			for(u=0;u<bins;u++) py[u] = 0;
			ybin = breaks(&d[j*p[1]], ybreak , param, dclon, duniques, fduniques, fbins);
			for(u=0;u<p[1];u++){
				yidx = 0;
				for(v=1;v<ybin;v++) if (d[j*p[1]+u] > ybreak[v]) yidx++;
				pxy[xindexs[u]+yidx*bins]+=1.;py[yidx]+=1.;
			}
			for(v=0;v<bins;v++) py[v]/=(double)p[1];
			for(u=0;u<bbb;u++) pxy[u]/=(double)p[1];
			ret[i+j*p[0]] = 0;
			for(u=0;u<bins;u++){//u=x;v=y
				for(v=0;v<bins;v++){
					ret[i+j*p[0]] += (pxy[u+v*bins]>0 && py[v]>0 && px[u]>0 )? pxy[u+v*bins]*(log2(pxy[u+v*bins])-log2(py[v])-log2(px[u])) :0;
				}
			}
			ret[j+i*p[0]] = ret[i+j*p[0]];
		}
	}
	UNPROTECT(3);
	return(mimat);
}

void mi_metric_cfun (double *d, int *p, double *ret){
//p[1] rows (dimension of Data), p[0] cols (amount of Data), p[2] bins, p[3] 0 = nxn matrix 1=triang matrix
// m the matrix where we store the results, according to the value of p[3] is the size it has to be. m is allocated before calling histo
	
	double *xbreak, *ybreak;
	double *px,*py,*pxy;
	double Hx,Hy,H;
	double *dclon,*duniques;
	
	int *fduniques,*fbins;
	int *param;
	int xbin,ybin,bins,*xindexs;
	int i,j,u,v;
	int xidx,yidx;
	int bbb;
	int lencor = 0;
	
	if (p[2] == 0) bins = (int)(1.5 +  log2 (p[1]));
	else bins = p[2];
	
	param = (int*) R_alloc(2,sizeof(int));
	bbb = (bins*bins);
			
	param[0] = p[1];
	param[1] = bins;
	
	xindexs = (int*) R_alloc(p[1],sizeof(int));
	ybreak = (double*) R_alloc((bins+1),sizeof(double));
	xbreak = (double*) R_alloc((bins+1),sizeof(double));
	pxy = (double*) R_alloc(bbb,sizeof(double));
	px = (double*) R_alloc(bins,sizeof(double));
	py = (double*) R_alloc(bins,sizeof(double));
	
	dclon = (double*) R_alloc(p[1],sizeof(double));
	duniques = (double*) R_alloc(p[1],sizeof(double));
	fduniques = (int*) R_alloc(p[1],sizeof(int));
	fbins = (int*) R_alloc(bins,sizeof(int));

	if(p[3]){
		for(i=0;i<p[0];i++){
			lencor = i*(i+1)/2;
			ret[i*p[0]+i-lencor] = 0.;
			xbin = breaks(&d[i*p[1]], xbreak , param, dclon, duniques, fduniques, fbins);
			for(u=0;u<bins;u++) px[u] = 0;
			for(u=0;u<p[1];u++){
				xindexs[u] = 0;
				for(v=1;v<xbin;v++) if (d[i*p[1]+u] > xbreak[v]) xindexs[u]++;
				px[xindexs[u]]+=1.;
			}
			for(v=0;v<bins;v++) {
				px[v]/=(double)p[1]; Hx-=(px[v]>0)?px[v]*log2(px[v]):0;
			}
			for(j=i+1;j<p[0];j++){
				for(u=0;u<bbb;u++) pxy[u] = 0;
				for(u=0;u<bins;u++) py[u] = 0;
				ybin = breaks(&d[j*p[1]], ybreak , param, dclon, duniques, fduniques, fbins);
				for(u=0;u<p[1];u++){
					yidx = 0;
					for(v=1;v<ybin;v++) if (d[j*p[1]+u] > ybreak[v]) yidx++;
					pxy[xindexs[u]+yidx*bins]+=1.;
					py[yidx]+=1.;
				}
				Hx = 0;Hy = 0;
				for(v=0;v<bins;v++) {
					py[v]/=(double)p[1]; Hy-=(py[v]>0)?py[v]*log2(py[v]):0;
				}
				
				for(u=0;u<bbb;u++) pxy[u]/=(double)p[1];
				ret[i+j*p[0]-lencor] = 0.;
				for(u=0;u<bins;u++){//u=x;v=y
					for(v=0;v<bins;v++){
			ret[i+j*p[0]-lencor] += (pxy[u+v*bins]>0 && py[v]>0 && px[u]>0 )? pxy[u+v*bins]*(log2(pxy[u+v*bins])-log2(py[v])-log2(px[u])) :0;
					}
				}
				H = Hx+Hy;if (H == 0) H = 1;
				ret[i+j*p[0]-lencor] *= 2;
				
				ret[i+j*p[0]-lencor]  = 1.-(ret[i+j*p[0]-lencor]/H);
			}
		}
	}
	else{
		for(i=0;i<p[0];i++){
			xbin = breaks(&d[i*p[1]], xbreak , param, dclon, duniques, fduniques, fbins);
			ret[i+i*p[0]] = 0;
			for(u=0;u<bins;u++) px[u] = 0;
			for(u=0;u<p[1];u++){
				xindexs[u] = 0; // stores x indexes so you don't have to recalculate for each y
				for(v=1;v<xbin;v++) if (d[i*p[1]+u] > xbreak[v]) xindexs[u]++;
				px[xindexs[u]]+=1.;
			}
			Hx = 0;
			for(v=0;v<bins;v++) {px[v]/=(double)p[1];Hx-=(px[v]>0)?px[v]*log2(px[v]):0;}
			for(j=i+1;j<p[0];j++){
				ret[j+j*p[0]] = 0;
				for(u=0;u<bbb;u++) pxy[u] = 0;
				for(u=0;u<bins;u++) py[u] = 0;
				ybin = breaks(&d[j*p[1]], ybreak , param, dclon, duniques, fduniques, fbins);
				for(u=0;u<p[1];u++){
					yidx = 0;
					for(v=1;v<ybin;v++) if (d[j*p[1]+u] > ybreak[v]) yidx++;
					pxy[xindexs[u]+yidx*bins]+=1.;py[yidx]+=1.;
				}
				Hy = 0;
				for(v=0;v<bins;v++) {py[v]/=(double)p[1];Hy-=(py[v]>0)?py[v]*log2(py[v]):0;}
				for(u=0;u<bbb;u++) pxy[u]/=(double)p[1];
				ret[i+j*p[0]] = 0;
				for(u=0;u<bins;u++){//u=x;v=y
					for(v=0;v<bins;v++){
						ret[i+j*p[0]] += (pxy[u+v*bins]>0 && py[v]>0 && px[u]>0 )? pxy[u+v*bins]*(log2(pxy[u+v*bins])-log2(py[v])-log2(px[u])) :0;
					}
				}
				H = Hx+Hy;if (H == 0) H = 1;ret[i+j*p[0]] *= 2;
				ret[i+j*p[0]] = ret[j+i*p[0]] =  1.-(ret[i+j*p[0]]/H) ;
			}
		}
	}
}

SEXP rbf_dot_multiscaled(SEXP Data, SEXP par,SEXP Sigma){
  //par[0] rows (dimension of Data), par[1] cols (amount of Data)

	double *d;
	int *p;

	int i,k,j,matlen,idx,kk,jj,kkk;
	double *m,dif,*sigma,done=1.,dtwo=2.;
	SEXP mimat;
	
	PROTECT(Data = coerceVector(Data, REALSXP));
	PROTECT(Sigma = coerceVector(Sigma, REALSXP));
	PROTECT(par = coerceVector(par, INTSXP));
	
	d = REAL(Data);p = INTEGER(par);sigma = REAL(Sigma);
// 	Rprintf("%i %i %lf\n",p[0],p[1],sigma[0]);

	PROTECT(mimat = allocMatrix(REALSXP, p[1],p[1]));
	m = REAL(mimat);
	for (k=0;k<p[1];k++){
// 		m[k*p[1]+k] = 0;
		kk = k*p[0];kkk = k*p[1];
		for (j=k;j<p[1];j++){
			jj = j*p[0];
			dif = 0.;m[kkk+j] = 0;
			for (i = 0;i<p[0];i++){
				dif = d[kk+i]- d[jj+i];
				m[kkk+j] += dif*dif;
			}
			m[j*p[1]+k] = m[kkk+j] = exp( -m[kkk+j]/(sigma[k]*sigma[j]) );
		}
	}

	UNPROTECT(4);
	return(mimat);
}
