#include <math.h>
#include "nrutil.h"
#include "structs.h"

#define NRANSI
#define TO 0
#define FROM 1

void mrqcof(float **y, float **sig, struct image *psf, float a[], int ia[],
    int ma, double **alpha, double beta[], double *chisq, 
    void (*funcs)(int, int, float *, float [], int[], int, 
    struct fitpars *fpar), struct fitpars *fpar, struct convpars *cpar);
void covsrt(double **covar, int ma, int ia[], int mfit);
void gaussj(double **a, int n, double **b, int m);
void copy_pars (float a[], int ia[], struct fitpars *fpar, int dir);
void par_monitor (int ma, float a[], int ia[], float atry[], double da[],
    struct fitpars *fpar);
void constraints (int ma, int mfit, float atry[], float a[], int ia[],
    struct cons *constr, double *da, float *aorig);
void jacobi (double **a, int n, double d[], double **v, int *nrot);

int mrqmin(float **y, float **sig, struct image *psf, int ma, double **covar, 
    double **alpha, double *chisq, void (*funcs)(int, int, float *, float [], 
    int [], int, struct fitpars *), float *alamda, struct fitpars *fpar, 
    struct cons *constr, struct convpars *cpar, double *sigma)
{
        void errcalc (double *sigma, double **alpha, int ma, int *ia, int mfit);

	int j,k,l;
	static int mfit;
	static double *da,**oneda,*beta,ochisq;
        static float *atry, *aorig;
	float *a;
	int *ia;

        /* Kludge -- copy from fpar->a to a */
        a = vector (1, ma);
        ia = ivector (1, ma);
        copy_pars (a, ia, fpar, FROM);

	if (*alamda < 0.0) {
		aorig=vector(1,ma);
 		copy_pars (aorig, ia, fpar, FROM);
		atry=vector(1,ma);
		beta=dvector(1,ma);
		da=dvector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]==1)
				mfit++;
		oneda=dmatrix(1,mfit,1,1);
		*alamda=0.01;
		mrqcof(y,sig,psf,a,ia,ma,alpha,beta,chisq,funcs,fpar,cpar);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}

	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda));
		oneda[j][1]=beta[j];
	}

	gaussj(covar,mfit,oneda,1);

	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {				/* Convergence!      */
                errcalc (sigma, alpha, ma, ia, mfit);   /* calculate uncert. */
		covsrt(covar,ma,ia,mfit);
		covsrt(alpha,ma,ia,mfit);
		free_dmatrix(oneda,1,mfit,1,1);
		free_dvector(da,1,ma);
		free_dvector(beta,1,ma);
		free_vector(atry,1,ma);
		free_vector (a, 1, ma);
		free_ivector (ia, 1, ma);
		return;
	}

        par_monitor (ma, a, ia, atry, da, fpar);

	constraints (ma, mfit, atry, a, ia, constr, da, aorig);

        copy_pars (atry, ia, fpar, TO);

	mrqcof(y,sig,psf,atry,ia,ma,covar,da,chisq,funcs,fpar,cpar);

	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=1;j<=mfit;j++) {
			for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
			beta[j]=da[j];
		}
		for (l=1;l<=ma;l++)
			a[l]=atry[l];

	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
        copy_pars (a, ia, fpar, TO);
	free_vector (a, 1, ma);
	free_ivector (ia, 1, ma);

	return (0);
}
#undef NRANSI


/************************************************************************/

float gammp (float a, float x);

void errcalc (double *sigma, double **alpha, int ma, int *ia, int mfit)
{
    void sortsig(double *sigma, int ma, int ia[], int mfit);
    float rootfind(float (*func)(float, float), float nu, float x1, float x2,
								float xacc);

    double *d, **v;
    int nrot, m, n;
    float delchi2;

    d = dvector (1, mfit);
    v = dmatrix (1, mfit, 1, mfit);

    for (m=1; m <= mfit; m++) {
	for (n=1; n <= mfit; n++)
	    v[n][m] = 0.;
    };

    jacobi (alpha, mfit, d, v, &nrot);  /* Find eigenvec and eigenval */

    /****************************************************************\
    *  This is the delta chi^2 contour that jointly bounds the 68%   *
    *  confidence region for all the mfit parameters, not just a     *
    *  single parameter.                                             *
    \****************************************************************/

/*    delchi2 = rootfind (&gammp, mfit, 0.8*mfit, mfit * 1.2, 0.05); */
    delchi2 = 1.;

    for (m=1; m<=mfit; m++) {
        sigma[m] = fabs(v[m][1]/sqrt(d[1]));
        for (n=1; n<=mfit; n++)
/*            sigma[m] += (fabs(v[m][n]) / sqrt(d[n])); */
            sigma[m] = FMAX(sigma[m], fabs(v[m][n]) / sqrt(d[n]));
    };

    /*****************************************************\
    *  Scale the uncertainties by the sqrt(delta chi^2).  *
    \*****************************************************/

    for (m=1; m<=mfit; m++)
        sigma[m] = sigma[m] * sqrt(delchi2);

    sortsig (sigma, ma, ia, mfit);
    free_dvector (d, 1, mfit);
    free_dmatrix(v, 1, mfit, 1, mfit);
}


#define SWAP(a,b) {swap=(a); (a)=(b); (b)=swap;}

void sortsig(double *sigma, int ma, int ia[], int mfit)
{
        int i,k;
        double swap;

        for (i=mfit+1; i<=ma; i++)
            sigma[i] = 0.;

        k = mfit;
        for (i=ma;i>=1;i--) {
            if (ia[i]) {
                SWAP (sigma[i], sigma[k]);
                k--;
	    };
        };
}

#undef SWAP


/************************************************************************/


#define JMAX 40
#define ONESIG 0.6827

float rootfind(float (*func)(float, float), float nu, float x1, float x2,
								float xacc) 
{
	void nrerror(char error_text[]);
	int j;
	float dx,f,fmid,xmid,rtb;

	f=(*func)(nu/2., x1/2.) - ONESIG;
	fmid=(*func)(nu/2., x2/2.) - ONESIG;
	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(nu/2., (xmid=rtb+(dx *= 0.5))/2.) - ONESIG;
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}
#undef JMAX
