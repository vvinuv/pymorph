#include <math.h>
#include "structs.h"
#include "nrutil.h"
#define JMAX 40

double LevMar (struct image *d, struct image *n, struct image *psf,
    struct fitpars *fpar, struct convpars *cpar, double **covar, int *ndof,
    int type);
void copy_structs (struct fitpars *new, struct fitpars *orig);
void delete_fitpars (struct fitpars *fpar);

float rtbis(double *chi2, float x1, float x2, float *parptr, 
      struct image *d, struct image *n, struct image *psf, struct
      fitpars *tptr, struct fitpars *fpar, struct convpars *cpar)
{
    void nrerror(char error_text[]);
    int j, ndof, flag, done;
    double dx,f,fmid,xmid,rtb;
    double chi2targ, yacc, **covar;
    float del;
    extern struct sampars sample;

    del = x2 - x1;
    done = 0;

    covar = dmatrix (1, sample.nobjs * NPARS, 1, sample.nobjs * NPARS);

    while (!done) {
	f = *chi2;

        *parptr = x2;
	fmid = LevMar (d, n, psf, tptr, cpar, covar, &ndof, 0);
        yacc = 0.1 * (*chi2)/ndof;

        if (fmid < *chi2) {
	    *chi2 = fmid;
	    copy_struct (fpar, tptr);
	};


        chi2targ = *chi2 + *chi2/ndof;

        f = chi2targ - f;
        fmid = chi2targ - fmid;

	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);

	flag = 0;
	for (j=1;j<=JMAX;j++) {
            copy_struct (tptr, fpar);
            xmid = *parptr = rtb+(dx *= 0.5);
	    fmid = LevMar (d, n, psf, tptr, cpar, covar, &ndof, 0);
            if (fmid < *chi2) {
	        *chi2 = fmid;
   	        copy_struct (fpar, tptr);
                chi2targ = *chi2 + *chi2/ndof;
	        x1 = xmid;
	        x2 = xmid + (dx *= 2.);
		flag = 1;
	    };
            if (flag == 1) break;
            fmid = chi2targ - fmid;

	    if (fmid <= 0.0) rtb=xmid;
	    if (fabs(fmid) <= yacc || fmid == 0.0) {
	        free_dmatrix (covar, 1, sample.nobjs * NPARS, 1, 
							sample.nobjs * NPARS);
	        return rtb;
	    };
	};
        if (j==JMAX)
	    done = 1;
    };
    nrerror ("Too many bisections in rtbis");
    return 0.0;
}
#undef JMAX
