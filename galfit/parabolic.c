#include <math.h>
#include "structs.h"
#include "nrutil.h"
#define MAXIT 30

double LevMar (struct image *d, struct image *n, struct image *psf,
    struct fitpars *fpar, struct convpars *cpar, double **covar, 
    double *sig, int *ndof, int type);
void copy_structs (struct fitpars *new, struct fitpars *orig);


float parabolic (double *chi2, float xmid, float del, float *parptr, struct
		image *d, struct image *n, struct image *psf, struct fitpars
		*tptr, struct fitpars *fpar, struct convpars *cpar)
{
    int j, ndof, done, skip;
    extern struct sampars sample;
    double **covar, f, chisq2, chisq1, chisq3, nchi2;
    double *dummy;
    double sigma;

    covar = dmatrix (1, sample.nobjs * NPARS, 1, sample.nobjs * NPARS);
    dummy = dvector (1, sample.nobjs * NPARS);

    done = 0;
    skip = 0;
    while (!done) {
        if (skip != 2) {
            copy_struct (tptr, fpar); 
            *parptr = xmid + del;
            nchi2=LevMar (d, n, psf, tptr, cpar, covar, dummy, &ndof, 0);
            if (nchi2 < *chi2) {
  	        xmid = *parptr;
                copy_struct (fpar, tptr); 
	        chisq3 = *chi2;   /*  /ndof;  */
                *chi2 = nchi2;   
                done = 0;
                skip = 1;
            } else {
                chisq1 = nchi2;   /*  /ndof;  */
                skip = 0;
	    };
        };

        if (skip != 1) {
            copy_struct (tptr, fpar); 
            *parptr = xmid - del;
            nchi2=LevMar (d, n, psf, tptr, cpar, covar, dummy, &ndof, 0);
            if (nchi2 < *chi2) {
		xmid = *parptr; 
                copy_struct (fpar, tptr); 
                chisq1 = *chi2;     /*  /ndof;  */ 
                *chi2 = nchi2; 
                done = 0; 
                skip = 2;
	    } else {
                chisq3 = nchi2;     /*  /ndof;  */
                skip = 0;
	    };
        };

        chisq2 = *chi2;    /*  /ndof;   */

        if (skip == 0) {
            copy_struct (tptr, fpar); 
            *parptr = (xmid-del) + del * ((chisq3 - chisq2)/
				(chisq1 - 2 * chisq2 + chisq3) + 0.5);

            if (*parptr > xmid+del || xmid-del > *parptr)
                nchi2 = LevMar (d, n, psf, tptr, cpar, covar, dummy, &ndof, 0);

            if (nchi2 < *chi2) {
                *chi2 = nchi2;
                copy_struct (fpar, tptr); 
                done = 0;
            } else
                done = 1;
        };
    };
    free_dvector (dummy, 1, sample.nobjs * NPARS);
    free_dmatrix (covar, 1,  sample.nobjs * NPARS, 1, sample.nobjs * NPARS);
    sigma = del * sqrt (2 / (chisq1 - 2 * chisq2 + chisq3));
    return (sigma);
}
#undef MAXIT

