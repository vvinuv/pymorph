#include <math.h>
#include "structs.h"
#include "nrutil.h"
#define MAXIT 30

double LevMar (struct image *d, struct image *n, struct image *psf,
    struct fitpars *fpar, struct convpars *cpar, double **covar, int *ndof,
    int type);
void copy_structs (struct fitpars *new, struct fitpars *orig);


float rtsec (double *chi2, float x1, float x2, float *parptr, struct image *d,
		struct image *n, struct image *psf, struct fitpars *tptr,
		struct fitpars *fpar, struct convpars *cpar)
{
    void nrerror(char error_text[]);
    int j, ndof, done, flag;
    double tf,fl,f,dx,swap,xl, chi2targ, yacc;
    extern struct sampars sample;
    double **covar, rts, del;

    covar = dmatrix (1, sample.nobjs * NPARS, 1, sample.nobjs * NPARS);
    del = (double)(x2 - x1);
    done = 0;

/*
    *parptr = x1;
    fl=LevMar (d, n, psf, tptr, cpar, covar, &ndof, 0);
    if (fl < *chi2) {
	*chi2 = fl;
	copy_struct (fpar, tptr);
    };

*/

    while (!done) {
        fl = *chi2;

        copy_struct (tptr, fpar); 
        *parptr = x2;
        f=LevMar (d, n, psf, tptr, cpar, covar, &ndof, 0);

        if (f < *chi2) {
   	    *chi2 = f;
	    copy_struct (fpar, tptr);
        };

        chi2targ = *chi2 + ndof; /* *(*chi2/ndof); */

        fl = chi2targ - fl;
        f = chi2targ - f;

        if (fabs(fl) < fabs(f)) {
 	    rts=x1;
	    xl=x2;
	    swap=fl;
	    fl=f;
	    f=swap;
        } else {
 	    xl=x1;
	    rts=x2;
        }
        flag = 0;
 
	yacc = 0.1 * ndof; /* *(*chi2/ndof); */
        for (j=1;j<=MAXIT;j++) {
            dx = (xl-rts)*f/(f-fl);
            if (dx > 0)
                if (dx > fabs(3 * del)) 
		    dx = fabs(3 * del);
            else if (dx < 0)
                if (fabs(dx) > fabs(3 * del))
                    dx = -fabs(3 * del);
	    xl=rts;
	    fl=f;
	    rts += dx;
            copy_struct (tptr, fpar);
            *parptr = rts;
            tf =  LevMar (d, n, psf, tptr, cpar, covar, &ndof, 0);
            if (tf < *chi2) {
                *chi2 = tf;
	        copy_struct (fpar, tptr);
                chi2targ = *chi2 + ndof; /* * (*chi2/ndof); */
		x1 = rts;
		x2 = x1 + del;
                flag = 1;
  	    };
            if (flag == 1) break;
  	    f= chi2targ - tf;

	    if (fabs(f) < yacc || dx <= 1.e-3) {
 	        free_dmatrix (covar, 1, sample.nobjs * NPARS, 1, 
						sample.nobjs * NPARS);
 	        return rts;
	    };
	};
        if (j == MAXIT) done = 1;
    };

   
    nrerror("Maximum number of iterations exceeded in rtsec");
    return 0.0;
}
#undef MAXIT

