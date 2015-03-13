#include <stdio.h>
#include <math.h>
#include <curses.h>
#include <string.h>
#include "structs.h"
#include "nrutil.h"

float rtbis (double *chi2, float xmin, float xmax, float *parptr, 
	        struct image *d, struct image *n, struct image *psf, struct
		fitpars *tptr, struct fitpars *fpar, struct convpars *cpar);

float rtsec (double *chi2, float xmin, float xmax, float *parptr, 
	        struct image *d, struct image *n, struct image *psf, struct
		fitpars *tptr, struct fitpars *fpar, struct convpars *cpar);

float parabolic (double *chi2, float xmid, float del, float *parptr, 
	        struct image *d, struct image *n, struct image *psf, struct
		fitpars *tptr, struct fitpars *fpar, struct convpars *cpar);
void delete_fitpars (struct fitpars *fpar);
void copy_struct (struct fitpars *new, struct fitpars *fpar);

int geterr (double *chi2, int ndof, struct image *d, struct image *n, 
	struct image *psf, struct fitpars *fpar, struct convpars *cpar, 
	float *uncert, double *sigma)
{

    float getdel (int i, struct fitpars *fpar, double cov);

    extern struct sampars sample;
    struct fitpars *tptr, *new, *orig;
    int errcond, newchi2, i, k, j=1;
    double origchi2;
    float del;

    orig = fpar;

    /*******************************************************************\
    *  Make a copy of the entire fitpar structure.  The parameters in   *
    *  this structure will be tweaked in order to determine the         *
    *  uncertainties.                                                   *
    \*******************************************************************/

    tptr = (struct fitpars *) malloc (sizeof (struct fitpars));
    new = tptr;

    while (orig != NULL) {
        if (orig->next != NULL) {
            tptr->next = (struct fitpars *) malloc (sizeof (struct fitpars));
            tptr = tptr->next;
	} else
	    tptr->next = NULL;
        orig = orig->next;
    };

    tptr = new;
    orig = fpar;

    /**********************************************************************\
    *  Now iterate through all the parameters and traverse the chi^2       *
    *  surface by holding one parameter fixed and allowing others to vary  *
    *  until chi^2 = chi^2best + chi^2best/NDOF is reached.                *
    \**********************************************************************/

    origchi2 = *chi2;
    while (orig != NULL) {
        copy_struct (new, fpar);
        for (i=1; i <= NPARS; i++) {
            if (orig->ia[i] == 1) {
		tptr->ia[i] = 0;
		orig->ia[i] = 0;
/*		del = getdel (i, tptr, covar[i][i]); */
		del = sigma[j];
                uncert[j] = parabolic (chi2, tptr->a[i], del,
				&(tptr->a[i]), d, n, psf, new,
						fpar, cpar);
		tptr->ia[i] = 1;
		orig->ia[i] = 1;
            } else
   	        uncert[j] = 0.;
            j++;
        };
	tptr = tptr->next;
        orig = orig->next;
    };

    j--;
    delete_fitpars (new->next);

/*    if (*chi2/ndof <= 0.99 * origchi2/ndof )
        return (0);
    else
        return (1);  */
    return (1);
}


float getdel (int i, struct fitpars *fpar, double c)
{

    if (i == 1 || i == 2) 
        return (sqrt(c));

    if (i == 3)
	return (sqrt(c));

    if (i >= 4 && i <= 8)
        return (sqrt(c));

    if (i == 9)
	return (sqrt(c));

    if (i == 10.)
	return (sqrt(c));
}
