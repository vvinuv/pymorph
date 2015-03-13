#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <curses.h>
#include "nrutil.h"
#include "structs.h"
#include "debug.h"

#define MAX 100
#define QUICKEND 10
#define HISTORY 10
#define MAXALAMDA 1.e10
#define MINALAMDA 1.e-10

void printpar (int lowx, int lowy, struct fitpars *par);
void writefits(char *phead, struct image *img, char *object, int add);
void funcs(int, int, float *, float [], int [], int, struct fitpars *fpar);
int keypoll (struct fitpars *fpar, int *ITERMAX, int lowx, int lowy, 
    float chisq, int ndof);
int mrqmin(float **model, float **sig, struct image *psf, int ma, 
    double **covar, double **alpha, double *chisq, 
    void (*funcs)(int, int, float *, float [], int [],
    int, struct fitpars *), float *alamda, struct fitpars *fpar, 
    struct cons *constr, struct convpars *cpar, double *sigma);
double chi2calc (struct fitpars *fpar, float **data, float **sigma,
    float **model, long naxes[], int *ndof);

struct inpars input;

double LevMar (struct image *d, struct image *s, struct image *psf, 
    struct fitpars *fpar, struct cons *constr, struct convpars *cpar, 
    double **covar, double *sig, int *ndof, int type)
{

    extern int (*pfunc) (const char *, ...);
    extern char device[];

    float oldchisq, alamda, dalamda, palamda, chi2nu, delchi2, avghist,
        olddelchi2;
    double **alpha;
    int countdown, niter, errcond, cflag, nfixed=0, i, totpars, nfree=0;
    unsigned long ndata;
    int lowx, highx, lowy, highy, quit, ITERMAX, hist, imax;
    struct fitpars *ptr;
    double dlim, chisq, history[HISTORY+1];

    extern struct sampars sample;
    extern struct image model;
    extern struct derivs df;

    sscanf (d->imgsect, "[%d:%d,%d:%d]", &lowx, &highx, &lowy, &highy);

    /******************************************************\
    *  Figure out how many free parameters are held fixed  *
    *  for statistical calculations                        *
    \******************************************************/
    ptr = fpar;
    while (ptr != NULL) {
        for (i=1; i <= NPARS; ++i) {
	    if ( ptr->ia[i] == 0 ) nfixed++;
            if ( ptr->ia[i] == 1 ) nfree++;
	    if (ptr->ia[i] == -1 ) ptr->ia[i] = 0;
	};
	ptr = ptr->next;
    };

    if (type == 1) {
        pfunc ("\n================================================================================\n");
        pfunc ("Initial parameters:\n");
        printpar (lowx, lowy, fpar);
        pfunc ("================================================================================\n");
        dlim = 1.e-4;
    }

    if (input.create) {
        pfunc ("\n\n");
        pfunc (">>> INPUT DATA IMAGE NOT FOUND.  Just creating the model specified.  <<< \n"); 
        pfunc (">>> The exposure time of the output image is assumed to be 1 second. <<<");
    };

    if (strncmp (device, "curses", 6) == 0)
        refresh();

    totpars = sample.nobjs * NPARS;

    /********************************************************************\
    * Set up various parameters used to control the exit of the fitting. *
    \********************************************************************/

    alpha = dmatrix (1, totpars, 1, totpars);
    chisq = 1.e27;
    ITERMAX = MAX;
    countdown = ITERMAX;
    niter = 0;
    cflag = 0;
    alamda = -1.;
    dalamda = 0.;
    palamda = 0.;
    quit = 0;

    /************************************************************\
    *                                                            *
    * Iterate to find the minimum until exit conditions are met. *
    *                                                            *
    \************************************************************/

    ndata = d->naxes[1] * d->naxes[2];

    if (nfree == 0) input.create = 2;
    hist = 1;
    while (1 && !input.create) {

        /**********************************\
        * Checking for exiting conditions  *
        \**********************************/

        if (countdown == -1 || nfree == 0 || quit)
            break;

	niter++;
	oldchisq = chisq;

	errcond = mrqmin (d->z, s->z, psf, totpars, covar, alpha, 
			&chisq, &funcs, &alamda, fpar, constr, cpar, sig);

	if (errcond) 
	    exit (errcond);

        *ndof = ndata - sample.nmask - nfree;
	chi2nu = chisq / *ndof;
	delchi2 = oldchisq - chisq;

	dalamda = alamda - palamda;
	palamda = alamda;

	/********************************************************\
        * Keep a history of delta chi^2 to see if the trend over *
	* HISTORY iterations is small enough to quit.            *
	\********************************************************/

	if (hist == HISTORY) 
	    hist = 1;
	if (niter != 1)           /* First iteration doesn't have a delchi2 */
	    history[hist++] = delchi2/sqrt(chi2nu);

	if (niter <= HISTORY)
	    imax = niter-1;
	else 
	    imax = HISTORY;
	
	avghist = 0.;
	for (i = 1; i <= imax; i++) 
	    avghist += (history[i] / imax);

	/******************************************************************\
	* Figure out when to exit the subroutine.  Countdown starts when:  *
        *    1.  Delta chi^2 (normalized by sqrt(chi2/nu) over history of  *
        *        10 iterations is less than 5, or                          *
        *    2.  delta Chi2/nu is less than dlim.                          *
        *                                                                  *
        * Countdown aborts when:					   *
        *    1.  The above are no longer true.                             *
        *    2.  When the delta chi^2 starts to be more negative, i.e.     *
        *        chi^2 is decreasing at a more rapid rate.                 *
	\******************************************************************/

	if (cflag == 0 && (avghist <= 5. || delchi2/chisq < dlim) &&
	    (delchi2 < olddelchi2)) {

	    countdown = QUICKEND;          /* Start early quit countdown */
	    cflag = 1;
	} else if (((delchi2 / sqrt(chi2nu) >= 5. && avghist >= 5.) ||
	    delchi2/chisq > dlim) || (delchi2 > olddelchi2)){

	    cflag = 0;
	    countdown = ITERMAX - niter;   /* Abort early quit countdown */
	};

	if (delchi2/chisq >= 1e-5)
            olddelchi2 = delchi2;
   
        /***************************************\
        * Display the current status of the fit *
        \***************************************/

	pfunc ("\n");
	pfunc ("Iteration : %-2d    Chi2nu: %-11.3e   ", niter, chi2nu);
	pfunc ("dChi2/Chi2: %-10.2e  alamda: %-10.0e\n", 
				(chisq - oldchisq)/chisq, alamda);
        printpar (lowx, lowy, fpar);
	pfunc ("COUNTDOWN = %d \n", countdown);

	countdown--;

        /***********************************\
        * Poll the keyboard for user input  *
        \***********************************/

        if (strncmp (device, "regular", 7) != 0)
            quit = keypoll(fpar, &ITERMAX, lowx, lowy, chisq, *ndof);
	if (ITERMAX != MAX && ITERMAX == niter)
	    quit = 1;

    };

    alamda = 0.;
    if (nfree > 0 && !input.create)
        mrqmin (d->z, s->z, psf, totpars, covar, alpha, &chisq, &funcs, 
	    &alamda, fpar, constr, cpar, sig);
    else
        chisq = chi2calc (fpar, d->z, s->z, model.z, model.naxes, ndof);

    /****************************************************\
    *  DONE!  Now announce what the exit conditions are  *
    *  if the user decided to quit early or if there is  *
    *  a problem with the input file.                    *
    \****************************************************/

    if (nfree == 0) {
        pfunc ("\n");
        pfunc ("-- No free parameters to fit.  Creating images only.\n\n");
    };

    if ( niter >= ITERMAX ) {
        pfunc ("\n\n The fitting ended because the maximum number of iterations \n ");
        pfunc ("has been reached, not necessarily because of convergence!\n\n");
    };

    free_dmatrix(alpha, 1, totpars, 1, totpars);

    return (chisq);
}


