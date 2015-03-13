#include <stdio.h>
#include <string.h>
#include <curses.h>
#include "structs.h"
#include "nrutil.h"

#define DONE 1

struct image *mask;      /* stamp-sized region of bad pixel mask */
struct image model, cmodel, tpsf, sincwin;
struct derivs df;   /* spaces reserved for derivative images */
struct inpars input;
float **y1a, **y2a, **y12a;
float xskycent, yskycent;
struct bicubic_grid **cg;

struct image *sigma_img (struct image *data);
double LevMar (struct image *d, struct image *n, struct image *psf, 
    struct fitpars *fpar, struct cons *constr, struct convpars *cpar, 
    double **covar, double *sigma, int *ndof, int type);
struct bicubic_grid **get_cgrid (struct image *psf);
void fprintpar (int lowx, int lowy, struct fitpars *par, FILE *ftype, 
    char newout[], double *sigma, float *coerr, double chisq, int ndof);
void annihilate_work_arrays (struct fitpars *fpar, struct image *, 
    struct derivs *);
void create_work_arrays(struct fitpars *fpar, struct image *, struct derivs *);
int geterr (double *chi2, int ndof, struct image *d, struct image *n, struct
    image *psf, struct fitpars *fpar, struct convpars *cpar, float
    *coerr, double *sigma);
void outmenu (struct fitpars *fpar, char newout[], float xoffset, 
    float yoffset, float chisq, int ndof);
void outmodel (struct fitpars *fpar, struct convpars *cpar, struct image *d, 
    struct image *model, struct derivs *df, double chisq, int ndof,
    int offx, int offy);
double chi2calc (struct fitpars *fpar, float **data, float **sigma, 
    float **model, long naxes[], int *ndof);
int (*pfunc) (const char *, ...);
void free_cgrid (struct bicubic_grid **m, long nrl, long nrh, long ncl, 
    long nch);
void assign_err (struct fitpars *fpar, double *sigma, float chi2nu);

void search_min (struct image *data, struct image *sigma, struct image *psf,
    struct image *badpix, struct fitpars *fpar, struct cons *constr,
    struct convpars *cpar) 
{
    void offset_center (struct fitpars *fpar, struct convpars *cpar, int xoff,
							int yoff);
    struct image *to_mini (struct image *);

    FILE *logfile;
    struct image *d, *s;  /* sub-region of data and noise image to fit */
    int i, ndof, lowx, lowy, highx, highy, optimized;
    double **covar, *sig;
    float *coerr;
    double chisq;
    long tnaxes[3];
    char newout[11];
    struct fitpars *fptr;
    extern struct sampars sample;
    extern char device[];

    /* When doing more than one object in an image, can modify this *\
    \* subroutine to loop through all the objects                   */

 /************************************************************************\
 *                                                                        *
 *    Set up all the image matrices for fitting once and for all          *
 *                                                                        *
 \************************************************************************/

    /*-----------------------------------------------------------------\
    |  First cut the section of the image that we want to fit from     |
    |  the original data/sigma images into new "mini-images."  In      |
    |  general the mini-images will be larger than the convolution     |
    |  region.  The mini-images are regions the user specified         |
    |  to fit with no convolution padding added.                       |
    \-----------------------------------------------------------------*/

    d = to_mini (data);
    if (sigma->err == 0)
        s = to_mini (sigma);
    else
	s = sigma_img (d);
    mask = to_mini (badpix);

    /*-------------------------------------------------\
    |  Offset the object centers and the convolution   |
    |  box relative to the stamp-sized image           |
    \-------------------------------------------------*/

    sscanf (d->imgsect, "[%d:%d,%d:%d]", &lowx, &highx, &lowy, &highy);
    offset_center (fpar, cpar, -(lowx-1), -(lowy-1));

    xskycent = (d->naxes[1] + 1)/ 2.;
    yskycent = (d->naxes[2] + 1)/ 2.;

    /*------------------------------------------\
    |  Prepare for bicubic interpolation of PSF |
    \------------------------------------------*/

    if (strncmp (psf->name, "none", 4) != 0) { 
        tpsf.naxes[1] = psf->naxes[1];
        tpsf.naxes[2] = psf->naxes[2];

	tpsf.z = matrix (1, tpsf.naxes[2], 1, tpsf.naxes[1]);
	cg = get_cgrid (psf);
    }

    /*------------------------------------------------------------\
    |  Create storage space for the derivative and model images.  |
    |  Store them as global.  These are real (as opposed to       |
    |  complex).  They're the same size as the postage sized      |
    |  image and will NOT be 2^N in size and will not have any    |
    |  padding.  "model" and "df" will have the size of the       |
    |  entire fitting region.                                     |
    \------------------------------------------------------------*/

    model = *d;
    create_work_arrays (fpar, &model, &df);

 /************************************************\
 *  Create covariance matrix and then minimize!   *
 \************************************************/

    covar = dmatrix (1, sample.nobjs * NPARS, 1, sample.nobjs * NPARS);
    coerr = vector (1, sample.nobjs * NPARS);
    sig = dvector (1, sample.nobjs * NPARS);

    optimized = 0;
    while (!optimized && !input.create) {
        chisq = LevMar (d, s, psf, fpar, constr, cpar, covar, sig, &ndof, 1);

        /*------------------------------------------------------\
        |  Now obtain uncertainties using the parabolic method  |
        \------------------------------------------------------*/

/*        if (input.coerr)
            optimized = geterr (&chisq, ndof, d, s, psf, fpar, cpar, coerr, 
								sig); 
        else   */
	    optimized = 1;
    };

    assign_err (fpar, sig, chisq/ndof);

 /*********\
 *  DONE!   \*************************************************\
 *  Now create model, output the data, model, residual,       *
 *  images into a 3-layered block.                            *
 \************************************************************/

    pfunc ("\n\nFit summary is now being saved into `fit.log'.\n\n");
    if (strncmp (device, "curses", 6) == 0)
        refresh();

    outmodel (fpar, cpar, d, &model, &df, chisq, ndof, (lowx-1), (lowy-1));

    /*  if only creating an output image but not fitting  */

    if (input.create)
        chisq = chi2calc (fpar, d->z, s->z, model.z, model.naxes, &ndof);

 /*********************************************************\
 *  Restore x and y offsets, then output the initial and   *
 *  final parameters into files.                           *
 \*********************************************************/

    xskycent += (lowx-1);
    yskycent += (lowy-1);
    offset_center (fpar, cpar, (lowx-1), (lowy-1));
    outmenu (fpar, newout, 0., 0., chisq, ndof);

    if ( (logfile = fopen ("fit.log", "a")) == (FILE *) 0)
        logfile = fopen ("fit.log", "w");
    fprintpar (1, 1, fpar, logfile, newout, sig, coerr, chisq, ndof);
    fclose (logfile);

 /***********************************************\
 *  Now free up all the memory grabbed up.       *
 \***********************************************/

    annihilate_work_arrays (fpar, &model, &df);
    free_dmatrix (covar, 1, sample.nobjs * NPARS, 1, sample.nobjs * NPARS);
    free_matrix (mask->z, 1, mask->naxes[2], 1, mask->naxes[1]);
    free_matrix (d->z, 1, d->naxes[2], 1, d->naxes[1]);
    free_matrix (s->z, 1, s->naxes[2], 1, s->naxes[1]);
    free_dvector (sig, 1, sample.nobjs * NPARS);
    if (strncmp (psf->name, "none", 4) != 0)
        free_cgrid (cg, 0, psf->naxes[2]+1, 0, psf->naxes[1]+1);

    free (mask);
    free (s);
    free (d);

    return;
}


/*****************************************************************************/


struct image *to_mini (struct image *img)
{
    extern int (*pfunc)(const char *, ...);
    extern char device[];
    unsigned long i, j, xmin, xmax, ymin, ymax;
    struct image *new;

    new = (struct image *) calloc (1, sizeof (struct image));
    if (new == (struct image *) 0) {
	pfunc ("Error in allocating memory for mini-image....");
        if (strncmp(device, "curses", 6)== 0)
            refresh();
	exit (1);
    };

    *new = *img;
    sscanf (img->imgsect, "[%d:%d,%d:%d]", &xmin, &xmax, &ymin, &ymax);
    new->naxes[1] = xmax - xmin + 1;
    new->naxes[2] = ymax - ymin + 1;

    new->z = matrix(1, new->naxes[2], 1, new->naxes[1]);

    xmax = IMIN (xmax, img->naxes[1]);   /*  In case the image size is      */
    ymax = IMIN (ymax, img->naxes[2]);   /*  smaller than the fitting size  */
					 /*  for whatever reason....        */

    for (j=ymin; j <= ymax; ++j){
	for (i=xmin; i <= xmax; ++i)
	    new->z[j-ymin+1][i-xmin+1] = img->z[j][i];
    }; 
    return (new);
}


/*****************************************************************************/

void offset_center (struct fitpars *fpar, struct convpars *cpar, int xoff,
							int yoff)
{
    while (fpar != NULL) {
        if (strncmp (fpar->objtype, "sky", 3) != 0) {
            fpar->a[1] += xoff;
            fpar->a[2] += yoff;
        };
        fpar = fpar->next;
    };
    cpar->cent[1] += xoff;
    cpar->cent[2] += yoff;
}

