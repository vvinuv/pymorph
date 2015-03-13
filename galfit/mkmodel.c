/***************************************************************\
*  This subroutine makes calls functions that make the models   *
*  and convolve the models with the PSF.                        *
*                                                               *
*  The model making happens in two steps.  First, models will   *
*  be created centered on the values specified in fpar->a[1]    *
*  and fpar->a[2] for all objects.  This model will not be      *
*  used in convolution because doing so would broaden out the   *
*  objects at the center.  While this is ok for galaxies it's   *
*  rather bad for fitting point sources.  Therefore, a second   *
*  model is created which is centered precisely on the pixel    *
*  center.                                                      *
*                                                               *
*  Note that regions inside "df->imgsect" will be convolved by  *
*  the PSF, but regions outside it will not.  First, to do      *
*  this correctly, one has to worry about how border padding    *
*  will influence the convolution within "df->imgsect".         *
*                                                               *
*  Specifically, since we're using FFT to do convolution, we    *
*  need to extend the image so that each side is 2^N in         *
*  length and width.  If we pad the extension with 0, then      *
*  around the edges of the model, convolution will effectively  *
*  convolve the model flux values with 0, hence corrupting the  *
*  model there.  But, instead of padding with 0, we pad a       *
*  square annulus around that region with real model data,      *
*  then only the region inside the annulus will be corrupted,   *
*  because the thickness of the annulus will be exactly half    *
*  the width and length of the convolution PSF.  The corrupted  *
*  region (annulus) will then be thrown away.                   *
*                                                               *
*  Therefore we devise the following scheme:                    *
*    1) The model inside the region "df->imgsect" will          *
*       be created in a separate array.                         *
*    2) That array will be extended on all sides by 1/2 the     *
*       length/width of the PSF, and that region will be        *
*       filled with model values.                               *
*    3) The rest of the image will be extended so that each     *
*       side has 2^N pixels in length and width.  This region   *
*       can be filled with 0 (or whatever, actually) padding.   *
*    4) Convolution takes place in this 2-D array.              *
*    5) Once convolution is done, copy this array back into     *
*       the final model image, being careful to leave out the   *
*       PSF extended annulus and the 0 padding.                 *
*                                                               *
\***************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "mymath.h"
#include "structs.h"
#include "nrutil.h"
#include "debug.h"

void create_conv_arrays (struct fitpars *fpar, struct convpars *cpar,
    struct image *model, struct derivs *df);
void convolve (struct image *model, double **fftpsf, long complex_naxes[]);
void annihilate_work_arrays (struct fitpars *fpar, struct image *, 
    struct derivs *df);
void copy_convregion (struct convpars *cpar, float **model,
    struct image convolved_img, int os);
double **psfprep (struct image psf, long model_naxes[], int add0pad,
    long complex_naxes[]);
void shift_psf (struct image psf, struct image *tpsf, float dx, float dy);
void objselect (float a[], int ia[], char objtype[], struct image *model, 
    struct derivs *df);
void clearmodel (struct image *model, struct derivs *df, struct fitpars *fpar);
void delete_derivs (struct derivs *df);

void mkmodel (struct image *model, struct derivs *df, struct fitpars *fpar, 
    struct convpars *cpar, int output) 
{

    void copy_model (struct fitpars *fpar, struct image *model, 
        struct derivs *df);
    extern struct image psf, tpsf, *kernel;
    extern struct inpars input;

    float ta[NPARS+1], dx, dy, xorig, yorig, xos, yos;
    unsigned long ndata;
    struct image tmodel;
    struct derivs tdf, *tptr, *dptr;
    struct fitpars *fptr, tfptr;
    double **fftpsf, **fftkern;
    int i, j, xmin, xmax, ymin, ymax, os;
    long complex_naxes[3];

    clearmodel (model, df, fpar);

    if (output)
        for (i=1; i<= NPARS; i++)
            tfptr.ia[i] = 0;

    /*****************************************************************\
    *  Create models for the entire fitting region.  The models       *
    *  created here are centered on the parameter centroids in        *
    *  fpar->a[1] and fpar->a[2], as opposed to exactly on the pixel  *
    *  center (below, in the next block of model generation).  This   *
    *  two step process is necessary so we don't broaden out the      *
    *  model profiles unnecessarily at the center-most pixel.         *
    \*****************************************************************/

    fptr = fpar;
    dptr = df;
    while (fptr != NULL) {
        if (output && fptr->outtype == 0)
            objselect (fptr->a, tfptr.ia, fptr->objtype, model, dptr);
        else if (!output)
            objselect (fptr->a, fptr->ia, fptr->objtype, model, dptr);
	fptr = fptr->next;
	dptr = dptr->next;
    };

    /*****************************************************************\
    *  Now create a temporary model image of *just* the convolution   *
    *  region to be used in convolution.  The temporary model is      *
    *  created exactly on pixel center, as opposed to the models      *
    *  that are not used in convolution (above), which are centered   *
    *  on actual object centers.  The fractional pixel shift in the   *
    *  convolved model comes about by convolving the temporary        *
    *  model with a shifted PSF.  The convolution region in the       *
    *  temporary model is then copied back to the unconvolved model.  *
    *								      *
    *  Notice the model used that's being convolved is created to     *
    *  be larger than the convolution box, by the size of the PSF,    *
    *  or one pixel larger if the size is an odd number.	      *
    \*****************************************************************/

    if ( strncmp (psf.name, "none", 4) != 0) {

        /**************************************************************\
        *  tmodel is just a dummy, pretty much.  The important models  *
        *  are stored in df->dp[0] matrices.                           *
        \**************************************************************/
    
        tmodel.dp[1] = model->dp[1];
        tmodel.dp[2] = model->dp[2];
        tmodel.muzpt = model->muzpt;
        tmodel.magzpt = model->magzpt;
        tmodel.naxes[1] = 0;
        tmodel.naxes[2] = 0;

        create_conv_arrays (fpar, cpar, model, &tdf);
        clearmodel (&tmodel, &tdf, fpar);

        /***************************************************\
        *  Cycle through all the objects and create models  *
        \***************************************************/

        fptr = fpar;
        tptr = &tdf;
	dptr = df;

	os = input.sampfac;

        while (fptr != NULL) {
            if (strncmp(fptr->objtype, "sky", 3)!=0 && 
		strncmp(fptr->objtype, "psf", 3)!=0) {

                sscanf (tptr->imgsect, "[%d:%d,%d:%d]", &xmin, &xmax, &ymin,
		    &ymax); 
	        strcpy (tmodel.imgsect, tptr->imgsect);

                tmodel.naxes[1] = tptr->naxes[1];
                tmodel.naxes[2] = tptr->naxes[2];

                for (i=1; i<= NPARS; i++) 
	            ta[i] = fptr->a[i];

		xorig = fptr->a[1] - xmin + 1;
		yorig = fptr->a[2] - ymin + 1;

	        ta[1] = NINT(xorig*os - 0.5*(os-1)) + NINT (cpar->psfsz[1]/2.);
                ta[2] = NINT(yorig*os - 0.5*(os-1)) + NINT (cpar->psfsz[2]/2.);
		ta[4] = ta[4] * os;

		/*************************************************\
		*  Special treatment for King and Nuker profiles  *
		\*************************************************/

		if (strncmp(fptr->objtype, "king", 4) == 0)
		    ta[5] = ta[5] * os;
		if (strncmp(fptr->objtype, "king", 4) == 0 || 
		    strncmp(fptr->objtype, "nuker", 5) == 0)
		    ta[3] = ta[3] + 5 * log10(os);

                if (output && fptr->outtype == 0 && tptr->naxes[1] != 0 && 
		    tptr->naxes[2] != 0)

                    objselect (ta, tfptr.ia, fptr->objtype, &tmodel, tptr);
                else if (!output && tptr->naxes[1] != 0 && tptr->naxes[2] != 0)
                    objselect (ta, fptr->ia, fptr->objtype, &tmodel, tptr);

            /*****************************************************\
            *  Now do convolution, cycle through all the objects  *
            \*****************************************************/

                for (i=0; i <= NPARS; i++) {
                    if (tmodel.naxes[1] != 0 && tmodel.naxes[2] != 0) {
                        if ((fptr->ia[i] == 1 && !output) || i == 0) {
	                    /*  Only need to shift the PSF once  */
	                    if (i == 0) { 
				xos = fptr->a[1] * os - 0.5*(os-1);
				yos = fptr->a[2] * os - 0.5*(os-1);
	                        dx = xos - NINT(xos);
	                        dy = yos - NINT(yos);

	                        shift_psf (psf, &tpsf, dx, dy);
			        fftpsf = psfprep (tpsf, tmodel.naxes, 0, 
					 complex_naxes);
			    };

			    if (dptr->naxes[1] !=0 && dptr->naxes[2] !=0) {
	                        /*  convolve.c takes struct image.  */
	                        tmodel.z = tptr->dpm[i]; 
                                convolve (&tmodel, fftpsf, complex_naxes); 

				/*  Charge diffusion kernel is also  *\
			        \*  applied in here.                 */
	                        copy_convregion (cpar, dptr->dpm[i], tmodel,
				    os);
			    };
			};
		    };
		}; 

                if (tmodel.naxes[1] > 0 && tmodel.naxes[2] > 0) 
                    free_dmatrix (fftpsf, 0, complex_naxes[2]-1, 0, 
			complex_naxes[1]*2-1);
	    };

	    fptr = fptr->next;
	    tptr = tptr->next;
	    dptr = dptr->next;
	};
    };

    /*****************************************************************\
    *  Now add all the subcomponents together into one net model      *
    *  image, merging together all the convolution pieces.  The       *
    *  derivative images are all taken care of by copy_convregion     *
    *  above.  The charge diffusion is also applied in there.         *
    \*****************************************************************/

    copy_model (fpar, model, df);

    if ( strncmp (psf.name, "none", 4) != 0) {
        tmodel.naxes[1] = 0;
        tmodel.naxes[2] = 0;

        annihilate_work_arrays (fpar, &tmodel, &tdf);
	delete_derivs (tdf.next);
    };

}


/*
    sprintf (tmodel.name, "tmodel.fits");
    writefits ("tmodel.fits", &tmodel, "tmodel", 0);
*/

/****************************************************************************/

void copy_model (struct fitpars *fpar, struct image *model, struct derivs *df)
{
    struct derivs *dptr;
    struct fitpars *fptr;
    int ix, iy, xmin, xmax, ymin, ymax;

    fptr = fpar;
    dptr = df;

    while (fptr != NULL) {
	for (iy = 1; iy <= model->naxes[2]; iy++)
	    for (ix = 1; ix <= model->naxes[1]; ix++)
		    model->z[iy][ix] += dptr->dpm[0][iy][ix];

        fptr = fptr->next;
	dptr = dptr->next;
    };

}
