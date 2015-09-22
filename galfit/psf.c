#include <math.h>
#include "mymath.h"
#include "structs.h"
#include "nrutil.h"

#define SHIFT 0.1

void shift_psf (struct image psf, struct image *tpsf, float dx, float dy);
double **psfprep (struct image psf, long model_naxes[], int add0pad,
    long complex_naxes[]);
void convolve (struct image *model, double **fftpsf, long complex_naxes[]);

void psfunc (float a[], int ia[], struct image *model, struct derivs *df)
{

    extern struct image psf, tpsf, *kernel;
    extern struct inpars input;
    struct image xp_psf, xm_psf, yp_psf, ym_psf, tmodel;
    int i, ix, iy, x, y, xx, yy, os;
    float ta[NPARS+1], dfdI, dIdmag, flux, dx, dy, xpsflow, ypsflow;
    long kern_complex_naxes[3];
    static double norm = 0.;
    double totcnts, **fftkern;

    os = input.sampfac;

    if (norm == 0.) {
        for (iy=1; iy <= psf.naxes[2]; iy++) {
	    for (ix=1; ix <= psf.naxes[1]; ix++)
	      norm += psf.z[iy][ix];
        };
    };

    for (i=1; i <= 3; i++)
        ta[i] = a[i];

    totcnts = pow (10., (model->magzpt - a[3])/2.5);
    ta[3] = totcnts / norm;

    dIdmag = ta[3] * log(10.) / -2.5;

    /*****************\
    *  Shift the PSF  *
    \*****************/

    ta[1] = a[1] * os - 0.5 * (os-1);
    ta[2] = a[2] * os - 0.5 * (os-1);

    dx = ta[1] - NINT(ta[1]);
    dy = ta[2] - NINT(ta[2]);

    shift_psf (psf, &tpsf, dx, dy);

    xp_psf.naxes[1] = xm_psf.naxes[1] = yp_psf.naxes[1] = ym_psf.naxes[1] = tpsf.naxes[1];
    xp_psf.naxes[2] = xm_psf.naxes[2] = yp_psf.naxes[2] = ym_psf.naxes[2] = tpsf.naxes[2];
    if (ia[1] == 1) {
        xp_psf.z = matrix (1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        xm_psf.z = matrix (1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        shift_psf (psf, &xp_psf, dx+SHIFT, 0.);      /*  dx */
        shift_psf (psf, &xm_psf, dx-SHIFT, 0.);      /* -dx */
    };

    if (ia[2] == 1) {
        yp_psf.z = matrix (1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        ym_psf.z = matrix (1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        shift_psf (psf, &yp_psf, 0, dy+SHIFT);       /*  dy */
        shift_psf (psf, &ym_psf, 0, dy-SHIFT);       /* -dy */
    };

    xpsflow = ta[1] - psf.naxes[1]/2;
    ypsflow = ta[2] - psf.naxes[2]/2; 

    for (iy=1; iy < psf.naxes[2]; iy++) {
        for (ix=1; ix < psf.naxes[1]; ix++) {
	    x= NINT((ix + xpsflow - 1 + 0.5 * (os-1))/os);
	    y= NINT((iy + ypsflow - 1 + 0.5 * (os-1))/os);

	    /*  Make sure the shifted PSF fits within the image  */

	    if (x >= 1 && x <= df->naxes[1] && y>= 1 && y <= df->naxes[2]) {

	        dfdI = tpsf.z[iy][ix];

           	if (ia[1] == 1) {
                    df->dpm[1][y][x] += ta[3] * (xp_psf.z[iy][ix] - 
		        xm_psf.z[iy][ix])/(2. * SHIFT/os);
            	};
	        
        	if (ia[2] == 1) {
                    df->dpm[2][y][x] += ta[3] * (yp_psf.z[iy][ix] - 
			ym_psf.z[iy][ix])/(2. * SHIFT/os);
                };
          
                if (ia[3] == 1)
                    df->dpm[3][y][x] += dfdI * dIdmag;

                flux = ta[3] * tpsf.z[iy][ix];
		df->dpm[0][y][x] += flux;          /*  shifted PSF  */
	    };
	};
    };

    if (kernel != NULL) {
        fftkern = psfprep (*kernel, df->naxes, 1, kern_complex_naxes);
	tmodel.naxes[1] = df->naxes[1];
	tmodel.naxes[2] = df->naxes[2];

        tmodel.z = df->dpm[0]; 
        convolve (&tmodel, fftkern, kern_complex_naxes);

	if (ia[1] == 1) {
	    tmodel.z = df->dpm[1];
	    convolve (&tmodel, fftkern, kern_complex_naxes);
	};

	if (ia[2] == 1) {
	    tmodel.z = df->dpm[2];
	    convolve (&tmodel, fftkern, kern_complex_naxes);
	};

	if (ia[3] == 1) {
	    tmodel.z = df->dpm[3];
	    convolve (&tmodel, fftkern, kern_complex_naxes);
	};

        free_dmatrix (fftkern, 0, kern_complex_naxes[2]-1, 0, 
	     kern_complex_naxes[1]*2-1);                        
    };

    if (ia[1] == 1) {
        free_matrix (xp_psf.z, 1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        free_matrix (xm_psf.z, 1, tpsf.naxes[2], 1, tpsf.naxes[1]);
    }

    if (ia[2] == 1) {
        free_matrix (yp_psf.z, 1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        free_matrix (ym_psf.z, 1, tpsf.naxes[2], 1, tpsf.naxes[1]);
    };
}

