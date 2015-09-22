#include <math.h>
#include "mymath.h"
#include "structs.h"
#include "nrutil.h"

#define SHIFT 0.01

void shift_psf (struct image psf, struct image *tpsf, float dx, float dy);

void psfunc (float a[], int ia[], struct image *model, struct derivs *df)
{
    extern struct image psf, tpsf;
    struct image xp_psf, xm_psf, yp_psf, ym_psf;
    int i, ix, iy, xlow, ylow, x, y;
    float ta[NPARS+1], dfdI, dIdmag, flux, dx, dy;
    static double norm = 0.;
    double totcnts;

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

    dx = a[1] - NINT(a[1]);
    dy = a[2] - NINT(a[2]);

    shift_psf (psf, &tpsf, dx, dy);

    xp_psf.naxes[1] = xm_psf.naxes[1] = yp_psf.naxes[1] = ym_psf.naxes[1] = tpsf.naxes[1];
    xp_psf.naxes[2] = xm_psf.naxes[2] = yp_psf.naxes[2] = ym_psf.naxes[2] = tpsf.naxes[2];
    if (ia[1] == 1) {
        xp_psf.z = matrix (1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        xm_psf.z = matrix (1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        shift_psf (psf, &xp_psf, SHIFT, 0.);      /*  dx */
        shift_psf (psf, &xm_psf, -SHIFT, 0.);     /* -dx */
    };

    if (ia[2] == 1) {
        yp_psf.z = matrix (1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        ym_psf.z = matrix (1, tpsf.naxes[2], 1, tpsf.naxes[1]);
        shift_psf (psf, &yp_psf, 0, SHIFT);       /*  dy */
        shift_psf (psf, &ym_psf, 0, -SHIFT);      /* -dy */
    };

    xlow = NINT(a[1]) - (int)(psf.naxes[1]/2.);
    ylow = NINT(a[2]) - (int)(psf.naxes[2]/2.);
    for (iy=1; iy <= psf.naxes[2]; iy++) {
        for (ix=1; ix <= psf.naxes[1]; ix++) {
	    x = ix + xlow - 1;
	    y = iy + ylow - 1;

	    /*  Make sure the shifted PSF is entirely within the image  */

	    if (x >= 1 && x <= df->naxes[1] && y >= 1 && 
	        y <= df->naxes[2]) {

	        dfdI = tpsf.z[iy][ix];

            	if (ia[1] == 1) {
                    df->dpm[1][y][x] = ta[3] * (xp_psf.z[iy][ix] - 
				  xm_psf.z[iy][ix])/(2. * SHIFT);
            	};
	        
        	if (ia[2] == 1) {
                    df->dpm[2][y][x] = -ta[3] * (yp_psf.z[iy][ix] - 
				  ym_psf.z[iy][ix])/(2. * SHIFT);
                };
          
                if (ia[3] == 1)
                    df->dpm[3][y][x] = dfdI * dIdmag;

                flux = ta[3] * tpsf.z[iy][ix];
		df->dpm[0][y][x] = flux;             /*  shifted PSF  */
	    };
	};
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

