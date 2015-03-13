#include <math.h>
#include <stdlib.h>
#include "structs.h"
#include "mymath.h"
#include "nrutil.h"
#include "const.h"
#include "debug.h"

#define OVERSAMPR 10.
#define MINR 1.e-15

float ratio (float c);
float dfridr (float (*ratio)(float n), float x, float h, float *err);

void exponential (float a[], int ia[], struct image *model, struct derivs *df)
{
    extern struct inpars input;

    unsigned long i, ix, iy;
    int j, nsubsamp, os;
    float xd, yd, gamma, delta, am, ta[NPARS+1], cosPA, sumdfdI0,
          sinPA, xf, yf, xc, yc, flux, ta4sqr, delota8,
	  gammaprime, gexp, dfdI0, dfdr, deltaprime, c, absgammac, 
          absdeltac, fabsrc, area_ratio, theta, step, drsdmag,
          totcnts, dc=0.1, dratio_dc,  err,
	  dfdn, dfdre, dfdq, dfdc, dcdmag, dndmag, dredmag, dqdmag;
    float dI0dmag, dI0dc, dI0drs, dI0dq;

    os = input.sampfac;

    for (i=1; i <= NPARS; i++)
        ta[i] = a[i];

    c = ta[10] + 2.;

    /*****************************************************************\
    *                                                                 *
    *  Pre-calculate some things we will frequently use later on but  *
    *  which needs to be calculated one time only.                    *
    *                                                                 *
    \*****************************************************************/

    /* a[3] is the galaxy total brightness while ta[3] */
    /* is the count rate at the disk scalelength.      */

    /* Notice the last term ta[8] in the first equation here.      */
    /* That's because the distance from the center is defined      */
    /* as  r = (|x|^c + |y/q|^c)^(1/c), where q is the axis ratio. */
    /* This would make the magnitude come out right.               */

    /* The area_ratio is the ratio of a pure ellipse area          */
    /* to that of a disky or boxy area.  Since we can calculate    */
    /* the flux of an elliptical model, the flux of a disky/boxy   */
    /* model is simply the ratio of the two areas.                 */

    area_ratio = ratio(c);
    dratio_dc = dfridr (ratio, c, dc, &err);

    ta4sqr = ta[4] * ta[4];
    totcnts = pow (10., (model->magzpt - a[3])/2.5);
    ta[3] = totcnts / 2. / PI / (ta4sqr * ta[8] / area_ratio);
    ta[9] = -1. * ta[9] / 180. * PI + PI/2.;    /* "Intuitive" orientation */
    cosPA = cos (ta[9]);
    sinPA = sin (ta[9]);

    /**********************************************************\
    *  Calculate some of the cross coupling terms between I0   *
    *  with all the other parameters.                          *
    \**********************************************************/

    dI0dmag = ta[3] * log(10.) / -2.5;

    dI0dq = -ta[3] / ta[8]; 
    dI0drs = -2 * ta[3] / ta[4]; 
    dI0dc = ta[3] / area_ratio * dratio_dc;

/* printf ("dI0dmag=%f   dI0dq=%f   dI0drs=%f   dI0dc=%f\n", dI0dmag, dI0dq, dI0drs, dI0dc); */

    /*******************************************\
    *                                           *
    *  Now fill the image and derivative grids  *
    *                                           *
    \*******************************************/

    for (iy=1; iy <= df->naxes[2]; iy++) {
        for (ix=1; ix <= df->naxes[1]; ix++) {

	    /* calculate the elliptical center distance */
	    xd = ix - ta[1];
	    yd = iy - ta[2];
	    gamma = xd * cosPA - yd * sinPA;
            absgammac = FMAX(MINR, pow(fabs(gamma), c));
	    delta = xd * sinPA + yd * cosPA;
            delota8 = delta/ta[8];
            absdeltac = FMAX(MINR, pow(fabs(delota8), c));
            fabsrc = FMAX(MINR, absgammac + absdeltac);

	    am = FMAX(MINR, pow(fabsrc, 1./c));

	    /******************************************************\
	    *                                                      *
            *  If we're near the center of the galaxy we'll        *
            *  subdivide the sampling.                             *
	    *                                                      *
	    \******************************************************/

            if (am < OVERSAMPR) {
                if (ta[4] * ta[8] <= 1 && am <= os) {

                    /* If scalelength is small, oversample by a lot more.  */

                    nsubsamp = IMIN(100, (int)(1./fabs(ta[4] * ta[8]) * 
			       20 / os));
                } else  if (am <= 3. * os)
                    nsubsamp = 20 / os;
                else {
                    /* If scalelength is large, bigger step */

                    nsubsamp = IMIN(100, NINT(20./am/os));
                };

	    } else
	        nsubsamp = 1;

            nsubsamp = IMAX (1, nsubsamp);

	    step = 1./nsubsamp;
	    xf = ix - 0.5 + 1./(2 * nsubsamp);
	    yf = iy - 0.5 + 1./(2 * nsubsamp);

	    sumdfdI0 = 0.;
	    xc = 0.;
	    while (xc <= nsubsamp - 1) {
	        yc = 0.;
	        xd = xf + xc * step - ta[1];
	        while (yc <= nsubsamp - 1) {
	            yd = yf + yc * step - ta[2];
		    gamma = xd * cosPA - yd * sinPA;
                    absgammac = FMAX (MINR, pow(fabs(gamma), c));
		    delta = xd * sinPA + yd * cosPA;
	            delota8 = delta/ta[8];
                    absdeltac = FMAX (MINR, pow(fabs(delota8), c));
	            fabsrc = FMAX(MINR, absgammac + absdeltac);
		    am = FMAX(MINR, pow(absgammac + absdeltac, 1./c));
	            gexp = -am / ta[4];
	            dfdI0 = exp(gexp);

		    flux = ta[3] * dfdI0;
		    dfdr = flux / -ta[4];
        	    theta = dfdr * am / (absgammac + absdeltac) / c;

		    if (ia[1] == 1) {
	                gammaprime = -c * pow(fabs(gamma), c-1) * cosPA;
	                deltaprime = -c * pow(fabs(delota8), c-1)/ta[8]*sinPA;
	                if (gamma < 0)
	                    gammaprime = -1 * gammaprime;  /* Derivative of  */
	                if (delta < 0)                     /* absolute value */
	                    deltaprime = -1 * deltaprime;

		        df->dpm[1][iy][ix] += theta*(gammaprime + deltaprime);
		    };

		    if (ia[2] == 1) {
	                gammaprime =  c * pow(fabs(gamma), c-1) * sinPA;
	                deltaprime =  -c * pow(fabs(delota8), c-1)/ta[8]*cosPA;
	                if (gamma < 0)
	                    gammaprime = -1 * gammaprime;  /* Derivative of  */
	                if (delta < 0)                     /* absolute value */
	                    deltaprime = -1 * deltaprime;

		        df->dpm[2][iy][ix] += theta*(gammaprime + deltaprime);
		    };

		    if (ia[3] == 1)
		        df->dpm[3][iy][ix] += dfdI0 * dI0dmag;   

		    if (ia[4] == 1)
		        df->dpm[4][iy][ix] += (flux*am/ta4sqr + dfdI0*dI0drs);

		    if (ia[8] == 1)
		        df->dpm[8][iy][ix] += (theta * -c * absdeltac /ta[8]
						         + dfdI0 * dI0dq);

		    if (ia[9] == 1) {
	                gammaprime = c * pow(fabs(gamma), c-1) * 
							-(xd*sinPA + yd*cosPA);
	                deltaprime = c * pow(fabs(delta), c-1)/ pow(ta[8],c) * 
                                                         (xd*cosPA - yd*sinPA);
	                if (gamma < 0)
	                    gammaprime = -1 * gammaprime;
	                if (delta < 0)
	                    deltaprime = -1 * deltaprime; 

		        df->dpm[9][iy][ix] -= PI/180. * theta * (gammaprime + 
								deltaprime);
		    };

	            if (ia[10] == 1) {
                        gammaprime = log(FMAX(MINR, fabs(gamma))) * absgammac;
                        deltaprime = log(FMAX(MINR, fabs(delota8)))*absdeltac;
                        df->dpm[10][iy][ix] += (dfdr *am* (-1./c/c * log(fabsrc)
                                 		+ 1./c/ fabsrc * (gammaprime + 
                                                 deltaprime)) + dfdI0*dI0dc);

	    	    };

		    sumdfdI0 += dfdI0;
		    yc = yc + 1.;
	        };
	        xc = xc+1.;
	    };

	    for (j = 1; j <= 10; j++) {
	        if (ia[j] == 1)
		    df->dpm[j][iy][ix] = df->dpm[j][iy][ix]/ nsubsamp/nsubsamp;
	    };

	    sumdfdI0 = sumdfdI0 / nsubsamp / nsubsamp;
	    flux = ta[3] * sumdfdI0;
	    df->dpm[0][iy][ix] = flux;               /* ... for convlution   */
        };
    };

/* #if DEBG
    sprintf (model->name, "model-0.fits");    Output image
    writefits ("test.fits", model, 0);
#endif  */

}

