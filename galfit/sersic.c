#include <math.h>
#include <stdlib.h>
#include "structs.h"
#include "mymath.h"
#include "nrutil.h"
#include "const.h"
#include "debug.h"

#define OVERSAMPR 16.
#define MINR 1.e-15

void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
float gammln(float xx);
void writefits(char *phead, struct image *img, char *object, int add);
float ratio (float c);

void sersic (float a[], int ia[], struct image *model, struct derivs *df)
{
    float gamm2nln(float n);
    float getinterp(float n);

    extern float Y2[];    
    extern float MU[];
    extern float TWOn[];
    extern ninterp;
    extern struct inpars input;

    float sconst=3.60729;
    float mu=7.66925;
    unsigned long i, ix, iy;
    int j, nsubsamp, os;
    float xd, yd, gamma, delta, am, xdf, ydf, af, ta[NPARS+1], xf, yf, cosPA,
	  sinPA, xc, yc, rr, h, err, invp5, invp5m2, p4invp5, fabsrc, magfac,
	  muinvp5, twodp4, mufac, theta, dgammlndn, dmudp5, flux, c,
	  delota8, rc, area_ratio, absgammac, absdeltac, gammaprime,
	  deltaprime, dfdr, step, dIedmag, dratio_dc, dc = 0.1, dfdIe, 
	  totcnt, dredmag, sumdfdIe, drdc;
    float dIedn, dIedq, dIedc, dIedre;

    float dfridr (float (*somefunc)(float n), float x, float h, float *err);

    os = input.sampfac;

    for (i=1; i <= NPARS; i++)
        ta[i] = a[i];

    c = ta[10] + 2.;

    if (ta[5] >= 5.) 
        mu = 1.9992 * ta[5] - 0.3271;
    else
        splint (TWOn, MU, Y2, ninterp, 2*ta[5], &mu);
    sconst = ta[5] * exp(gammln(2.*ta[5])) * exp(mu) * pow(mu, -2.*ta[5]);

    /*****************************************************************\
    *                                                                 *
    *  Pre-calculate some things we will frequently use later on but  *
    *  which needs to be calculated one time only.                    *
    *                                                                 *
    \*****************************************************************/

    /* a[3] is the galaxy total brightness while ta[3] */
    /* is the count rate at the effective radius.      */

    /* Notice the last term ta[8] in the first equation here.      */
    /* That's because the distance from the center is defined      */
    /* as  r = (|x|^c + |y/q|^c)^(1/c), where q is the axis ratio. */
    /* This would make the magnitude come out right.               */

    /* The area_ratio is the ratio of a pure ellipse area          */
    /* to that of a disky or boxy area.  Since we can calculate    */
    /* the flux of an elliptical model, the flux of a disky/boxy   */
    /* model is simply the ratio of the two areas.                 */

    area_ratio = ratio (c);
    dratio_dc = dfridr (ratio, c, dc, &err);

    totcnt = pow (10., (model->magzpt - a[3])/ 2.5) ;
    ta[3] = totcnt / 2 / PI / (sconst * ta[4] * ta[4] * ta[8] / area_ratio); 

/*  ta[9] = -1. * a[9] / 180. * PI;    This is for Sextractor    */

    ta[9] = -ta[9] / 180. * PI + PI/2.;    /* "Intuitive" orientation */
    cosPA = cos (ta[9]);
    sinPA = sin (ta[9]);

    if (ia[5] == 1) {
	if (ta[5] > 5.) 
	    dmudp5 = 1.9992;
	else {
	    h = 0.1;	/* Hardwire inital step size for derivative */
	    dmudp5 = dfridr (getinterp, ta[5], h, &err);
	};
	if (ta[5] > 3.)
	    h = log(ta[5]) * 0.1;  /* Hardwire initial step size for deriv. */
	else
	    h = 0.1;
	dgammlndn = dfridr (gamm2nln, ta[5], h, &err);
    };

    invp5 = 1./ta[5];
    invp5m2 = invp5 - 2.;
    p4invp5 = pow(ta[4], invp5);
    magfac = log(10.) / -2.5;
    muinvp5 = mu * invp5 / pow(ta[4], invp5+1.);
    twodp4 = 2./ta[4];
    mufac = mu * invp5 * invp5;

    /**********************************************************\
    *  Calculate some of the cross coupling terms between Ie   *
    *  with all the other parameters.                          *
    \**********************************************************/

    dIedmag = ta[3] * log(10.) / -2.5;

    dIedre = -2 * ta[3] / ta[4];
    dIedq = -ta[3] / ta[8]; 
    dIedc = ta[3] / area_ratio * dratio_dc;
    dIedn = -ta[3] * (invp5 + dgammlndn + dmudp5 - 
				 2*(log(mu) + ta[5] / mu * dmudp5) );

/* printf ("dIedmag=%f   dIedq=%f   dIedre=%f   dIedc=%f\n", dIedmag, dIedq, dIedre, dIedc); */

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
                if (ta[4] * ta[8] <= 1. && am <= os) {

		    /* If scalelength is small, oversample by a lot more.  */

                    nsubsamp = IMIN(100, (int)(1./fabs(ta[4]*ta[8])*2*
			       OVERSAMPR/os));
                } else  if (am <= 4. * os)
                    nsubsamp = 2 * OVERSAMPR / os;
	        else {
	            /* If scalelength is large, bigger step */

                    nsubsamp = IMIN(100, NINT(2*OVERSAMPR/am/os));
	        };
	    } else
	        nsubsamp = 1;

	    nsubsamp = IMAX (1, nsubsamp);

	    xf = ix - 0.5 + 1./(2 * nsubsamp);
	    yf = iy - 0.5 + 1./(2 * nsubsamp);
            step = 1./nsubsamp;

            sumdfdIe = 0.;
	    xc = 0.;
 	    while (xc <= nsubsamp - 1) {
	        yc = 0.;
	        xd = xf + xc * step - ta[1];
	        while (yc <= nsubsamp - 1) {
	            yd = yf + yc * step - ta[2];
	            gamma = xd * cosPA - yd * sinPA;
                    absgammac = FMAX(MINR, pow(fabs(gamma), c));
	            delta = xd * sinPA + yd * cosPA;
	            delota8 = delta/ta[8];
                    absdeltac = FMAX(MINR, pow(fabs(delota8), c));
     	            fabsrc = FMAX(MINR, absgammac + absdeltac);
	            am = FMAX(MINR, pow (absgammac + absdeltac, 1./c));
	            rr = pow(am/ta[4], 1./ta[5]);
	            dfdIe = exp(-mu*(rr-1.));
                    flux = ta[3] * dfdIe;
                    dfdr = flux * -mu / ta[5] / ta[4] * pow(am/ta[4], 
								1/ta[5] - 1);
   	            theta =  dfdr * am / (absgammac + absdeltac) / c;

		    if (ia[1] == 1) {
	                gammaprime = -c * pow(fabs(gamma), c-1) * cosPA;
	                deltaprime = -c * pow(fabs(delota8), c-1) / ta[8] * 
									sinPA;
	                if (gamma < 0)
	                    gammaprime = -1 * gammaprime;  /* Derivative of  */
	                if (delta < 0)                     /* absolute value */
	                    deltaprime = -1 * deltaprime;
		        df->dpm[1][iy][ix] += theta*(gammaprime + deltaprime);
	            };

		    if (ia[2] == 1) {
	                gammaprime =  c * pow(fabs(gamma), c-1) * sinPA;
	                deltaprime =  -c * pow(fabs(delota8), c-1) / ta[8] * cosPA;
	                if (gamma < 0)
	                    gammaprime = -1 * gammaprime;  /* Derivative of  */
	                if (delta < 0)                     /* absolute value */
	                    deltaprime = -1 * deltaprime;
		        df->dpm[2][iy][ix] += theta*(gammaprime + deltaprime);
	            };

		    if (ia[3] == 1) {
/*		        dredmag = 1./(-2.5/log(10.)*mu / ta[5]*
					pow(am/ta[4],1./ta[5])/ta[4]); */
		        df->dpm[3][iy][ix] += dfdIe * dIedmag;
		    };

		    if (ia[4] == 1)
		        df->dpm[4][iy][ix] += (flux * muinvp5 * pow(am, invp5)
						     	    + dfdIe * dIedre);
		    if (ia[5] == 1)
		        df->dpm[5][iy][ix] += (flux * (mufac * pow(am/ta[4], 
				invp5) * log(am/ta[4]) - dmudp5 * 
				(pow(am/ta[4], invp5) - 1.)) + dfdIe * dIedn);

		    if (ia[8] == 1)
		        df->dpm[8][iy][ix] += (theta * -c * absdeltac/ta[8] + 
							   + dfdIe * dIedq);

		    if (ia[9] == 1) {
	                gammaprime = c * pow(fabs(gamma), c-1) * 
						-(xd*sinPA + yd*cosPA);
	                deltaprime = c * pow(fabs(delta), c-1)/pow(ta[8],c) * 
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
	                deltaprime = log(FMAX(MINR, fabs(delota8))) *absdeltac;
		        drdc = am* (-1./c/c * log(fabsrc) + 1./c/ fabsrc * 
					(gammaprime + deltaprime));
        	        df->dpm[10][iy][ix] += (dfdr * drdc + dfdIe*dIedc);
		    };

		    sumdfdIe += dfdIe;
	            yc = yc + 1.;
	        };
	        xc = xc + 1.;
	    };

	    for (j = 1; j<= 10; j++) {
	        if (ia[j] == 1)
		    df->dpm[j][iy][ix] = df->dpm[j][iy][ix] /nsubsamp/nsubsamp;
	    };

            sumdfdIe = sumdfdIe / nsubsamp / nsubsamp;
	    flux = ta[3] * sumdfdIe;
	    df->dpm[0][iy][ix] = flux;               /* ... for convlution   */
        };
    };

/* #if DEBUG
    sprintf (model->name, "model-0.fits");   Output image 
    writefits ("test.fits", model, "", 0);
#endif  */

}


float gamm2nln(float n)
{

    float n2, g2nln;

    n2 = 2. * n;
    g2nln = gammln(n2);
    return (g2nln);
}


float getinterp(float n)
{
    extern float TWOn[], MU[], Y2[];
    extern int ninterp;
    float n2, y;

    n2 = 2. * n;
    splint (TWOn, MU, Y2, ninterp, n2, &y);
    return (y);

}
