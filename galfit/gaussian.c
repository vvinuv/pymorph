#include <math.h>
#include <stdlib.h>
#include "structs.h"
#include "mymath.h"
#include "nrutil.h"
#include "const.h"

#define MINSZ 10.0
#define OVERSAMPR 10.0
#define MINR 1.e-15

float ratio (float c);
float dfridr (float (*ratio)(float n), float x, float h, float *err);

void gaussian (float a[], int ia[], struct image *model, struct derivs *df)
{
    extern struct inpars input;

    unsigned long i;
    int ix, iy, j, nsubsamp, os;
    float ta[NPARS+1]; 
    float xc, yc, flux, rsqr, ta4sqr, dIdmag, step, 
	  dfdr, xf, yf, dfdI, gexp, totcnts, dFWdsig, gamma, delta, am,
          xd, yd, cosPA, sinPA, delota8, ta8sqr, fabsrc, absgammac, 
          absdeltac, theta, area_ratio, c, gammaprime, deltaprime, tempsqr,
          totcounts, dfdsig, err, dc=0.1, dratio_dc, sumdfdI; 
    float dIdq, dIdc, dIdsig;

    os = input.sampfac;

    for (i=1; i <= NPARS; i++)
        ta[i] = a[i];

    c = ta[10] + 2.;

    ta[4] = fabs(ta[4]);

    dFWdsig = 2.3548200;             /*  2 * sqrt (2. * log(2.));  */
    ta[4] = ta[4] / dFWdsig;  			/* FWHM -> sigma */
    ta4sqr = ta[4] * ta[4];
    ta8sqr = ta[8] * ta[8];

    /**********************************************************\
    *  Notice the last term ta[8] in the third equation.       *
    *  That's because the distance from the center is defined  *
    *  as  r = (|x|^c + |y/q|^c)^(1/c), where q is the axis    *
    *  ratio.  This would make the magnitude come out right.   *
    \**********************************************************/

    area_ratio = ratio(c);
    dratio_dc = dfridr (ratio, c, dc, &err);

    totcnts = pow (10., (model->magzpt - a[3])/2.5);
    ta[3] =  totcnts / 2. / PI / (ta4sqr * ta[8] / area_ratio);
    dIdmag = ta[3] * log(10.) / -2.5;

    dIdsig = -2 * ta[3] / ta[4];
    dIdq = -ta[3] / ta[8]; 
    dIdc = ta[3] / area_ratio * dratio_dc;

    ta[9] = -1. * ta[9] / 180. * PI + PI/2.;    /* "Intuitive" orientation */
    cosPA = cos (ta[9]);
    sinPA = sin (ta[9]);

    for (iy=1; iy <= df->naxes[2]; iy++) {
        for (ix=1; ix <= df->naxes[1]; ix++) {

	    /* calculate the elliptical center distance */
	    xd = ix - ta[1];
	    yd = iy - ta[2];
            gamma = xd * cosPA - yd * sinPA;
            absgammac = FMAX(MINR, pow(fabs(gamma), c));
            delta = xd * sinPA + yd * cosPA;
            delota8 = delta / ta[8];
            absdeltac = FMAX(MINR, pow(fabs(delota8), c));
            fabsrc = FMAX(MINR, absgammac + absdeltac);
	    
            am = FMAX(MINR, pow(fabsrc, 1./c));
	    
	    rsqr = am * am;
	    gexp = -1. * rsqr / ta4sqr / 2.;
	    dfdI = exp(gexp);
	    
	    /***************************************************\
	    *                                                   *
	    *  Create the model.  If we're near the center      *
            *  of the model we'll subdivide the sampling.       *
	    *                                                   *
	    \***************************************************/
	    
            if (am < OVERSAMPR) {
                if (ta[4] * ta[8] <= 1. && am <= os) {
	        
                    /* If scalelength is small, oversample by a lot more.  */
	        
                    nsubsamp = IMIN(100,(int)(1./fabs(ta[4]*ta[8])*2*
			       OVERSAMPR/os));
                } else  if (am <= 3. * os)
                    nsubsamp = 2 * OVERSAMPR / os;
                else {
                    /* If scalelength is large, bigger step */
		    
                    nsubsamp = IMIN(100, NINT(20./am/os));
                };	    
	    } else
	        nsubsamp = 1;

            nsubsamp = IMAX (1, nsubsamp);

            xf = ix - 0.5 + 1./(2 * nsubsamp);
            yf = iy - 0.5 + 1./(2 * nsubsamp);
            step = 1./nsubsamp;
	    
	    sumdfdI = 0.;
	    xc = 0.;
	    while (xc <= nsubsamp - 1) {
	        yc = 0.;
	        xd = xf + xc * step - ta[1];
	        while (yc <= nsubsamp - 1) {
	      	    yd = yf + yc * step - ta[2];
                    gamma = xd * cosPA - yd * sinPA;
                    absgammac = FMAX(MINR, pow(fabs(gamma), c));
                    delta = xd * sinPA + yd * cosPA;
                    delota8 = delta / ta[8];
                    absdeltac = FMAX(MINR, pow(fabs(delota8), c));
	            fabsrc = FMAX(MINR, absgammac + absdeltac);
                    am = FMAX(MINR, pow(absgammac + absdeltac, 1./c));
		    rsqr = am * am;
                    gexp = -1. * rsqr / ta4sqr / 2.;

                    dfdI = exp(gexp);
		    flux = ta[3] * dfdI;
                    dfdr = -ta[3] * dfdI * am / ta4sqr;   
                    theta = dfdr * am / (absgammac + absdeltac)/ c;  

                    if (ia[1] == 1) {
                        gammaprime = -c * pow(fabs(gamma), c-1) * cosPA;
                        deltaprime = -c * pow(fabs(delota8), c-1)/ta[8]*sinPA;
                        if (gamma < 0)
                            gammaprime = -1 * gammaprime;  /* Deriv. of      */
                        if (delta < 0)                     /* absolute value */
                            deltaprime = -1 * deltaprime;
	    
                        df->dpm[1][iy][ix] += theta*(gammaprime + deltaprime);
                    };
	        
                    if (ia[2] == 1) {
                        gammaprime =  c * pow(fabs(gamma), c-1) * sinPA;
                        deltaprime =  -c * pow(fabs(delota8), c-1) / ta[8] * cosPA;
                        if (gamma < 0)
                            gammaprime = -1 * gammaprime;  /* Deriv. of      */
                        if (delta < 0)                     /* absolute value */
                            deltaprime = -1 * deltaprime;
		        
                        df->dpm[2][iy][ix] += theta*(gammaprime + deltaprime);
                    };
          
                    dfdsig = -flux * 2 * gexp/ta[4];
                    if (ia[3] == 1)
	                df->dpm[3][iy][ix] += dfdI * dIdmag;
	      
	            if (ia[4] == 1)
                        df->dpm[4][iy][ix] += (dfdsig + dfdI * dIdsig)/dFWdsig;

    	            if (ia[8] == 1) 
                  	    df->dpm[8][iy][ix] += (theta * -c * absdeltac /
	          					ta[8] + dfdI * dIdq);
	          
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
                        deltaprime = log(FMAX(MINR, fabs(delota8)))* absdeltac;
                        df->dpm[10][iy][ix] += (dfdr *am* (-1./c/c* log(fabsrc)
                                 	     + 1./c/ fabsrc * (gammaprime + 
                                                deltaprime)) + dfdI*dIdc);
                    };

		    sumdfdI += dfdI;
	            yc = yc + 1.;
	        };
	        xc = xc + 1.;
	    };

            for (j=1; j<= NPARS; j++) {
                if (ia[j] ==1)
                    df->dpm[j][iy][ix] = df->dpm[j][iy][ix]/ nsubsamp/nsubsamp;
	    						
	    };

	    sumdfdI = sumdfdI / nsubsamp / nsubsamp;
            flux = ta[3] * sumdfdI;
	    df->dpm[0][iy][ix] = flux;               /*  ... for convlution  */
	};
    };
}

