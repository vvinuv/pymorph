#include <math.h>
#include <stdlib.h>
#include "structs.h"
#include "mymath.h"
#include "nrutil.h"
#include "const.h"

#define MINSZ 3.0
#define MAXRAD 1000.
#define OVERSAMPR 10.
#define MINR 1.e-15

float ratio (float c);

void king (float a[], int ia[], struct image *model, struct derivs *df)
{
    extern struct inpars input;

    unsigned long i;
    int j, ix, iy, nsubsamp, os;

    float xd, yd, am, ta[NPARS+1], step, delota8, xf, yf, xc, yc, ra, rb,
	  rdiff, dfdA, flux, rsqr, rcsqr, rtsqr, dIdr, cosPA, sinPA, gamma,
	  delta, c, gammaprime, deltaprime, absgammac, absdeltac, fabsrc,
	  theta, err, coeff, sumdfdA, zpt, magfac;

    float dfridr (float (*ratio)(float n), float x, float h, float *err);

    os = input.sampfac;
    
    for (i=1; i <= NPARS; i++)
        ta[i] = a[i];

    c = ta[10] + 2.;

    ta[4] = fabs(ta[4]);
    rcsqr = ta[4] * ta[4];
    rtsqr = ta[5] * ta[5];

    ta[9] = -1. * ta[9] / 180. * PI + PI/2.;   /* "Intuitive" orientation */
    cosPA = cos (ta[9]);
    sinPA = sin (ta[9]);

    zpt = (ta[3] - model->muzpt)/-2.5;
    ta[3] = pow (10., zpt);
    magfac = log(10.)/-2.5;

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
	    *  Create the model psf.  If we're near the center     *
            *  of the PSF we'll subdivide the sampling.            *
	    *                                                      *
	    \******************************************************/
	    
            if (am < OVERSAMPR) {            
                if (ta[4] * ta[8] <= 1. && am <= os) {
	    
                    /* If scalelength is small, oversample by a lot more.  */
	    
                    nsubsamp = IMIN(100, (int)(1./fabs(ta[4] * ta[8]) * 
			       40 / os));
                } else  if (am <= 3. * os)
                    nsubsamp = 20/os;
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
	    
	    sumdfdA = 0.;
	    xc = 0.;

	    /*  Truncate the fit at the tidal radius  */	    

	    while (xc <= nsubsamp - 1 && am < ta[5]) {
	        yc = 0.;
	        xd = xf + xc*step - ta[1];
	        while (yc <= nsubsamp - 1) {
		    yd = yf + yc*step - ta[2];
                    gamma = xd * cosPA - yd * sinPA;
                    absgammac = FMAX(MINR, pow (fabs(gamma), c));
                    delta = xd * sinPA + yd * cosPA;
                    delota8 = delta/ta[8];
                    absdeltac = FMAX(MINR, pow (fabs(delota8), c));
	            fabsrc = FMAX(MINR, absgammac + absdeltac);
                    am = FMAX(MINR, pow(absgammac + absdeltac, 1./c));
	            rsqr = am * am;
	        
		    ra = 1+rsqr/rcsqr;
		    rb = 1+rtsqr/rcsqr;
		    rdiff =  1./pow (ra, 1./ta[6]) - 1./pow(rb, 1./ta[6]);
		    dfdA = pow (rdiff, ta[6]);
	            flux = ta[3] * dfdA;

		    coeff = ta[6] * ta[3] * pow (rdiff, ta[6] - 1);
	            dIdr = coeff * -2 * am / (ta[6] * rcsqr * 
			   pow (ra, 1/ta[6] + 1));
                    theta = dIdr * am / (absgammac + absdeltac) / c;
	        
	    	    if (ia[1] == 1) {
                        gammaprime = -c * pow(fabs(gamma), c-1) * cosPA;
                        deltaprime = -c * pow(fabs(delota8), c-1) /ta[8]*sinPA;
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
	                df->dpm[3][iy][ix] += flux * magfac;
		    
	            if (ia[4] == 1) 
	                df->dpm[4][iy][ix] += coeff * 2/ta[6] / rcsqr / ta[4]* 
					      (rsqr / pow(ra, 1. / ta[6] + 1) -
					       rtsqr/ pow(rb, 1. / ta[6] + 1));
		    
	            if (ia[5] == 1) 
	                df->dpm[5][iy][ix] += coeff * 2 * ta[5] / rcsqr / 
					      pow(rb, 1./ta[6] + 1);
		    
	            if (ia[6] == 1) 
	                df->dpm[6][iy][ix] += flux*(log(rdiff)+1/ta[6]/rdiff* 
					      (log (ra) / pow(ra, 1./ta[6]) -
					       log (rb) / pow(rb, 1./ta[6])));
		    
                    if (ia[8] == 1)
                        df->dpm[8][iy][ix] += (theta * -c * absdeltac /ta[8] );
		    
                    if (ia[9] == 1) {
                        gammaprime = c * pow(fabs(gamma), c-1) * 
		    					-(xd*sinPA + yd*cosPA);
                        deltaprime = c * pow(fabs(delta), c-1) / 
                                          pow(ta[8],c) * (xd*cosPA - yd*sinPA);
                        if (gamma < 0)
                            gammaprime = -1 * gammaprime;
                        if (delta < 0)
                            deltaprime = -1 * deltaprime; 

                        df->dpm[9][iy][ix] -= PI/180. * theta * (gammaprime + 
								deltaprime);
                    };

                    if (ia[10] == 1) {
                        gammaprime = log(FMAX(MINR, fabs(gamma))) * absgammac;
                        deltaprime = log(FMAX(MINR, fabs(delota8))) * 
								   absdeltac;
                        df->dpm[10][iy][ix] += (dIdr *am* (-1./c/c* log(fabsrc)
                                      		+ 1./c/ fabsrc * (gammaprime + 
                                                    deltaprime)));
                    };
		    
	            sumdfdA += dfdA;
		    yc = yc + 1.;
	        };	    
	        xc = xc + 1.;
	    };	    
		    
            for (j = 1; j<= 10; j++) {
                if (ia[j] == 1)
                        df->dpm[j][iy][ix] = df->dpm[j][iy][ix] /nsubsamp/
								    nsubsamp;
            };	    
	    
	    sumdfdA = sumdfdA / (nsubsamp * nsubsamp);
	    flux = sumdfdA * ta[3];
	    df->dpm[0][iy][ix] = flux;               /*  ... for convlution  */
        };
    };		    
}	    
	    
	    
	    
	    
	    
	    
	    
