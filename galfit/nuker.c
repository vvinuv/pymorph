#include <math.h>
#include <stdlib.h>
#include "structs.h"
#include "nrutil.h"
#include "const.h"
#include "debug.h"

#define OVERSAMPR 5.
#define peaksamp 50
#define dr 0.05 / peaksamp
#define MINR 1.e-15

void writefits(char *phead, struct image *img, char *object, int add);
float rdist (float x, float y, float a[]);
float egrid (float (*nuker)(float a[], float x, float y), float xmin, 
                        float xmax, float ymin, float ymax, float a[]);

void nuker (float a[], int ia[], struct image *model, struct derivs *df)
{
    extern struct inpars input;

    unsigned long i, ix, iy;
    float nukefunc (float ta[], float x, float y);
    float xd, yd, gamma, delta, am, xdf, ydf, af, ta[NPARS+1], cosPA,
          sinPA, xf, yf, xc, yc, gammaf, deltaf, magfac, flux, B, C, D,
	  E, F, a3zpt, ta4r, rta4, rta4ta5, ta4rta7, powDE, 
	  dIdr, ta8sqr, r, delota8, cee, theta, gammaprime, deltaprime,
          absgammac, absdeltac, fabsrc, absgammafc, absdeltafc, nsubsamp, 
          step;
    int j, os;

/*  Ib=ta[3], rb=ta[4], alpha=ta[5], beta=ta[6], gamma=ta[7]  */

    os = input.sampfac;

    a[5] = FMAX (a[5], 0.01);
    for (i=1; i <= NPARS; i++)
        ta[i] = a[i];

    cee = ta[10] + 2.;

    /*****************************************************************\
    *                                                                 *
    *  Pre-calculate some things we will frequently use later on but  *
    *  which needs to be calculated one time only.                    *
    *                                                                 *
    \*****************************************************************/

    /* a[3] is the galaxy surface brightness at break radius   */
    /* while ta[3] is the number count at that radius.         */

    B = pow(2, (ta[6] - ta[7]) /ta[5]);
    E = ( ta[7] - ta[6] ) / ta[5];
    a3zpt = (ta[3] - model->muzpt)/-2.5;
    ta8sqr = ta[8] * ta[8];
    
    ta[3] = pow(10., a3zpt);
    ta[9] = -1. * ta[9] / 180. * PI + PI/2.;    /* "Intuitive" orientation */
    cosPA = cos (ta[9]);
    sinPA = sin (ta[9]);

    magfac = log(10.) / -2.5;

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
            absgammac = FMAX(MINR, pow (fabs(gamma), cee));
	    delta = xd * sinPA + yd * cosPA;
            delota8 = delta/ta[8];
            absdeltac = FMAX(MINR, pow(fabs(delota8), cee));
            fabsrc = FMAX(MINR, absgammac + absdeltac);

	    am = FMAX(MINR, pow(fabsrc, 1./cee));

	    /******************************************************\
	    *                                                      *
            *  If we're near the center of the galaxy we'll        *
            *  subdivide the sampling.                             *
	    *                                                      *
	    \******************************************************/

	    if (am < OVERSAMPR * os)
                flux = egrid (nukefunc, ix-0.5, ix+0.5, iy-0.5, iy+0.5, ta);

	    /****************************************\
	    *                                        *
	    *  Don't need to do fine sampling here   *
	    *                                        *
	    \****************************************/

            ta4r = ta[4]/am;
            rta4 = am/ta[4];
	    ta4rta7 = pow(ta4r, ta[7]);
	    C = ta[3] * ta4rta7;
	    rta4ta5 = pow(rta4, ta[5]);
	    D = 1. + rta4ta5;
	    powDE = pow(D, E);
	    F = C * powDE;

	    if (am >= OVERSAMPR * os)
	        flux = B * F;

	    df->dpm[0][iy][ix] = flux;               /*  ... for convlution  */

            dIdr = flux / am * (-ta[7] + E / D * rta4ta5 * ta[5]);
            theta = dIdr * am / (absgammac + absdeltac) / cee;

            if (ia[1] == 1) {
                gammaprime = -cee * pow(fabs(gamma), cee-1) * cosPA;
                deltaprime = -cee * pow(fabs(delota8), cee-1) / ta[8] * sinPA;
                if (gamma < 0)
                    gammaprime = -1 * gammaprime;   /* Derivative of   */
                if (delta < 0)                      /* absolute value  */
                    deltaprime = -1 * deltaprime;

                df->dpm[1][iy][ix] = theta*(gammaprime + deltaprime);
            };

            if (ia[2] == 1) {
                gammaprime =  cee * pow(fabs(gamma), cee-1) * sinPA;
                deltaprime =  -cee * pow(fabs(delota8), cee-1) / ta[8] * cosPA;
                if (gamma < 0)
                    gammaprime = -1 * gammaprime;   /* Derivative of   */
                if (delta < 0)                      /* absolute value  */
                    deltaprime = -1 * deltaprime;

                df->dpm[2][iy][ix] = theta*(gammaprime + deltaprime);
            };

	    if (ia[3] == 1)
	        df->dpm[3][iy][ix] = B * ta4rta7 * powDE * magfac * ta[3];
	    if (ia[4] == 1)
	        df->dpm[4][iy][ix] = flux * (ta[7] / ta[4] - E / D * ta[5] *
           	  			rta4ta5 / ta[4]);
	    if (ia[5] == 1)
	        df->dpm[5][iy][ix] = flux * E/ta[5] * (log (2.) - log (D) +
				     E * ta[5]/D * rta4ta5 * log (rta4) );
	    if (ia[6] == 1)
	        df->dpm[6][iy][ix] = flux * (log (2.)/ta[5] + log(D) * -1./ta[5]);
	    if (ia[7] == 1)
	        df->dpm[7][iy][ix] = flux * (-log (2.)/ta[5] + log(ta4r) + 
							log(D)/ta[5]);
            if (ia[8] == 1)
                df->dpm[8][iy][ix] = theta * (-cee * absdeltac /ta[8]);

            if (ia[9] == 1) {
                gammaprime = cee * pow(fabs(gamma), cee-1) * 
						-(xd*sinPA + yd*cosPA);
                deltaprime = cee * pow(fabs(delta), cee-1) / pow(ta[8],cee) * 
						(xd*cosPA - yd*sinPA);
                if (gamma < 0)
                    gammaprime = -1 * gammaprime;
                if (delta < 0)
                    deltaprime = -1 * deltaprime; 

                df->dpm[9][iy][ix] = -PI/180.*theta* (gammaprime + deltaprime);
            };

            if (ia[10] == 1) {
                gammaprime = log(FMAX (MINR, fabs(gamma))) * absgammac;
                deltaprime = log(FMAX (MINR, fabs(delota8))) * absdeltac;
                df->dpm[10][iy][ix] = dIdr * am * (-1./cee/cee * log(fabsrc)
                       		 + 1./cee/ fabsrc * (gammaprime + deltaprime));
            };
        };
    };

#if DEBUG
    sprintf (model->name, "model-0.fits");  /* Output image */
    writefits ("test.fits", model, "", 0);
#endif

}

float nukefunc (float ta[], float x, float y)
{
    float nuke, r, c, xd, yd, Ib;

    r = rdist (x, y, ta);
    nuke = pow (2, ( ta[6] - ta[7])/ta[5]) * ta[3] * pow((ta[4] / r), ta[7]) *
                pow(1 + pow((r/ta[4]), ta[5]), (ta[7]-ta[6])/ta[5]);
    return (nuke); 
}
