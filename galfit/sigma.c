#include <math.h>
#include <curses.h>
#include <stdio.h>
#include "nrutil.h"
#include "structs.h"

#define NCUTOFF 10.

void smooth (float **in, float **out, long naxes[], float sigma, float radius);
void mkimg (float **array, long naxes[], char *outname, char *tname);

struct image *sigma_img (struct image *d)
{
    int i, j, npix;
    float var, varp, varn, effgain, mean=0., pedestal=0.;
    struct image *sigma;
    extern char device[];

    sigma = (struct image *) calloc (1, sizeof (struct image));
    if (sigma == (struct image *) 0) {
        pfunc ("Error in allocating memory for mini-image....");
        if (strncmp(device, "curses", 6)== 0)
            refresh();
        exit (1);
    };

    npix = d->naxes[1] * d->naxes[2]; 
    sigma->naxes[1] = d->naxes[1];
    sigma->naxes[2] = d->naxes[2];

    sigma->z = matrix (1, sigma->naxes[2], 1, sigma->naxes[1]);
    effgain = FMAX(d->gain * d->ncombine, 1.);

    smooth (d->z, sigma->z, sigma->naxes, 2., 3.);

    /******************************************************************\ 
    *  Make sure the sky is not below 0.  Make sure it's at least 10.  *
    *  Mean will do -- doesn't have to be accurate.                    *
    \******************************************************************/

    for (j=1; j <= sigma->naxes[2]; ++j) {
        for (i=1; i <= sigma->naxes[1]; ++i)
	    mean += sigma->z[j][i];
    };
    mean = mean / sigma->naxes[1] / sigma->naxes[2];
    if (mean < 10.) pedestal = 10. - mean;

    /*********************************\
    *  Now create the variance array  *
    \*********************************/

    for (j=1; j <= sigma->naxes[2]; ++j) {
        for (i=1; i <= sigma->naxes[1]; ++i) {
            varp = (sigma->z[j][i] + pedestal) * effgain;

	    if (varp <= -3 * d->rdnoise) 
	        varp = 1.e30;

            varn = d->rdnoise * d->rdnoise * d->ncombine;
	    var = varp + varn;
	    if (var <= 1.)
	        var = 1.;

 	    sigma->z[j][i] = sqrt(var) / effgain; 
        };
    };

    return (sigma);

}


/*
    sprintf (sigma->name, "sigma.fits");
    mkimg (sigma->z, sigma->naxes, sigma->name, "sigma.fits");
*/
