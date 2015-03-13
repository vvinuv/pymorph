#include <math.h>
#include <curses.h>
#include <stdio.h>
#include "nrutil.h"
#include "structs.h"

#define NCUTOFF 10.

void smooth (float **in, float **out, long naxes[], float sigma, float radius);
void mkimg (float **array, long naxes[], char *outname, char *tname);

struct image *fakenoise (struct image *d)
{
    int i, j, npix;
    float var, varp, varn, effgain, mean=0., pedestal=0.;
    struct image *noise;
    extern char device[];

    noise = (struct image *) calloc (1, sizeof (struct image));
    if (noise == (struct image *) 0) {
        pfunc ("Error in allocating memory for mini-image....");
        if (strncmp(device, "curses", 6)== 0)
            refresh();
        exit (1);
    };

    npix = d->naxes[1] * d->naxes[2]; 
    noise->naxes[1] = d->naxes[1];
    noise->naxes[2] = d->naxes[2];

    noise->z = matrix (1, noise->naxes[2], 1, noise->naxes[1]);
    effgain = FMAX(d->gain * d->ncombine, 1.);

    smooth (d->z, noise->z, noise->naxes, 2., 3.);

    /******************************************************************\ 
    *  Make sure the sky is not below 0.  Make sure it's at least 10.  *
    *  Mean will do -- doesn't have to be accurate.                    *
    \******************************************************************/

    for (j=1; j <= noise->naxes[2]; ++j) {
        for (i=1; i <= noise->naxes[1]; ++i)
	    mean += noise->z[j][i];
    };
    mean = mean / noise->naxes[1] / noise->naxes[2];
    if (mean < 10.) pedestal = 10. - mean;

    /*********************************\
    *  Now create the variance array  *
    \*********************************/

    for (j=1; j <= noise->naxes[2]; ++j) {
        for (i=1; i <= noise->naxes[1]; ++i) {
            varp = (noise->z[j][i] + pedestal) * effgain;

	    if (varp <= -3 * d->rdnoise) 
	        varp = 1.e30;

            varn = d->rdnoise * d->rdnoise * d->ncombine;
	    var = varp + varn;
	    if (var <= 1.)
	        var = 1.;

 	    noise->z[j][i] = sqrt(var) / effgain; 
        };
    };

    return (noise);

}


/*
    sprintf (noise->name, "noise.fits");
    mkimg (noise->z, noise->naxes, noise->name, "noise.fits");
*/
