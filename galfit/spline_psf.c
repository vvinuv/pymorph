#include "nrutil.h"
#include "structs.h"

void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a);

float **spline_psf (struct image *psf, float *x1a, float *x2a)
{
    int i;
    float **y2a;

    for (i=1; i <= psf->naxes[1]; i++) 
	x2a[i] = i;
    for (i=1; i <= psf->naxes[2]; i++) 
	x1a[i] = i;

    y2a = matrix (1, psf->naxes[2], 1, psf->naxes[1]);
    splie2 (x1a, x2a, psf->z, psf->naxes[2], psf->naxes[1], y2a);
    return (y2a);
}
