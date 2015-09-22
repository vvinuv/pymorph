#include "structs.h"
#include "nrutil.h"

void psf_derivs (struct image *psf, float **y1a, float **y2a, float **y12a)
{
    int j, k;

    y1a = matrix (1, psf->naxes[2], 1, psf->naxes[1]);
    y2a = matrix (1, psf->naxes[2], 1, psf->naxes[1]);
    y12a = matrix (1, psf->naxes[2], 1, psf->naxes[1]);
    
    for (j = 1; j <= psf->naxes[2]; j++) {
	for (k = 1; k <= psf->naxes[1]; k++) {
	    if (j==1 || k==1 || j == psf->naxes[2] || k == psf->naxes[1]) {
		y1a[j][k] = 0.;
		y2a[j][k] = 0.;
		y12a[j][k] = 0.;
	    } else {
	        y1a[j][k] = (psf->z[j][k+1] - psf->z[j][k-1])/2.;
	        y2a[j][k] = (psf->z[j+1][k] - psf->z[j-1][k])/2.;
	        y12a[j][k] = (psf->z[j+1][k+1] - psf->z[j+1][k-1] - 
			      psf->z[j-1][k+1] + psf->z[j-1][k-1]) / 4.;
	    };
	};
    };
}
