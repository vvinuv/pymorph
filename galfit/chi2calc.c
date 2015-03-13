#include <math.h>
#include "structs.h"

double chi2calc (struct fitpars *fpar, float **data, float **sig, 
				float **model, long naxes[], int *ndof)
{
    extern unsigned long *pos; 
    extern struct image *mask;
    int i, ix, iy, nfree;
    double chi2 = 0.;
    float dy;

    *ndof = 0;
    for (iy=1; iy<=naxes[2]; iy++) {
        for (ix=1; ix<=naxes[1]; ix++) {
            if (mask->z[iy][ix] < 1.) {
                dy = data[iy][ix] - model[iy][ix];
                chi2 += dy*dy /(sig[iy][ix] * sig[iy][ix]);
                (*ndof)++;
	    };
	};
    };

    nfree = 0;
    while (fpar != NULL) {
        for (i=1; i <= NPARS; i++)
            if (fpar->ia[i] == 1) nfree++;
	fpar = fpar->next;
    };
    (*ndof) -= nfree;

    return (chi2);
}
