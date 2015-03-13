/**************************************************************\
*  This subroutine Gaussian smooths an image given a Gaussian  *
*  sigma and a smoothing radius.                               *
\**************************************************************/

#include <math.h>
#include "nrutil.h"
#include "mymath.h"

void smooth (float **in, float **out, long naxes[], float sigma, float radius)
{
    int i, j, k, l, xmin, ymin, xmax, ymax;
    float sumwt, sumflux, wt, rad;

    for (j=1; j <= naxes[2]; j++) {
	for (i=1; i <= naxes[1]; i++) {
	    xmin = IMAX (1, i - NINT(radius));
	    xmax = IMIN (i+NINT(radius), naxes[1]);
	    ymin = IMAX (1, j - NINT(radius));
	    ymax = IMIN (j+NINT(radius), naxes[2]);
	    sumwt = 0.;
	    sumflux = 0.;
	    for (l=ymin; l <= ymax; l++) {
		for (k=xmin; k <= xmax; k++) {
		    rad = sqrt((k-i)*(k-i)+(l-j)*(l-j));
		    if (rad <= radius) {
		        wt = exp (-rad*rad/2/sigma/sigma);
			sumwt += wt;
			sumflux += in[l][k] * wt;
		    };
		};	    
	    };
	    out[j][i] = sumflux / sumwt;
	};
    };
}
