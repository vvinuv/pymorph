#include <math.h>
#include "nrutil.h"
#include "structs.h"

float sigcheck (struct image sigma)
{
    int i, j, flag0 = 0;
    float min=1e10;

    for (j=1; j <= sigma.naxes[2]; j++)
	for (i=1; i <= sigma.naxes[1]; i++) {
	    if (sigma.z[j][i] != 0.)
	        min = FMIN(min, fabs(sigma.z[j][i]));
	    else
		flag0 = 1;
	};

    if (flag0) {
        for (j=1; j <= sigma.naxes[2]; j++)
	    for (i=1; i <= sigma.naxes[1]; i++) {
	        if (sigma.z[j][i] == 0.)
		    sigma.z[j][i] = min;
		else if (sigma.z[j][i] <= 0.)
		    sigma.z[j][i] = fabs(sigma.z[j][i]);
	    };
	return (min);
    } else
        return (flag0);

}
