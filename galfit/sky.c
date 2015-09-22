#include "structs.h"
#include "mymath.h"

void sky (float *a, int *ia, struct image *model, struct derivs *df)
{
    int i, j;
    extern float xskycent, yskycent;

    for (j=1; j <= model->naxes[2]; j++) {
        for (i=1; i <= model->naxes[1]; i++) {
	    df->dpm[0][j][i] += (a[1] + (i-xskycent) * a[2] + 
							(j-yskycent) * a[3]);
	    if (ia[1] == 1) df->dpm[1][j][i] = 1.;
	    if (ia[2] == 1) df->dpm[2][j][i] = i-xskycent;
	    if (ia[3] == 1) df->dpm[3][j][i] = j-yskycent;
        };
    };
}

