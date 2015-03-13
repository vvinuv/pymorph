#include "structs.h"

void addsky (unsigned long pos[], float sky, struct image *model)
{
    unsigned long i, j;

    for (j=1; j <= model->naxes[2]; j++)
        for (i=1; i <= model->naxes[1]; i++)
	    model->z[j][i] += sky;
}

