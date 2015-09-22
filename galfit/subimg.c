#include "structs.h"

void subimg (struct image *data, struct image *model, struct image *sub) 
{

    int ix, iy;

    for (iy = 1; iy <= model->naxes[2]; iy++)
        for (ix = 1; ix <= model->naxes[1]; ix++)
	    sub->z[iy][ix] = data->z[iy][ix] - model->z[iy][ix];

}
