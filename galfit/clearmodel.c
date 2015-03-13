#include <stdio.h>
#include "structs.h"

void clearmodel (struct image *model, struct derivs *df, struct fitpars *fpar)
{
    long i, j;
    int k;
    struct derivs *dfptr;
    struct fitpars *fptr;

    for (j=1; j <= model->naxes[2]; j++) {
        for (i=1; i <= model->naxes[1]; i++)
            model->z[j][i] = 0.;
    };

    dfptr = df;
    fptr = fpar;
    while (dfptr != NULL) {
        for (k=0; k<=10; k++) {
            if (fptr->ia[k] == 1 || k==0) {
                for (j=1; j <= dfptr->naxes[2]; j++) {
                    for (i=1; i <= dfptr->naxes[1]; i++)
                        dfptr->dpm[k][j][i] = 0.;
                };
            };
        };
        dfptr = dfptr->next;
        fptr = fptr->next;
    };
}

