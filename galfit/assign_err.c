#include <math.h>
#include "structs.h"

void assign_err (struct fitpars *fpar, double *sig, float chi2nu)
{
    int i = 1, j;

    while (fpar != NULL) {
	for (j = 1; j <= 10; j++)
	    fpar->sig[j] = sig[i++] * sqrt(chi2nu);
	fpar = fpar->next;
    };
}
