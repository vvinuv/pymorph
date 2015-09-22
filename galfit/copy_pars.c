#include "structs.h"

#define TO 0
#define FROM 1

void copy_pars (float a[], int ia[], struct fitpars *fpar, int dir)
{

    int i, k, j=0, count=0;

    while (fpar != NULL) {
        for (i = 1; i<= NPARS; i++) {
	    k = NPARS * count + i;
	    if (dir == TO) {
		fpar->a[i] = a[k];
	    } else if (dir == FROM ) {
		a[k] = fpar->a[i];
		ia[k] = fpar->ia[i];
	    };
	};
        fpar = fpar->next;
	count++;
    };

}

