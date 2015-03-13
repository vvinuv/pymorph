#include "structs.h"
#include <stdlib.h>
#include "debug.h"

void sersic (float a[], int ia[], struct image *model, struct derivs *df);
void nuker (float a[], int ia[], struct image *model, struct derivs *df);
void exponential (float a[], int ia[], struct image *model, struct derivs *df);
void moffat (float a[], int ia[], struct image *model, struct derivs *df);
void gaussian (float a[], int ia[], struct image *model, struct derivs *df);
void king (float a[], int ia[], struct image *model, struct derivs *df);
void psfunc (float a[], int ia[], struct image *model, struct derivs *df);
void sky (float a[], int ia[], struct image *model, struct derivs *df);

void objselect (float a[], int ia[], char objtype[], struct image *model, 
    struct derivs *df)
{


    if (strncmp (objtype, "devauc", 6)==0 ||
				strncmp (objtype, "sersic", 6)==0)
	sersic (a, ia, model, df);

    if (strncmp (objtype, "nuker", 5)==0 )
	nuker (a, ia, model, df);

    if (strncmp (objtype, "expdisk", 7)==0 )
	exponential (a, ia, model, df);

    if (strncmp (objtype, "moffat", 6)==0 )
	moffat (a, ia, model, df);

    if (strncmp (objtype, "gaussian", 8)==0 )
	gaussian (a, ia, model, df);

    if (strncmp (objtype, "king", 4)==0 )
	king (a, ia, model, df);

    if (strncmp (objtype, "psf", 3)==0 )
	psfunc (a, ia, model, df);

    if (strncmp (objtype, "sky", 3)==0)
	sky (a, ia, model, df);


}


