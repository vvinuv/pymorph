#include <math.h>
#include <stdlib.h>
#include "structs.h"

float pa (float instr_pa);

void assign_pars (struct inpars input, struct image *data, 
		  struct image *sigma, struct image *psf, 
		  struct image *badpix, struct fitpars *fpar, 
		  struct convpars *cpar)
{
    int i;
    
    strcpy (data->name, input.data);
    strcpy (sigma->name, input.sigma);
    strcpy (psf->name, input.psf);
    strcpy (badpix->name, input.badpix);

    while (fpar != NULL) {

        if (strncmp (fpar->objtype, "sersic", 6)==0) 
	{
            fpar->ia[6] = -1;		/*  code -1 means this parameter  */
            fpar->ia[7] = -1;		/*  is not being used.            */
	};

        if (strncmp (fpar->objtype, "devauc", 6)==0) 
        {
            fpar->a[5] = 4.;
	    fpar->ia[5] = 0;
	    fpar->ia[6] = -1;
	    fpar->ia[7] = -1;
	};

        if (strncmp (fpar->objtype, "expdisk", 7)==0) 
	{
	    fpar->ia[5] = -1;
	    fpar->ia[6] = -1;
	    fpar->ia[7] = -1;
	};

        if (strncmp (fpar->objtype, "gaussian", 8)==0) 
	{
            fpar->ia[5] = -1;
            fpar->ia[6] = -1;
            fpar->ia[7] = -1;
	};

        if (strncmp (fpar->objtype, "psf", 8)==0) 
	{
            fpar->ia[4] = -1;
            fpar->ia[5] = -1;
            fpar->ia[6] = -1;
            fpar->ia[7] = -1;
            fpar->ia[8] = -1;
            fpar->ia[9] = -1;
            fpar->ia[10] = -1;
	};

        if (strncmp (fpar->objtype, "king", 6)==0) 
	{
	    if (fpar->a[6] == 0.) {
		fpar->a[6] = 2.;
		fpar->ia[6] = 0;
	    };
            fpar->ia[7] = -1;
	};

        if (strncmp (fpar->objtype, "moffat", 6)==0) 
	{
            fpar->ia[6] = -1;
            fpar->ia[7] = -1;
	};

        if ( fpar->a[8] > 1.) {
            fpar->a[8] = 1./ fpar->a[8];
            fpar->a[9] = fpar->a[9] + 90.;
            fpar->a[4] = fpar->a[4] / fpar->a[8];
        };
        fpar->a[9] = pa (fpar->a[9]);

        /* Axis ratio can't be too small */
        if (fpar->ia[8] < 0.01) fpar->ia[8] = 0.01;

        if (strncmp (fpar->objtype, "sky", 3)==0) 
        {
	    fpar->ia[4] = -1;
	    fpar->ia[5] = -1;
	    fpar->ia[6] = -1;
	    fpar->ia[7] = -1;
	    fpar->ia[8] = -1;
	    fpar->ia[9] = -1;
	};

        for (i = 1; i<= NPARS; i++)
            if (fpar->ia[i] != 1 && fpar->ia[i] != 0 && fpar->ia[i] != -1)
                fpar->ia[i] = 0;

	fpar = fpar->next;
    };

    cpar->boxsz[1] = input.convbox[1];
    cpar->boxsz[2] = input.convbox[2];
    cpar->sampfac = input.sampfac;
 
    data->magzpt = input.magzpt;
    data->dp[1] = input.d[1];
    data->dp[2] = input.d[2];
    data->muzpt = data->magzpt + 2.5 * log10(data->dp[1] * data->dp[2]);
}
