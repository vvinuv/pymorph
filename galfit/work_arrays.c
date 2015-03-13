#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "nrutil.h"

void convregion (struct derivs *df, float a[], struct image *d, 
							struct convpars *cpar);

void create_work_arrays (struct fitpars *fpar, struct image *model, 
							struct derivs *df)
{
    int i;
    struct derivs *newdf;

    while (fpar != NULL) {
	for (i = 0; i<=NPARS; i++)
	    if (fpar->ia[i] == 1 || i==0) {

	        /*  Allocate space for the entire fitting region  */

	        df->naxes[1] = model->naxes[1];
	        df->naxes[2] = model->naxes[2];
	        df->dpm[i] = matrix (1, model->naxes[2], 1, model->naxes[1]);
	    };

	if (fpar->next != NULL) {
	    newdf = (struct derivs *) malloc ((size_t)(sizeof (struct derivs)));
	    if (!newdf) nrerror("Error in creating work array \n");
	    df->next = newdf;
	    df = df->next;
	};

	df->next = NULL;
	fpar = fpar->next;
    };

    if (model->naxes[2] > 0 && model->naxes[1] > 0)
        model->z = matrix (1, model->naxes[2], 1, model->naxes[1]);

}

/************************************************************************/

void create_conv_arrays (struct fitpars *fpar, struct convpars *cpar, 
					struct image *model, struct derivs *df)
{
    int i;
    struct derivs *newdf;

    while (fpar != NULL) {
	for (i = 0; i<=NPARS; i++) {
            if (strncmp(fpar->objtype, "sky", 3)!=0) {
	        if (fpar->ia[i] == 1 || i==0) {

 	            /*  Create images for convolution region only  */

	            convregion (df, fpar->a, model, cpar);
		    if (df->naxes[2] > 0 && df->naxes[1] > 0)
	                df->dpm[i] = matrix(1, df->naxes[2], 1, df->naxes[1]);
		}
	    } else {
	        df->naxes[1] = 0;
	        df->naxes[2] = 0;
	    };
	};

	if (fpar->next != NULL) {
	    newdf = (struct derivs *) malloc ((size_t)(sizeof (struct derivs)));
	    if (!newdf) nrerror("Error in creating work array \n");
	    df->next = newdf;
	    df = df->next;
	};

	df->next = NULL;
	fpar = fpar->next;
    };
}

/************************************************************************/

void annihilate_work_arrays (struct fitpars *fpar, struct image *model,
    struct derivs *df)
{
    int i;
    struct fitpars *fptr;
    struct derivs *dptr;
    extern struct sampars sample;

    fptr = fpar;
    dptr = df;

    while (fptr != NULL) {
        for (i=0; i <= NPARS; ++i) {
	    if (dptr->naxes[1] > 0 && dptr->naxes[2] > 0)
                if (fptr->ia[i]==1 || i==0) free_matrix (dptr->dpm[i], 1,
		    dptr->naxes[2], 1, dptr->naxes[1]);
	};

	fptr = fptr->next;
	dptr = dptr->next;
    };

    if (model->naxes[2] > 0 && model->naxes[1] > 0)
        free_matrix (model->z, 1, model->naxes[2], 1, model->naxes[1]);

}

/************************************************************************/

void delete_fitpars (struct fitpars *fpar)
{
    if (fpar != NULL) {
	delete_fitpars (fpar->next);
	free((struct fitpars *)fpar);
    };
}

/************************************************************************/

void delete_derivs (struct derivs *df)
{

    if (df != NULL) {
	delete_derivs (df->next);
	free((struct derivs *)df);
    };
}
