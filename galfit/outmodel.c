#include <stdio.h>
#include <string.h>
#include "nrutil.h"
#include "structs.h"

#define NEW 0
#define ADD 1

void writefits (char *phead, struct image *img, char *object, int add);
void writeheader (struct image *img, struct fitpars *fpar, double chisq,
				int ndof, int hdu, int offx, int offy);
void subimgs (struct image *data, struct image *model, struct image *sub);
void mkmodel (struct image *model, struct derivs *df, struct fitpars *fpar,
					struct convpars *cpar, int output);

void outmodel (struct fitpars *fpar, struct convpars *cpar, struct image *d, 
	 struct image *model, struct derivs *df, double chisq, int ndof,
							int offx, int offy)
{
    struct fitpars *fptr = fpar;
    struct fitpars *tmpptr;
    struct derivs *dfptr = df;
    struct derivs *tmpdf;
    struct image sub;
    extern struct inpars input;
    char obj[25], infile[80];
    extern unsigned long *pos, *cpos;

    /****************************************************\
    *  Output original image into the output image block *
    \****************************************************/

    /* Output original image into the output image block */

    if (input.create != 1) {
        strcpy (infile, d->name);
	strcpy (d->name, input.output);
        writefits (infile, d, "original image", NEW);
    };

    /**********************************************************\
    *  Create and output final models into 2nd image in block  *
    \**********************************************************/

    sprintf (obj, "model");
    strcpy (model->name, input.output);
    mkmodel (model, df, fpar, cpar, 1);
    if (input.create != 1) {
        writefits (infile, model, obj, ADD);
	writeheader (model, fpar, chisq, ndof, 3, offx, offy);

    };

    /****************************************************\
    *  Output residual image into third image in block   *
    \****************************************************/

    if (input.create != 1) {
        sub = *d;
	sub.z = matrix (1, d->naxes[2], 1, d->naxes[1]);
        subimg (d, model, &sub);
        sprintf (obj, "residual map");
        writefits (infile, &sub, obj, ADD);
        free_matrix(sub.z, 1, d->naxes[2], 1, d->naxes[1]);
    } else {
        writefits (infile, model, obj, NEW);
	writeheader (model, fpar, chisq, ndof, 1, offx, offy);
    };

    /**********************************************************\
    *  Create and output final model components into 4th +      *
    *  image in blocks                                          *
    \**********************************************************/
    do
      {
	sprintf (obj, fptr->objtype);
	strcpy (model->name, input.output);
	tmpptr = fptr->next;
	fptr->next = NULL;
	tmpdf = dfptr->next;
	dfptr->next = NULL;
	printf("making model %s\n", fptr);
	mkmodel (model, dfptr, fptr, cpar, 1);
	if (input.create != 1) {
	  writefits (infile, model, obj, ADD);
	  writeheader (model, fptr, chisq, ndof, 3, offx, offy);
	  };
	fptr->next = tmpptr;
	dfptr->next = tmpdf;
	fptr = tmpptr;
	dfptr = tmpdf;
      }while( fptr != NULL);
}
