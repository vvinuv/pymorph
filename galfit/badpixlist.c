#include <stdio.h>
#include "structs.h"
#include "nrutil.h"

#define STRLEN 80

int badpixlist (struct image *badpix)
{
    FILE *filename;
    
    int ix, iy;
    float x, y;
    char string[STRLEN];

    badpix->z = matrix (1, badpix->naxes[2], 1, badpix->naxes[1]);

    for (iy = 1; iy <= badpix->naxes[2]; iy++) {
        for (ix = 1; ix <= badpix->naxes[1]; ix++) {
            badpix->z[iy][ix] = 0.;
        };
    };

    if ((filename = fopen (badpix->name, "r")) == (FILE *)0)
	return (1);

    while (fscanf (filename, " %[^\n]", string) != EOF) {
	if (strncmp (string, "#", 1) != 0) {
	    sscanf (string, " %f %f", &x, &y);
	    ix = (unsigned long) x;
	    iy = (unsigned long) y;
	    if ( ix >= 1 && ix <= badpix->naxes[1] && iy >= 1 && 
						iy <= badpix->naxes[2] )
	        badpix->z[iy][ix] = 1.;
	};
    };
    fclose (filename);
    return (0);
}
