#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"
#include "structs.h"

void printerror(int status);

void groupfits(char *ename, char *gname) {

    fitsfile *gfptr, *mfptr;
    int status, ngroup, nhdu, i, *intnull;
    char *charnull;

    if (fits_open_file (&mfptr, ename, READONLY, &status))
	printerror (status);

    if (fits_get_num_hdus (mfptr, &nhdu, &status))
	printerror (status);

    printf ("%s %s %d\n", ename, gname, nhdu);

    intnull = (int *)0;
    charnull = (char *)0;

    if (fits_create_file (&gfptr, gname, &status))
	printerror (status);
    if (fits_create_group (gfptr, charnull, GT_ID_REF, &status))
	printerror (status);
    if (fits_add_group_member (gfptr, mfptr, 1, &status))
	printerror (status);

    for (i=1; i < nhdu; ++i) {
        if (fits_add_group_member (gfptr, mfptr, 1, &status))
	    printerror (status);
    }

    fits_close_file (gfptr, &status);
    fits_close_file (mfptr, &status);
    return;
}
