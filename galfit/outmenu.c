#include <stdlib.h>
#include <stdio.h>
#include "structs.h"

void menu (struct inpars *input, struct fitpars *fpar, FILE *logfile, 
			float xoffset, float yoffset, float chisq, int ndof);

char *outmenu (struct fitpars *fpar, char newout[], float xoffset, 
				float yoffset, float chisq, int ndof)
{
    int i = 1;
    FILE *logfile;
    extern struct inpars input;

    sprintf (newout, "galfit.01");
    while ( (logfile = fopen (newout, "r")) != (FILE *) 0 ) {
        sprintf (newout, "galfit.%02d", ++i);
        fclose (logfile);
    };
    logfile = fopen (newout, "w");
    input.create = 0;
    menu (&input, fpar, logfile, xoffset, yoffset, chisq, ndof);
    fclose(logfile);
}
