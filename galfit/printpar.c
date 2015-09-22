#include <stdio.h>
#include <math.h>
#include <curses.h>
#include "structs.h"

float pa (float instrpa);

void printpar (int lowx, int lowy, struct fitpars *fpar)
{
    extern char device[];
    extern int (*pfunc)(const char *, ...);
    extern struct inpars input;
    extern float xskycent, yskycent;
    float x, y;
    int noff = 0;

    while (fpar != NULL) {
        x = fpar->a[1] + lowx - 1.;
	y = fpar->a[2] + lowy - 1.;
        fpar->a[9] = pa(fpar->a[9]);     /* make sure  -90 < PA <= 90 */
        
        if (strncmp(fpar->objtype, "sersic", 6)==0 || 
                                    strncmp(fpar->objtype, "moffat", 6)==0)
            pfunc (" %-9s: (%6.2f, %6.2f)  %7.2f %8.2f %7.2f  %6.2f  %6.2f %6.2f\n",
		fpar->objtype, x, y, fpar->a[3], fpar->a[4], fpar->a[5],
				 fpar->a[8], fpar->a[9], fpar->a[10]);

        if (strncmp(fpar->objtype, "devauc", 6)==0 || 
				strncmp(fpar->objtype, "expdisk", 7)==0 ||
				strncmp(fpar->objtype, "gaussian", 8)==0)
            pfunc (" %-9s: (%6.2f, %6.2f)  %7.2f %8.2f    ----  %6.2f  %6.2f %6.2f\n",
		fpar->objtype, x, y, fpar->a[3], fpar->a[4], 
				fpar->a[8], fpar->a[9], fpar->a[10]);

        if (strncmp(fpar->objtype, "psf", 3)==0)
            pfunc (" %-9s: (%6.2f, %6.2f)  %7.2f     ----    ----    ----    ----   ----\n",
		fpar->objtype, x, y, fpar->a[3]);

        if (strncmp(fpar->objtype, "nuker", 5)==0)
            pfunc (" %-9s: (%6.2f, %6.2f)  %5.2f %7.2f %5.2f %4.2f %4.2f  %5.2f %5.2f %5.2f\n",
	        fpar->objtype, x, y, fpar->a[3], fpar->a[4], fpar->a[5], 
                fpar->a[6], fpar->a[7], fpar->a[8], fpar->a[9], fpar->a[10]);

        if (strncmp(fpar->objtype, "king", 4)==0)
            pfunc (" %-9s: (%6.2f, %6.2f)  %7.2f %6.2f %6.2f %5.2f %5.2f %6.2f %6.2f\n",
	        fpar->objtype, x, y, fpar->a[3], fpar->a[4], fpar->a[5], 
                fpar->a[6], fpar->a[8], fpar->a[9], fpar->a[10]);

        if (strncmp(fpar->objtype, "penny", 5)==0)
            pfunc (" %-9s: (%6.2f, %6.2f)  %5.2f %7.2f  %6.2f\n",
				fpar->objtype, x, y, 
					fpar->a[3], fpar->a[4], fpar->a[5]);

	if (strncmp(fpar->objtype, "sky", 3) ==0)
	    pfunc (" sky      : [%6.2f, %6.2f]  %7.2f   %-5.2e   %-5.2e\n", 
		xskycent + lowx - 1, yskycent + lowy - 1, fpar->a[1], fpar->a[2], fpar->a[3]);

	fpar = fpar->next;
	noff += (NPARS+1);
    };

    if (strncmp (device, "curses", 6) == 0)
        refresh();
    return;
}

