#include <stdio.h>
#include <curses.h>
#include "structs.h"

void outmenu (struct fitpars *fpar, char newout[], float xoffset, 
				float yoffset, float chisq, int ndof);

int keypoll (struct fitpars *fpar, int *ITERMAX, int lowx, int lowy,
						float chisq, int ndof)
		
{
    int quit = 0;
    extern int (*pfunc)(const char *, ...);
    extern char device[];
    extern struct inpars input;
    struct fitpars tfpar;
    char C, newout[11];
    float xoffset, yoffset;

    C = getch();

    if (C == 'n' || C == 'N') {
        pfunc ("\n");
        pfunc ("Enter a new maximum number of iterations ==> ");
        if (strncmp (device, "both", 4) ==0)
	    scanf (" %d", ITERMAX);
	else {
            refresh();
            nodelay (stdscr, FALSE);
	    scanw (" %d", ITERMAX);
            nodelay (stdscr, TRUE);
	};
    };

    if (C == 'o' || C == 'O') {
        xoffset = (float) lowx - 1.;
        yoffset = (float) lowy - 1.;
        outmenu (fpar, newout, xoffset, yoffset, chisq, ndof);
        pfunc ("\n");
        pfunc ("Outputting current parameters to %s.\n", newout);
    };

    if (C == 'p' || C == 'P') {
        pfunc ("\n");
        pfunc ("Paused....  Hit P to unpause.\n");
        noecho();
        while ((C = getch()) != 'p' && C != 'P')
	    ;
        pfunc ("Continuing....\n");
        echo();
        if (strncmp (device, "curses", 6) ==0)
	    refresh();
    };

    if (C == 'q' || C == 'Q') {
        quit = 1;
        pfunc ("\n\n");
        pfunc ("Quitting out early now....\n");
        pfunc ("The solution may not be fully optimized.\n\n");
        if (strncmp (device, "curses", 6) == 0)
            refresh ();
    };

    return (quit);
}
