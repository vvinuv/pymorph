#include <curses.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define NCHAR 200

#ifdef MACOSX
    #define SYSSTRING "\\ps -a > term-search"
    #define SLEEPSTR "sleep"
#elif defined(SUN)
    #define SYSSTRING "\\ps -A > term-search"
    #define SLEEPSTR "sleep"
#else
    #define SYSSTRING "\\ps -x > term-search"
    #define SLEEPSTR "sleep 31536000"
#endif

#ifdef PANTHER
    #define PREFIX "tty"
#else
    #define PREFIX ""
#endif

char *port_search ();

void initcurses (char device[])
{
    SCREEN *newxterm;
    FILE *xterm;
    int errn, port=0, ntries=0, i;
    char instring[200], *tty, dev[12];

    if (strncmp (device, "curses", 6) ==0) {
        initscr();
        nocbreak ();
        scrollok(stdscr, TRUE); 
        nodelay (stdscr, TRUE);
        echo ();
    } else if (strncmp (device, "both", 4) ==0) {

        /***********************************\
        *  Do a first pass to make sure no  *
        *  GALFIT window is already open    *
        \***********************************/

        tty = port_search ();

        /***********************************\
        *  If no previous GALFIT window,    *
        *  open one and search for its port *
        \***********************************/

        if (tty == NULL) {
            printf ("\n");
            errn = system ("xterm +sb -T \"GALFIT\" -geometry 21x6 -bg green1 -e sleep 31536000 &"); 
            while (tty == NULL && ntries++ < 20) {
                printf ("Trying to open up an interaction window....\n");
		system ("sleep 1");
                tty = port_search ();
            };
            if (tty == NULL) {
                printf ("\n");
                printf ("Interaction window not detected yet.  Quitting....\n");
                exit(0);
            };
	};

        strcpy (dev, tty);
        xterm = fopen (dev, "r+");
        newxterm = newterm("xterm", xterm, xterm); 

        refresh();
        nodelay (stdscr, TRUE);
        echo ();
        move (1, 0);

        printw ("  N = New max iter \n");
        printw ("  O = Output params \n");
        printw ("  P = Pause fit  \n");
        printw ("  Q = Quit ASAP \n");
        refresh (); 
        free (tty);
    };

    cbreak();

}

char *port_search ()
{
    int pid, oldport=0;
    FILE *ttyfile;
    char *tty, suffix[NCHAR], instring[NCHAR];

    system (SYSSTRING);

    ttyfile = fopen ("term-search", "r");

    tty = (char *)malloc (NCHAR);
    while (fgets(instring, NCHAR, ttyfile) != NULL) {
        if (strstr (instring, SLEEPSTR) != NULL) {
            sscanf (instring, "%d %s", &pid, &suffix);
            oldport = 1;
            sprintf (tty, "/dev/%s%s", PREFIX, suffix);
            break;
	};
    };

    fclose (ttyfile);
    system ("rm term-search");

    if (!oldport) {
        free (tty);
	return ((char *)NULL);
    } else
        return (tty);
}
