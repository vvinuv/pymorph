#include <stdio.h>
#include <string.h>
#include "math.h"
#include "fitsio.h"
#include "structs.h"
#include "nrutil.h"

#define STRLEN 100
void template (FILE *, struct inpars *, struct fitpars *fpar);
void getline (char buffer[], int n);
void menu (struct inpars *input, struct fitpars *fpar, FILE *ftype, 
			float xoffset, float yoffset, float chisq, int ndof);

int totnobj=0;

void read_input (struct inpars *input, struct fitpars *fpar)
{
    void checkpar (char *, struct inpars *, struct fitpars *fpar, int *objnum);

    FILE *pfile;
    char string[STRLEN];
    int i, quit=0, objnum;

    strcpy (fpar->objtype, "none");
    fpar->next = NULL;
    
    /*************************************************************\
    *  Read in a file that has all the parameters in the needed   *
    *  If there's no such file, then the user must enter them by  *
    *  hand via a menu.                                           *
    \*************************************************************/

    if (!strncmp(input->initparfile, "\0", 1)) {
        printf ("Enter template file name: ");
        getline (input->initparfile, STRLEN);
	input->interact = 1;
    };
    if (strncmp(input->initparfile, "\0", 1)) {
	pfile = fopen (input->initparfile, "r");
	if (pfile == (FILE *) 0) {
	    printf ("\nError reading template parameter file....\n");
	    printf ("Enter in the parameters manually. \n");
            strcpy (input->initparfile, "none");
	} else {
	    template (pfile, input, fpar);
	    fclose (pfile);
	}
    };

    /* Loop until all the required parameters are read in */

    objnum = totnobj;
    while ( !quit ){
        menu (input, fpar, stdout, 0., 0., 0., 0);
        sprintf (string, "");
        while (strncmp(string, "r", 1) !=0 && !quit && input->interact == 1){
            printf ("Currently modifying object number:  %d \n", objnum);
            printf ("Enter item, initial value(s), fit switch(es) ==> ");
            getline (string, STRLEN);
            checkpar (string, input, fpar, &objnum);
	    if (strncmp(string, "q", 1) == 0 || input->interact == 0) quit = 1;
	};
	if (input->interact == 0) quit= 1;
    };
    return;
}


/* ------------------------------------------------------------------------ *\
\* ------------------------------------------------------------------------ */


void checkpar (char *string, struct inpars *input, struct fitpars *fpar, 
							int *objnum)
{
    char com[STRLEN], *ptr, objtype[10];
    int intcom, delnum, box[5], i;
    struct fitpars *newobj, *ptrstart;
    extern struct sampars sample;
    FILE *pfile;

    void del_obj (struct fitpars *fpar, int delnum, struct sampars *sample);

    ptrstart = fpar;

    for (i=0; i < *objnum - 1; ++i) fpar = fpar->next;

    sscanf (string, " %s", com);

    /*****************************
    *    Convert to uppercase    *
    *****************************/

    intcom = 'A';
    ptr=com;
    while (intcom != '\0') {
        intcom = *ptr;
        if ( 'a' <= intcom && intcom <= 'z' )
            *ptr = *ptr + 'A' - 'a';
        ptr++;
    };

    if (strncmp (com, "A", 1)==0)
	sscanf (string, " %s %s", com, input->data);
    else if (strncmp (com, "B", 1)==0)
	sscanf (string, " %s %s", com, input->output);
    else if (strncmp (com, "C", 1)==0)
	sscanf (string, " %s %s", com, input->sigma);
    else if (strncmp (com, "D", 1)==0)
	sscanf (string, " %s %s %s", com, input->psf, input->kernel);
    else if (strncmp (com, "E", 1)==0)
	sscanf (string, " %s %d", com, &input->sampfac);
    else if (strncmp (com, "F", 1)==0) 
	sscanf (string, " %s %s", com, input->badpix);
    else if (strncmp (com, "G", 1)==0) 
	sscanf (string, " %s %s", com, &input->constraints);
    else if (strncmp (com, "H", 1)==0) {
	sscanf (string, " %s %d %d %d %d", com, &box[1], &box[2],
					     &box[3], &box[4]);
	sprintf (input->imgsect, "[%d:%d,%d:%d]", box[1], box[2], 
                                                            box[3], box[4]);
    } else if (strncmp (com, "I", 1)==0)
	sscanf (string, " %s %d %d", com, &input->convbox[1], 
							&input->convbox[2]);
    else if (strncmp (com, "J", 1)==0)
	sscanf (string, " %s %f", com, &input->magzpt);
    else if (strncmp (com, "K", 1)==0)
	sscanf (string, " %s %f %f", com, &input->d[1], &input->d[2]);
    else if (strncmp (com, "L", 1)==0)
	sscanf (string, " %s %d", com, &input->coerr);
    else if (strncmp (com, "M", 1)==0)
	sscanf (string, " %s %d", com, objnum);
    else if (strncmp (com, "O", 1)==0) {
	sscanf (string, " %s %s", com, &input->device);
        if (strncmp (input->device, "both", 4) != 0 &&
                             strncmp (input->device, "reg", 3) != 0 &&
                                strncmp (input->device, "cur", 3) != 0)
        strcpy (input->device, "regular");
    }
    else if (strncmp (com, "P", 1)==0)
	sscanf (string, " %s %d", com, &input->create);
    else if (strncmp (com, "S", 1)==0)
	sscanf (string, " %s %d", com, &input->interact);
    else if (strncmp (com, "X", 1)==0) {
	sscanf (string, " %s %d", com, &delnum);
        if (delnum <= totnobj) {
            fpar = ptrstart;
	    del_obj (fpar, delnum, &sample);
        };
	totnobj = sample.nobjs;
        if (*objnum > totnobj) *objnum = totnobj;
    }
    else if (strncmp (com, "Z", 1)==0)
	sscanf (string, " %s %d", com, &fpar->outtype);
    else if (strncmp (com, "0", 1)==0 || strncmp (com, "N", 1)==0) {

    /* Create a new link in the chain to hold parameters for a     */
    /* new object.                                                 */

	sscanf (string, " %s %s", com, objtype);

        while ( strncmp (objtype, "sersic", 6) != 0 && 
	        strncmp (objtype, "expdisk", 7) != 0 && 
	        strncmp (objtype, "nuker", 5) != 0 &&
		strncmp (objtype, "gaussian", 8) != 0 &&
		strncmp (objtype, "psf", 3) != 0 &&
		strncmp (objtype, "moffat", 6) != 0 &&
		strncmp (objtype, "devauc", 6) != 0 &&
		strncmp (objtype, "king", 4) != 0 &&
		strncmp (objtype, "sky", 3) != 0)       {
            printf ("\nEnter new object name: ");
  	    getline (objtype, 10);
	};

        while (fpar->next != NULL) fpar = fpar->next;

        if (strncmp(fpar->objtype, "none", 4) != 0) {
            newobj = (struct fitpars *) malloc(sizeof (struct fitpars));
	    fpar->next = newobj;
            fpar = fpar->next;

	    for (i = 0; i <= NPARS; i++) {    /* Initialize parameters */
		fpar->a[i] = 0.;
		fpar->ia[i] = 0;
	    };
            fpar->a[8] = 1.;     /* Axis ratio */
	    fpar->outtype = 0;
	};

	totnobj++;
	fpar->next = NULL;
        strcpy (fpar->objtype, objtype);
	sample.nobjs++;
	printf ("\n");
    } 
    else if (strncmp (com, "10", 2)==0)
	sscanf (string, " %s %f %d", com, &fpar->a[10], &fpar->ia[10]);
    else if (strncmp (com, "1", 1)==0 )
	if (strncmp (fpar->objtype, "sky", 3) == 0) 
	    sscanf (string, " %s %f %d", com, &fpar->a[1], &fpar->ia[1]);
	else
	    sscanf (string, " %s %f %f %d %d", com, &fpar->a[1], &fpar->a[2], 
			                        &fpar->ia[1], &fpar->ia[2]);
    else if (strncmp (com, "2", 1)==0 ) {
	if (strncmp (fpar->objtype, "sky", 3) == 0) 
	    sscanf (string, " %s %f %d", com, &fpar->a[2], &fpar->ia[2]);
    } 
    else if (strncmp (com, "3", 1)==0)
	sscanf (string, " %s %f %d", com, &fpar->a[3], &fpar->ia[3]);
    else if (strncmp (com, "4", 1)==0)
	sscanf (string, " %s %f %d", com, &fpar->a[4], &fpar->ia[4]);
    else if (strncmp (com, "5", 1)==0)
	sscanf (string, " %s %f %d", com, &fpar->a[5], &fpar->ia[5]);
    else if (strncmp (com, "6", 1)==0)
	sscanf (string, " %s %f %d", com, &fpar->a[6], &fpar->ia[6]);
    else if (strncmp (com, "7", 1)==0)
	sscanf (string, " %s %f %d", com, &fpar->a[7], &fpar->ia[7]);
    else if (strncmp (com, "8", 1)==0) {
	sscanf (string, " %s %f %d", com, &fpar->a[8], &fpar->ia[8]);
	if (fpar->a[8] == 1.)
	    fpar->a[8] = 0.999999;
    } else if (strncmp (com, "9", 1)==0)
	sscanf (string, " %s %f %d", com, &fpar->a[9], &fpar->ia[9]);

    else if (strncmp (com, "T", 1)==0) {
        printf ("Enter template file name: ");
        getline (string, STRLEN);
        if (strncmp(string, "\0", 1)) {
	    pfile = fopen (string, "r");
	    if (pfile == (FILE *) 0) {
	        printf ("\nError reading template parameter file....\n");
	        printf ("Enter in the parameters manually. \n");
	    } else {
		fpar = ptrstart;
	        template (pfile, input, fpar);
	        fclose (pfile);
		*objnum = totnobj;
	    };
        };
    }

    else if (strncmp (com, "#", 1)==0 || strncmp (com, "Q", 1) ==0 || 
	strncmp(com,"=",1)==0 || strncmp(com,"R",1)==0)    /* Do nothing */;
    else {
	printf ("\n");
	printf ("Your options are: \n");
	printf ("  [ ] Enter parameter number/letter followed by value(s)\n");
        printf ("  [M number] Modify object number parameters \n");
        printf ("  [N] Add new object \n");
	printf ("  [R] Redisplay the parameter menu \n");
	printf ("  [T] Read in initial parameter template file \n");
	printf ("  [X number] Delete object number \n");

	printf ("  [Q] Quit menu and go on to fitting \n");
	printf ("\n");
    };
    return;
}


/* ------------------------------------------------------------------------ *\
\* ------------------------------------------------------------------------ */


void template (FILE *pfile, struct inpars  *input, struct fitpars *fpar)
{
    void checkpar (char *, struct inpars *, struct fitpars *fpar, int *objnum);
    char string[STRLEN];

    while ( ! feof (pfile) ) {
        if ( fscanf (pfile, " %[^\n]", string) > 0 )
             checkpar (string, input, fpar, &totnobj);
    };
}


/* ------------------------------------------------------------------------ *\
\* ------------------------------------------------------------------------ */

void del_obj (struct fitpars *fpar, int delnum, struct sampars *sample) 
{

    int i=0, j;
    struct fitpars *fptr;

    if ( delnum == 1 ) {
        for (i = 1; i <= sample->nobjs-1; i++) {
            strcpy (fpar->objtype, fpar->next->objtype);
            for (j = 1; j<= NPARS; j++) {
                fpar->ia[j] = fpar->next->ia[j];
                fpar->a[j] = fpar->next->a[j];
            };
            fpar->outtype = fpar->next->outtype;
            if (i != sample->nobjs-1)
	        fpar = fpar->next;
        };
        free (fpar->next);
        fpar->next = NULL;
        if (sample->nobjs == 1) strcpy (fpar->objtype, "none");
    } else {
        for (i = 1; i < delnum-1; i++) fpar=fpar->next;
        fptr = fpar->next;
        fpar->next = fptr->next;
        free (fptr);
    };

    --(sample->nobjs); 
}
