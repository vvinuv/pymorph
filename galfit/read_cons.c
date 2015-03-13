#include "structs.h"
#include <stdio.h>
#include <stdlib.h>
#include "string.h"

int findpos (char *haystack, char *needle);
int blankstr (char *instring);
void conv2lower (char *par);

int read_cons (struct inpars input, struct cons *constr)
{
    int parconv (char *string);

    FILE *cfile;
    char instring[80];
    char ops[] = "*/+-";
    char *p;
    char tok1[10], tok2[10], tok3[10], tok4[10], tok5[10], tok6[10];
    int op_pos, com_pos, ncons, pos;
    float temp;

    constr->next =  NULL;

    cfile = fopen (input.constraints, "r");
    if (cfile == (FILE *) 0)
	return (1);

    ncons = 0;
    while (fgets (instring, 80, cfile) != NULL) {
        /* Don't be fooled when the op symbol comes after the comment field */
        com_pos = findpos (instring, "#");
        op_pos = findpos (instring, ops);

        /***************************************************************\
        *  Figure out whether the constraint is to relate 2 parameters  *
        *  or to relate just a single one.                              *
        \***************************************************************/

	if (ncons >= 1 && !blankstr(instring)) {
	    constr->next = (struct cons *) malloc (sizeof (struct cons));
	    constr = constr->next;
	    constr->next = NULL;
	};

        if (com_pos != 1 && !blankstr(instring)) {
 	    ncons++;

	    /* Read in the first 6 tokens and parse them for an operator */

            if ((p = strtok (instring, " \t")) != NULL)
 	        strcpy (tok1, p);
	    else
		return (1);

	    if ((p = strtok (NULL, " \t")) != NULL)
	        strcpy (tok2, p);
	    else
		return (1);

	    if ((p = strtok (NULL, " \t")) != NULL)
	        strcpy (tok3, p);
	    else
		return (1);

            if ((p = strtok (NULL, " \t")) != NULL)
  	        strcpy (tok4, p);
	    else
		return (1);

            p = strtok (NULL, " \t");
	    if (p != NULL)
	 	strcpy (tok5, p);
            p = strtok (NULL, " \t");
	    if (p != NULL)
	    	strcpy (tok6, p);

	    pos = findpos (tok1, ops);     /* Find operator in first token */
	    if (pos != 0) {		   /* If found, dice it up.        */
	    	constr->op = tok1[pos-1];
	    	p = strtok (tok1, ops);
	 	constr->comp[1] = atoi(p);
	    	p = strtok (NULL, "#+/-* \t");
	    	if (p != NULL && !blankstr (p)) {
	            constr->comp[2] = atoi(p);
	            constr->par = parconv (tok2);
	    	    constr->cval[1] = atof (tok3);
	     	    constr->cval[2] = atof (tok4);
	    	} else {
		    constr->comp[2] = atoi(tok2);
		    constr->par = parconv (tok3);
          	    constr->cval[1] = atof (tok4);
	            constr->cval[2] = atof (tok5);
	        };
	    } else {                     /* Look for operator in 2nd token.  */
	        constr->comp[1] = atoi (tok1);
	        pos = findpos (tok2, ops);    

	        if (pos != 0) {	         /* Found operator -- dice it up.    */
	            constr->op = tok2[pos-1];
	            p = strtok (tok2, "#+/-* \t");
	            if (p != NULL && !blankstr (p)) {
	                constr->comp[2] = atoi(p);
	 	        constr->par = parconv (tok3);
	  	        constr->cval[1] = atof (tok4);
		        constr->cval[2] = atof (tok5);
		    } else {			
		        constr->comp[2] = atoi(tok3);
	   	        constr->par = parconv (tok4);
	  	        constr->cval[1] = atof (tok5);
		        constr->cval[2] = atof (tok6);
		    };
	    	};	    
	    };

	    if (pos == 0) {		/* No operator found. */
	        constr->op = 0;
	        constr->comp[1] = atoi (tok1);
	        constr->comp[2] = 0;
	        constr->par = parconv (tok2);
 	        constr->cval[1] = atof (tok3);
		if (findpos (tok4, "to") != 0) {
		    constr->op = 't';           /* Constrain by min max */
		    constr->cval[2] = atof (tok5); 
		} else
	            constr->cval[2] = atof (tok4);
	    };

	    /* Make sure the constraints are in the right order */
	    if (constr->cval[1] > constr->cval[2]) {
		temp = constr->cval[1];
		constr->cval[1] = constr->cval[2];	
		constr->cval[2] = temp;
	    };
        };
    };
    
    return (0);
}


/********************************************************\
*  Convert parameters from names like mag, re, etc.      *  
*  to parameter numbers.                                 *
\********************************************************/

int parconv (char *par)
{
   /*  Check to see whether the parameter is a name or a number  */
    if (*par >= 48 && *par <= 57)
        return (atoi(par));
    else {
        conv2lower (par);
        if (strncmp (par, "x", 1) == 0 || strncmp (par, "sky", 3) == 0) 
	    return (1);
        if (strncmp (par, "y", 1) == 0) return (2);
        if (strncmp (par, "mag", 3) == 0) return (3);
        if (strncmp (par, "re", 2) == 0 || strncmp (par, "rs", 2) == 0 ||
                                                strncmp (par, "rb", 2) == 0)
            return (4);
        if (strncmp (par, "n", 1) == 0 || strncmp (par, "alpha", 5) == 0)
            return (5);
        if (strncmp (par, "beta", 4) == 0) return (6);
        if (strncmp (par, "gamma", 5) == 0) return (7);
        if (strncmp (par, "q", 1) == 0 || strncmp(par, "b/a", 3) == 0) 
            return (8);
        if (strncmp (par, "pa", 2) == 0) return (9);
        if (strncmp (par, "c", 1) == 0) return (10);
    }

    return (-1);
}

