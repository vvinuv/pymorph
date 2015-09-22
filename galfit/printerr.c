#include <stdio.h>
#include "fitsio.h"

void printerror(int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/


    if (status)
    {
       fits_report_error(stderr, status);    /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

