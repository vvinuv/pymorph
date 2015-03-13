#include <string.h>
#include "structs.h"

/*******************************\
*  Copy the parameter structs   *
\*******************************/

void copy_struct (struct fitpars *to, struct fitpars *from)
{
    int i;

    while (from != NULL) {
        for (i = 1; i <= NPARS; i++) {
            strncpy (to->objtype, from->objtype, 8);
            to->a[i] = from->a[i];
            to->ia[i] = from->ia[i];
        };
        to = to->next;
        from = from->next;
    };
}
