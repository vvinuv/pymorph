#include <ctype.h>

/*****************************************************************\
*  Return the location of the character in "haystack" that first  *
*  matches a character in "needle".  For example:		  *
*     findpos ("abcde", "ce")					  *
#  will return the value 3.					  *
\*****************************************************************/

int findpos (char *haystack, char *needle)
{
    char *pos;
    int npos;

    npos = 0;
    while (*haystack != '\0') {
        npos++;
        pos = needle;
        while (*pos != '\0') {
            if (*pos == *haystack)
                return (npos);
            ++pos;
        };
        ++haystack;
    };
    return (0);
}


/************************************************************\
*  Check to see if the string is blank or comment.  If it's  *
*  either, return true (1), else return false (0).           *
\************************************************************/

int blankstr (char *instring)
{
    while (*instring != '\0' && *instring != '#') {
        if (*instring != ' ' && *instring != '\t' && *instring != '\n')
            return (0);
        ++instring;
    };
    return (1);
}


/**************************************\
*  Convert a string to all lower case  *
\**************************************/

void conv2lower (char *par)
{
    while (*par != '\0') {
        *par = tolower(*par);
        ++par;
    }
}
