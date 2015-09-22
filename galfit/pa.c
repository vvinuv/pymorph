#include <stdlib.h>

float pa (float instrpa)
{
    float newpa;
    int intpart;

    intpart = (int) (instrpa / 180.);
    newpa = instrpa - (intpart * 180.);

    if ( newpa <= -90.)
	newpa = newpa + 180.;
    else if ( newpa > 90. ) 
        newpa = newpa - 180.;

    return (newpa);
}
