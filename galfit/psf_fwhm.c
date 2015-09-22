/*****************************************************\
*  Roughly determine the FWHM of the input PSF image  *
\*****************************************************/

#include "structs.h"
#include "mymath.h"

float psf_fwhm (struct image psf)
{
    int i, j, rowc=1, colc=1;
    float lpix, hpix, lpflux, hpflux, hm, fwhm, max;

    max = psf.z[1][1];

    for (j=1; j <= psf.naxes[2]; j++) {
	for (i=1; i <= psf.naxes[1]; i++) {
	    if (psf.z[j][i] > max) {
    		colc = i;
		rowc = j;
		max = psf.z[j][i];
	    };
	};
    };

    hm = psf.z[rowc][colc] / 2.;

    i = 1;
    while (psf.z[rowc][i++] <= hm);
    lpix = i-2.;
    lpix += (hm - psf.z[rowc][(int)lpix]) / (psf.z[rowc][(int)lpix+1] - 
						psf.z[rowc][(int)lpix]);

    i = psf.naxes[1];
    while (psf.z[rowc][i--] <= hm);
    hpix = i+2.;
    hpix -= (hm - psf.z[rowc][(int)hpix]) / (psf.z[rowc][(int)hpix-1] - 
						psf.z[rowc][(int)hpix]);

    fwhm = hpix - lpix;

    i = 1;
    while (psf.z[i++][colc] <= hm);
    lpix = i-2.;
    lpix += (hm - psf.z[(int)lpix][colc]) / (psf.z[(int)lpix+1][colc] - 
						psf.z[(int)lpix][colc]);

    i = psf.naxes[2];
    while (psf.z[i--][colc] <= hm);
    hpix = i+2.;
    hpix -= (hm - psf.z[(int)hpix][colc]) / (psf.z[(int)hpix-1][colc] - 
						psf.z[(int)hpix][colc]);

    fwhm += (hpix - lpix);
    fwhm = fwhm / 2.;
    if (fwhm < 3.)
        fwhm -= 0.5;         /* On average this naive estimate is 0.5 pix off */
    return (fwhm);
}
