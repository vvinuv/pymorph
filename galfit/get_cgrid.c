#include "structs.h"
#include "nrutil.h"
#define NR_END 1
#define FREE_ARG char*

void bcucof(float y[], float y1[], float y2[], float y12[], float d1,
    float d2, float **c);
void get_derivs (struct image *psf, int i, int j, float **y1a, float **y2a,
    float **y12a, float yi[], float y1i[], float y2i[], float y12i[]);

struct bicubic_grid **get_cgrid (struct image *psf)
{
    struct bicubic_grid **new_cgrid(long nrl, long nrh, long ncl, long nch);

    int j, k, l, m;
    float **y1a, **y2a, **y12a;
    float yi[5], y1i[5], y2i[5], y12i[5];
    struct bicubic_grid **cg;

    y1a = matrix (0, psf->naxes[2]+1, 0, psf->naxes[1]+1);
    y2a = matrix (0, psf->naxes[2]+1, 0, psf->naxes[1]+1);
    y12a = matrix (0, psf->naxes[2]+1, 0, psf->naxes[1]+1);

    /*  Calculate first and second derivatives  */
    for (j = 0; j <= psf->naxes[2]+1; j++) {
	for (k = 0; k <= psf->naxes[1]+1; k++) {
	    if (j<=1 || k<=1 || j >= psf->naxes[2] || k >= psf->naxes[1]) {
		y1a[j][k] = 0.;
		y2a[j][k] = 0.;
		y12a[j][k] = 0.;
	    } else {
	        y1a[j][k] = (psf->z[j][k+1] - psf->z[j][k-1])/2.;
	        y2a[j][k] = (psf->z[j+1][k] - psf->z[j-1][k])/2.;
	        y12a[j][k] = (psf->z[j+1][k+1] - psf->z[j+1][k-1] - 
			      psf->z[j-1][k+1] + psf->z[j-1][k-1]) / 4.;
	    };
	};
    };

    /* Get bicubic interpolation coefficients for the PSF */

    cg = new_cgrid (0, psf->naxes[2]+1, 0, psf->naxes[1]+1);

    for (k = 0; k <= psf->naxes[2]+1; k++) {
	for (j = 0; j <= psf->naxes[1]+1; j++) {
	    cg[k][j].c = matrix (1, 4, 1, 4);
            get_derivs (psf, j, k, y1a, y2a, y12a, yi, y1i, y2i, y12i);
            bcucof (yi, y1i, y2i, y12i, 1., 1., cg[k][j].c);
        };
    };

    free_matrix (y1a, 0, psf->naxes[2]+1, 0, psf->naxes[1]+1);
    free_matrix (y2a, 0, psf->naxes[2]+1, 0, psf->naxes[1]+1);
    free_matrix (y12a, 0, psf->naxes[2]+1, 0, psf->naxes[1]+1);

    return (cg);
}


/*==========================================================================*/


struct bicubic_grid **new_cgrid(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        struct bicubic_grid **m;

        /* allocate pointers to rows */
        m=(struct bicubic_grid **) malloc((size_t)((nrow+NR_END)*
                                                sizeof(struct bicubic_grid *)));
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(struct bicubic_grid *) malloc((size_t)((nrow*ncol+NR_END)*
                                                sizeof(struct bicubic_grid)));
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

/*==========================================================================*/

void free_cgrid(struct bicubic_grid **m, long nrl, long nrh, long ncl, 
    long nch)
/* free a float matrix allocated by matrix() */
{
	int i, j;
	for (j=nrl; j <= nrh; j++)
	    for (i=ncl; i <= nch; i++)
	        free_matrix (m[j][i].c, 1, 4, 1, 4);

        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}


/*==========================================================================*/

void get_derivs (struct image *psf, int x, int y, float **y1a, float **y2a, 
    float **y12a, float yi[], float y1i[], float y2i[], float y12i[])
{
    yi[1] = yi[2] = yi[3] = yi[4] = 0.;
    y1i[1] = y1i[2] = y1i[3] = y1i[4] = 0.;
    y2i[1] = y2i[2] = y2i[3] = y2i[4] = 0.;
    y12i[1] = y12i[2] = y12i[3] = y12i[4] = 0.;

    if (x == 0) {
	if (y >= 1 && y <= psf->naxes[2]) {
            yi[2] = psf->z[y][x+1];
            y1i[2] = y1a[y][x+1];
            y2i[2] = y2a[y][x+1];
            y12i[2] = y12a[y][x+1];
	}
	if (y < psf->naxes[2]) {
            yi[3] = psf->z[y+1][x+1];
            y1i[3] = y1a[y+1][x+1];
            y2i[3] = y2a[y+1][x+1];
            y12i[3] = y12a[y+1][x+1];
	};
    } else if (y == 0) {
	if (x >= 1 && x <= psf->naxes[1]) {
            yi[4] = psf->z[y+1][x];
            y1i[4] = y1a[y+1][x];
            y2i[4] = y2a[y+1][x];
            y12i[4] = y12a[y+1][x];
	};
	if (x < psf->naxes[1]) {
            yi[3] = psf->z[y+1][x+1];
            y1i[3] = y1a[y+1][x+1];
            y2i[3] = y2a[y+1][x+1];
            y12i[3] = y12a[y+1][x+1];
 	};
    } else if (x == psf->naxes[1]) {
	if (y >= 1 && y <= psf->naxes[2]) {
            yi[1] = psf->z[y][x];
            y1i[1] = y1a[y][x];
            y2i[1] = y2a[y][x];
            y12i[1] = y12a[y][x];
	};
	if (y < psf->naxes[2]) {
            yi[4] = psf->z[y+1][x];
            y1i[4] = y1a[y+1][x];
            y2i[4] = y2a[y+1][x];
            y12i[4] = y12a[y+1][x];
	};
    } else if (y == psf->naxes[2]) {
	if (x >= 1 && x <= psf->naxes[1]) {
            yi[1] = psf->z[y][x];
            y1i[1] = y1a[y][x];
            y2i[1] = y2a[y][x];
            y12i[1] = y12a[y][x];
	};
	if (x < psf->naxes[1]) {
            yi[2] = psf->z[y][x+1];
            y1i[2] = y1a[y][x+1];
            y2i[2] = y2a[y][x+1];
            y12i[2] = y12a[y][x+1];
	};
    } else if (x >= 1 && x < psf->naxes[1] && y >= 1 && y < psf->naxes[2]) {
        yi[1] = psf->z[y][x];
        yi[2] = psf->z[y][x+1];
        yi[4] = psf->z[y+1][x];
        yi[3] = psf->z[y+1][x+1];

        y1i[1] = y1a[y][x];
        y1i[2] = y1a[y][x+1];
        y1i[3] = y1a[y+1][x+1];
        y1i[4] = y1a[y+1][x];

        y2i[1] = y2a[y][x];
        y2i[2] = y2a[y][x+1];
        y2i[3] = y2a[y+1][x+1];
        y2i[4] = y2a[y+1][x];

        y12i[1] = y12a[y][x];
        y12i[2] = y12a[y][x+1];
        y12i[3] = y12a[y+1][x+1];
        y12i[4] = y12a[y+1][x];
    };

}
