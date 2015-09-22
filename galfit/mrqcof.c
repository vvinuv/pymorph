#include "nrutil.h"
#include "structs.h"
#include "debug.h"
#define NRANSI

void writefits(char *phead, struct image *img, char *object, int add);
void mkimg (float **array, long naxes[], char *outname, char *tname);
void mkmodel (struct image *model, struct derivs *df, 
	      struct fitpars *fpar, struct convpars *cpar, int output);

void mrqcof(float **y, float **sig, struct image *psf, float a[], 
	int ia[], int ma, double **alpha, double beta[], double *chisq,
	void (*funcs)(int x, int y, float *, float [], int [], int, 
	struct fitpars *), struct fitpars *fpar, struct convpars *cpar)
{
	extern struct image *mask;
        extern struct sampars sample;
        extern struct image model;
        extern unsigned long *pos;
        extern struct derivs df;

	struct fitpars *fptr;
	struct derivs *dptr;
        char name[80];

	int ix, iy, j,k,l,m,mfit=0;
	float ymod,wt,sig2i,dy,*dyda;

	dyda=vector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]==1) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}

	mkmodel (&model, &df, fpar, cpar, 0);

#if CKIMG
    sprintf (model.name, "gal+psf-model.fits");            /*  Output image */
    writefits ("test.fits", &model, "gal + PSF", 0);

    fptr = fpar;
    dptr = &df; 
    while (fptr != NULL) {
        for (j=0; j<= NPARS; j++){
            if (fptr->ia[j] == 1 || j==0) {
                sprintf (name, "deriv%d.fits", j);
                mkimg (dptr->dpm[j], dptr->naxes, name, "small.fits");
            };
        };
        fptr = fptr->next;
        dptr = dptr->next;
    };
#endif

	sample.nmask = 0;
	*chisq=0.0;
	for (iy=1; iy<= model.naxes[2]; iy++) {
	    for (ix=1;ix<=model.naxes[1];ix++) {
	        if (mask->z[iy][ix] < 1.) {        /* Calculate chi^2 only   */
		    (*funcs)(ix,iy,&ymod,dyda,ia,ma,fpar);  /* for unflagged */
		    sig2i=1.0/(sig[iy][ix]*sig[iy][ix]);           /* pixels */
		    dy=y[iy][ix]-ymod;
		    for (j=0,l=1;l<=ma;l++) {
			if (ia[l]==1) {
			    wt=dyda[l]*sig2i;
			    for (j++,k=0,m=1;m<=l;m++)
				if (ia[m]==1) alpha[j][++k] += wt*dyda[m]; 
			    beta[j] += dy*wt;
			};
		    };
		    *chisq += dy*dy*sig2i;
	        } else
		    sample.nmask++;	   /* count up the number of flagged pixels */
	    };
	};
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_vector(dyda,1,ma);
}
#undef NRANSI


/* --------------------------------------------------------------------- *\
\* --------------------------------------------------------------------- */


void funcs(int x, int y, float *ymodel, float dy[], int ia[], int ma,
						struct fitpars *fpar)
{
    int i, j, objnum=0;
    extern struct derivs df;
    extern struct image model;
    struct derivs *dptr;

    *ymodel = model.z[y][x];
    dptr = &df;
    while (fpar != NULL) {
	for (i=1; i <= NPARS; ++i){
	    j = objnum * NPARS + i;
            if (ia[j] == 1)
	        dy[j] = dptr->dpm[i][y][x];
	};
	dptr = dptr->next;
	fpar = fpar->next;
	objnum++;
    };

}

