
void copymat (float **fmat, double **dmat, long naxes[], int op)
{
    int i, j;

    for (j=1; j <= naxes[2]; j++) {
	for (i=1; i <= naxes[1] * 2; i++) {
	    if (op == 1)
	        dmat[j-1][i-1] = fmat[j][i];
	    else if (op == -1)
		fmat[j][i] = dmat[j-1][i-1];
	};
    };
}
