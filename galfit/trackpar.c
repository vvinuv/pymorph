
void trackpar (int **track, int ia[], int *mfit, int ma)
{
    int j, k, l;

    for (*mfit=0,j=1;j<=ma;j++) {
        if (ia[j]==1) {
            (*mfit)++;
            if (j%10 != 0) {
                k = (int)(j/10)+1;
                l = j%10;
            } else {
                k = j/10;
                l = 10;
            };
            track[l][k] = *mfit;
        };
    };
}
