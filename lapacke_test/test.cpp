#include <stdio.h>
#include "lapacke.h"

int main()
{
    /*int N = 4;
    double A[16] = { 1, 2 , 3 , 1,
         4, 2 , 0 , 2,
        -2, 0 ,-1 , 2,
         3, 4 , 2 ,-3};
    double B[8] = { 6, 2 , 1 , 8,
        1, 2 , 3 , 4};
    int ipiv[4];
    int n = N;
    int nrhs = 2;
    int lda = N;
    int ldb = 2;*/

    int N = 3;
    double A[9] = { 1, 3, 5,
                    2, 4, 7,
                    -3, 2, 5
                  };
    double B[3] = { 2, -1, -5};
    int ipiv[3];
    int n = N;
    int nrhs = 1;
    int lda = N;
    int ldb = 1;

    int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, lda, ipiv, B, ldb);
    printf("info:%d\n", info);
    if (info == 0)
    {
        int i = 0;
        int j = 0;
        for (j = 0; j < nrhs; j++)
        {
            printf("x%d\n", j);
            for (i = 0; i < N; i++)
                printf("%.6g \t", B[i + j * N]);
            printf("\n");
        }
    }

    return 0;
}
