//2016-04-04
#include "mapping.h"
#include "pointers.h"
#include "matrix_N.h"
#include "matrix_C.h"
#include "cg_sites.h"
#include "fg_atoms.h"
#include "geom.h"
#include "neighbor.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Generating the N matrix, no need to cache the dcdr matrices

Matrix_N::Matrix_N(Mapping *map) : Pointers(map) {}
Matrix_N::~Matrix_N()
{
    cleanup();
}

void Matrix_N::init()
{
    fg_num  = fg_atoms->fg_num;
    cg_num  = cg_sites->cg_num;

    // Initialize M
    N = new double ***[3];
    for (int idim = 0; idim < 3; idim++)
    {
        N[idim] = new double **[3];
        for (int jdim = 0; jdim < 3; jdim++)
        {
            N[idim][jdim] = new double*[cg_num];
            for (int I = 0; I < cg_num; I++)
            {
                N[idim][jdim][I] = new double[fg_num];
                for (int j = 0; j < fg_num; j++) N[idim][jdim][I][j] = 0.0;
            }
        }
    }
}

void Matrix_N::generate_N()
{
    double   **r      = fg_atoms->r;
    double   **R      = cg_sites->R;
    double   **C      = matrix_C->C;
    double ***  dw     = matrix_C->dw;
    double   **W      = matrix_C->W;
    double    *w_sum  = matrix_C->w_sum;
    double   **dw_sum = matrix_C->dw_sum;

    // Notation dcdr[dim][I][j][i] = /partial C_I_i /partial r_dim_j
    for (int I = 0; I < cg_num; I++)
    {
        for (int j = 0; j < fg_atoms->fg_num; j++)
        {
            for (int idim = 0; idim < 3; idim++)
                for (int jdim = 0; jdim < 3; jdim++)
                    N[idim][jdim][I][j] = 0.0;
        }

        for (int jj = 0; jj < neighbor->numNeighbors[I]; jj++)
        {
            int j = neighbor->list[I][jj];

            // Note:: dC_I_i/dr_j is not a function of i!
            double dcdx = dw[I][j][0] / W[I][j] - dw_sum[j][0] / w_sum[j];
            double dcdy = dw[I][j][1] / W[I][j] - dw_sum[j][1] / w_sum[j];
            double dcdz = dw[I][j][2] / W[I][j] - dw_sum[j][2] / w_sum[j];

            //printf("%e %e %e\n", dcdx,dcdy,dcdz);

            for (int ii = 0; ii < neighbor->numNeighbors[I]; ii++)
            {
                int i = neighbor->list[I][ii];

                if (i != j)
                {
                    for (int jdim = 0; jdim < 3; jdim++)
                    {
                        /* Notation:: N[d1][d2][I][j], d1->dimension of CG particles, d2->dimension of FG particles
                         * N[d1][d2][I][j] = sumover(i){r_i_d1 * dC_I_i/dr_j_d2} */
                        double coord = r[i][jdim];
                        wrap_coord_relative(coord, R[I][jdim], fg_atoms->L);
                        N[jdim][0][I][j] +=  C[I][i] * C[I][j] * dcdx * (coord - R[I][jdim]);
                        N[jdim][1][I][j] +=  C[I][i] * C[I][j] * dcdy * (coord - R[I][jdim]);
                        N[jdim][2][I][j] +=  C[I][i] * C[I][j] * dcdz * (coord - R[I][jdim]);
                        //N[jdim][0][I][j] +=  C[I][i] * C[I][j] * dcdx * r[i][jdim];
                        //N[jdim][1][I][j] +=  C[I][i] * C[I][j] * dcdy * r[i][jdim];
                        //N[jdim][2][I][j] +=  C[I][i] * C[I][j] * dcdz * r[i][jdim];
                    }
                }
                else
                {
                    for (int jdim = 0; jdim < 3; jdim++)
                    {
                        /* Notation:: N[d1][d2][I][j], d1->dimension of CG particles, d2->dimension of FG particles
                         * N[d1][d2][I][j] = sumover(i){r_i_d1 * dC_I_i/dr_j_d2} */
                        double coord = r[i][jdim];
                        wrap_coord_relative(coord, R[I][jdim], fg_atoms->L);
                        N[jdim][0][I][j] += (C[I][i] - 1.0) * C[I][j] * dcdx * (coord - R[I][jdim]);
                        N[jdim][1][I][j] += (C[I][i] - 1.0) * C[I][j] * dcdy * (coord - R[I][jdim]);
                        N[jdim][2][I][j] += (C[I][i] - 1.0) * C[I][j] * dcdz * (coord - R[I][jdim]);
                        //N[jdim][0][I][j] +=  (C[I][i] - 1.0) * C[I][j] * dcdx * r[i][jdim];
                        //N[jdim][1][I][j] +=  (C[I][i] - 1.0) * C[I][j] * dcdy * r[i][jdim];
                        //N[jdim][2][I][j] +=  (C[I][i] - 1.0) * C[I][j] * dcdz * r[i][jdim];
                    }
                }
            }
            //printf("%12.8lf\t", N[0][0][I][j]);
        }
        //printf("\n");
    }
}

void Matrix_N::compute()
{
    generate_N();
}

void Matrix_N::cleanup()
{
    for (int idim = 0; idim < 3; idim++)
    {
        for (int jdim = 0; jdim < 3; jdim++)
        {
            for (int I = 0; I < cg_num; I++)
            {
                delete[] N[idim][jdim][I];
            }
            delete[] N[idim][jdim];
        }
        delete[] N[idim];
    }
    delete[] N;
    printf("cleaning up N...\n");
}

