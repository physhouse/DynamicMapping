#include "mapping.h"
#include "pointers.h"
#include "vector_CNv.h"
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

Vector_CNv::Vector_CNv(Mapping *map) : Pointers(map) {}
Vector_CNv::~Vector_CNv()
{
    cleanup();
}

void Vector_CNv::init()
{
    fg_num  = fg_atoms->fg_num;
    cg_num  = cg_sites->cg_num;

    // Initialize M
    CNv = new double *[3];
    for (int idim = 0; idim < 3; idim++)
    {
        CNv[idim] = new double[cg_num];
    }
}

void Vector_CNv::generate_CNv()
{
    double   **r      = fg_atoms->r;
    double   **R      = cg_sites->R;
    double   **C      = matrix_C->C;
    double ***  dw     = matrix_C->dw;
    double   **W      = matrix_C->W;
    double    *w_sum  = matrix_C->w_sum;
    double *v = fg_atoms->v;

    for (int I = 0; I < cg_num; I++)
    {
        for (int idim = 0; idim < 3; idim++)
        {
            CNv[idim][I] = 0.0;
        }
    }

    // To compute (C + N) * v, we start by looping over the fine-grained j,
    // then choose a neighboring CG site J, then a second neighboring site I,
    // and compute the influence of the Jj ownership relationship on each I
    // according to
    // (Nv)_IJj = R_Ij * c_Ij * (\delta_IJ - f_Jj) * -(d log w(R_Jj) / d R_Jj \cdot v_j)
    // First one gets v_j into cache, then calculates -(d log w(R_Jj) / d R_Jj \cdot v_j),
    // then calculates -R_Ij * c_Ij * (\delta_IJ - f_Jj) * (d log w(R_Jj) / d R_Jj \cdot v_j).
    // Cv is accumulated along the way as well.
    // See Eq. 17 of Overlap_Matrix_Physics_III.pdf

    //for (int j = 0; j < cg_sites->cg_num; j++) {
    for (int j = 0; j < fg_atoms->fg_num; j++) {
        // Access v_j
        double vj[3];
        for (int jdim = 0; jdim < 3; jdim++) 
        {
            vj[jdim] = v[jdim * fg_num + j];
        }
        // Calculate -(d log w(R_Jj) / d R_Jj \cdot v_j) for each neighbor J
        // This is -(1 / w_Jj) * (sum_jdim dw_Jjjdim * v_jjdim) 
        for (int cgneigh = 0; cgneigh < neighbor->numFgNeighbors[j]; cgneigh++)
        {
            int J = neighbor->fgList[j][cgneigh];
            double vj_dot_logprox_Jj = 0.0;
            for (int jdim = 0; jdim < 3; jdim++) 
            {
                vj_dot_logprox_Jj += dw[J][j][jdim] * vj[jdim];
            }
            vj_dot_logprox_Jj /= W[J][j];

            // Prepare f_Jj also.
            double f_Jj = W[J][j] / w_sum[j];
            
            // Calculate -R_Ij * c_Ij * (\delta_IJ - f_Jj) * (d log w(R_Jj) / d R_Jj \cdot v_j)
            // for each I.
            for (int cgneigh2 = 0; cgneigh2 < neighbor->numFgNeighbors[j]; cgneigh2++) 
            {
                int I = neighbor->fgList[j][cgneigh2];

                if (I == J) {
                    // Add the contribution from Nv. Delta_IJ = 1.
                    double N_IJj_scalar = C[I][j] * (1 - f_Jj) * vj_dot_logprox_Jj;
                    for (int Idim = 0; Idim < 3; Idim++) {
                        double coord = r[j][Idim];
                        wrap_coord_relative(coord, R[I][Idim], fg_atoms->L);
                        CNv[Idim][I] +=  (R[I][Idim] - coord) * N_IJj_scalar;
                    }
                    // Add the contribution from Cv
                    for (int Idim = 0; Idim < 3; Idim++) {
                        CNv[Idim][I] += C[I][j] * vj[Idim];
                    }
                } else {
                    // Add the contribution from Nv. Delta_IJ = 0.
                    double N_IJj_scalar = - C[I][j] * f_Jj * vj_dot_logprox_Jj;
                    for (int Idim = 0; Idim < 3; Idim++) {
                        double coord = r[j][Idim];
                        wrap_coord_relative(coord, R[I][Idim], fg_atoms->L);
                        CNv[Idim][I] +=  (R[I][Idim] - coord) * N_IJj_scalar;
                    }
                }
            }
        }
    }
}

void Vector_CNv::compute()
{
    generate_CNv();
}

void Vector_CNv::cleanup()
{
    for (int idim = 0; idim < 3; idim++)
    {
        delete[] CNv[idim];
    }
    delete[] CNv;
    printf("cleaning up CNv...\n");
}

