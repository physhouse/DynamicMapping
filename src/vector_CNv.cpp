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
    double    **dw_sum  = matrix_C->dw_sum;
    double *v = fg_atoms->v;

    // Reset the (C + N) v vector to zero before adding
    // all of the elements together.
    for (int I = 0; I < cg_num; I++)
    {
        for (int idim = 0; idim < 3; idim++)
        {
            CNv[idim][I] = 0.0;
        }
    }

    // We calculate the N_Ij . v_j elements according to the formula
    // N_Ij . vj = (R_I - r_j) C_Ij [( dlog wIj dr_j - dlog Wj dr_j) . v_j]
    for (int j = 0; j < fg_atoms->fg_num; j++) {
        // Access v_j.
        double vj[3];
        for (int jdim = 0; jdim < 3; jdim++) vj[jdim] = v[jdim * fg_num + j];
        for (int cgneigh = 0; cgneigh < neighbor->numFgNeighbors[j]; cgneigh++)
        {
            int I = neighbor->fgList[j][cgneigh];
            double N_Ij_scalar = 0;
            // Calculate C_Ij [( dlog wIj dr_j - dlog Wj dr_j) . v_j]
            for (int jdim = 0; jdim < 3; jdim++)
            {
                N_Ij_scalar += C[I][j] * (dw[I][j][jdim] / W[I][j] - dw_sum[j][jdim] / w_sum[j]) * vj[jdim];
            }
            // Add the correct vectors to (C + N) v.
            for (int Idim = 0; Idim < 3; Idim++)
            {
                // Add the Cv vector term for this Ij pair.
                CNv[Idim][I] += C[I][j] * vj[Idim];
                // Add the Nv vector term for this Ij pair.
                double coord = r[j][Idim];
                wrap_coord_relative(coord, R[I][Idim], fg_atoms->L);
                CNv[Idim][I] +=  (R[I][Idim] - coord) * N_Ij_scalar;
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

