#include "mapping.h"
#include "pointers.h"
#include "matrix_C.h"
#include "fg_atoms.h"
#include "cg_sites.h"
#include <cmath>

// constructor and destructor

Matrix_C::Matrix_C(Mapping *map) : Pointers(map) {}
Matrix_C::~Matrix_C()
{
    cleanup();
}

// Calculation functions
double Matrix_C::w_Ij(double *R, double *r)
{
    double distance = 0.0;
    for (int dim = 0; dim < 3; dim++)
    {
        double r_dim = R[dim] - r[dim];
        if (r_dim > 0.5 * L) distance += (r_dim - L) * (r_dim - L);
        else if (r_dim < -0.5 * L) distance += (r_dim + L) * (r_dim + L);
        else distance += r_dim * r_dim;
    }
    return weight(distance);
}

double Matrix_C::weight(double r)
{
    double r0 = r / rcut;
    return exp(- r0 * r0);
}

void Matrix_C::weight_deriv(double *R, double *r, double *dw_vec)
{
    double distance = 0.0;
    for (int dim = 0; dim < 3; dim++)
    {
        double r_dim = R[dim] - r[dim];
        if (r_dim > 0.5 * L) distance += (r_dim - L) * (r_dim - L);
        else if (r_dim < -0.5 * L) distance += (r_dim + L) * (r_dim + L);
        else distance += r_dim * r_dim;
    }

    double r0 = distance / rcut;
    double dwdr = - 2.0 * (r0 / rcut) * exp(-r0 * r0);

    for (int idim = 0; idim < 3; idim++)
    {
        double r_dim = R[idim] - r[idim];
        if (r_dim > 0.5 * L) r_dim -= L;
        else if (r_dim < -0.5 * L) r_dim += L;

        dw_vec[idim] = dwdr * r_dim / distance;
    }
    //printf("dw : %lf   %lf   %lf\n", dw_vec[0], dw_vec[1], dw_vec[2]);
}

// Initialization
void Matrix_C::init()
{
    fg_num = fg_atoms->fg_num;
    cg_num = cg_sites->cg_num;

    W = new double*[cg_num];
    C = new double*[cg_num];
    dw = new double **[cg_num];
    dw_sum = new double*[fg_num];

    for (int i = 0; i < cg_num; i++)
    {
        W[i] = new double[fg_num];
        C[i] = new double[fg_num];
        dw[i] = new double*[fg_num];

        for (int j = 0; j < fg_num; j++)
        {
            dw[i][j] = new double[3];
            W[i][j]  = 0.0;
            C[i][j]  = 0.0;
            for (int idim = 0; idim < 3; idim++) dw[i][j][idim] = 0.0;
        }
    }

    for (int i = 0; i < fg_num; i++)
    {
        dw_sum[i] = new double[3];
        for (int idim = 0; idim < 3; idim++) dw_sum[i][idim] = 0.0;
    }
    w_sum = new double[fg_num];
    memset(w_sum, 0, sizeof(double)*fg_num);

    L    = fg_atoms->L;
    rcut = cg_sites->rcut;

    compute();
}

//cleanup
void Matrix_C::cleanup()
{
    delete []w_sum;
    for (int i = 0; i < cg_num; i++)
    {
        for (int j = 0; j < fg_num; j++) delete[] dw[i][j];
        delete[] dw[i];
        delete[] W[i];
        delete[] C[i];
    }
    for (int i = 0; i < fg_num; i++) delete[] dw_sum[i];
    delete[] dw_sum;
    delete[] W;
    delete[] C;
    delete[] dw;
}

// Generating the matrices
void Matrix_C::matrixGenerator()
{
    // Generating the W matrix
    double **R = cg_sites->R;
    double **r = fg_atoms->r;
    for (int i = 0; i < cg_num; i++)
        for (int j = 0; j < fg_num; j++)
        {
            W[i][j] = w_Ij(R[i], r[j]);
            weight_deriv(R[i], r[j], dw[i][j]);
        }
    // Generating the w_sum matrix
    for (int i = 0; i < fg_num; i++)
    {
        w_sum[i] = 0.0;
        dw_sum[i][0] = 0.0;
        dw_sum[i][1] = 0.0;
        dw_sum[i][2] = 0.0;
        for (int j = 0; j < cg_num; j++)
        {
            w_sum[i] += W[j][i];
            dw_sum[i][0] += dw[j][i][0];
            dw_sum[i][1] += dw[j][i][1];
            dw_sum[i][2] += dw[j][i][2];
        }
    }

    // Generating the C matrix
    for (int i = 0; i < cg_num; i++)
    {
        double norm = 0.0;
        for (int j = 0; j < fg_num; j++)
        {
            C[i][j] = W[i][j] / w_sum[j];
            norm   += C[i][j];
        }

        for (int j = 0; j < fg_num; j++)
        {
            C[i][j] /= norm;
        }
    }
}

void Matrix_C::compute()
{
    matrixGenerator();
}
