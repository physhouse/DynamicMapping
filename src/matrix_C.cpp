#include "mapping.h"
#include "pointers.h"
#include "matrix_C.h"
#include "fg_atoms.h"
#include "cg_sites.h"
#include "geom.h"
#include "neighbor.h"
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
    double sq_distance = calc_sq_distance(R, r, L);
    return weight(sqrt(sq_distance));
}

double Matrix_C::weight(double r)
{
    /* Gaussian */
    //double r0 = r / rcut;
    //return exp(- r0 * r0);

    /* tanh */
    //return 0.5 - 0.5 * tanh((r - rcut) / (rcut * 0.66));
    double sigma = rcut * sigmaOverRcut;
    return 0.5 - 0.5 * tanh((r - rcut) / sigma);
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

    distance = sqrt(distance);

    //Gaussian
    //double r0 = distance / rcut;
    //double dwdr = - 2.0 * (r0 / rcut) * exp(-r0 * r0);

    //double tanhr = tanh((distance - rcut) / (rcut * 0.66));
    //double dwdr = -0.5 * (1.0 - tanhr * tanhr) / (rcut * 0.66);

    double sigma = rcut * sigmaOverRcut;
    double tanhr = tanh((distance - rcut) / sigma);
    double dwdr = -0.5 * (1.0 - tanhr * tanhr) / sigma;

    for (int idim = 0; idim < 3; idim++)
    {
        double r_dim = R[idim] - r[idim];
        wrap_coord_relative(r_dim, 0, L);
        dw_vec[idim] = dwdr * r_dim / distance;
    }
}

// Initialization
void Matrix_C::init()
{
    fg_num = fg_atoms->fg_num;
    cg_num = cg_sites->cg_num;

    W = new double*[cg_num];
    C = new double*[cg_num];
    F = new double[cg_num];
    sumC = new double[fg_num];
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
    memset(w_sum, 0, sizeof(double) * fg_num);
    memset(sumC, 0, sizeof(double) * fg_num);

    L    = fg_atoms->L;
    rcut = cg_sites->rcut;
    sigmaOverRcut = cg_sites->sigmaOverRcut;

    // set up the gaussian distribution
    width = 4.8;

    // Initialize the mass
    M0 = (double)fg_num / (double)cg_num;
    for (int i = 0; i < cg_num; i++)
    {
        F[i] = M0;
    }

    compute();
}

//cleanup
void Matrix_C::cleanup()
{
    delete [] w_sum;
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
    delete[] F;
    delete[] dw;
}

// Generating the matrices
void Matrix_C::matrixGenerator()
{
    // Generating the W matrix
    double **R = cg_sites->R;
    double **r = fg_atoms->r;
    int     *numFgNeighbors = neighbor->numFgNeighbors;
    int    **fgList = neighbor->fgList;

    for (int i = 0; i < cg_num; i++)
    {
        for (int j = 0; j < neighbor->numNeighbors[i]; j++)
        {
            int index = neighbor->list[i][j];
            W[i][index] = w_Ij(R[i], r[index]);
            weight_deriv(R[i], r[index], dw[i][index]);
        }
    }

    // Generating the w_sum matrix
    for (int i = 0; i < fg_num; i++)
    {
        w_sum[i] = 0.0;
        dw_sum[i][0] = 0.0;
        dw_sum[i][1] = 0.0;
        dw_sum[i][2] = 0.0;
        for (int j = 0; j < numFgNeighbors[i]; j++)
        {
            int index = fgList[i][j];
            w_sum[i] += W[index][i];
            dw_sum[i][0] += dw[index][i][0];
            dw_sum[i][1] += dw[index][i][1];
            dw_sum[i][2] += dw[index][i][2];
        }
    }

    // Generating the C matrix
    for (int i = 0; i < cg_num; i++)
    {
        double norm = 0.0;
        for (int j = 0; j < neighbor->numNeighbors[i]; j++)
        {
            int index = neighbor->list[i][j];
            C[i][index] = W[i][index] / w_sum[index];
            norm   += C[i][index];
        }

        for (int j = 0; j < fg_num; j++)
        {
            C[i][j] /= norm;
        }
    }

    for (int i = 0; i < fg_num; i++)
    {
        sumC[i] = 0.0;
        for (int j = 0; j < numFgNeighbors[i]; j++)
        {
            int index = fgList[i][j];
            sumC[i] += C[index][i];
        }
    }

    generateMass();
}

double Matrix_C::massWeights(double M)
{
    return exp(- (M - M0) * (M - M0) / (2 * width * width));
}

double Matrix_C::generateWeights()
{
    double **R = cg_sites->R;
    double **r = fg_atoms->r;
    int     *numNeighbors = neighbor->numNeighbors;
    int     *numFgNeighbors = neighbor->numFgNeighbors;
    int    **list = neighbor->list;
    int    **fgList = neighbor->fgList;

    double error = 0.0;
    double alpha = 0.05;

    for (int I = 0; I < cg_num; I++)
    {
        for (int i = 0; i < numNeighbors[I]; i++)
        {
            int index = list[I][i];
            double temp = W[I][index];
            W[I][index] = (1.0 - alpha) * temp + alpha * weightsVarying(R[I], r[index], F[I]);
            error += (W[I][index] - temp) * (W[I][index] - temp);
        }
    }

    double test_sum = 0.0;
    for (int i = 0; i < fg_num; i++)
    {
        w_sum[i] = 0.0;
        for (int j = 0; j < numFgNeighbors[i]; j++)
        {
            int index = fgList[i][j];
            w_sum[i] += W[index][i];
        }
        test_sum += w_sum[i];
    }
    return sqrt(error);
}

void Matrix_C::generateMass()
{
    double *M = cg_sites->M;
    int      *numNeighbors = neighbor->numNeighbors;
    int    **list = neighbor->list;

    for (int I = 0; I < cg_num; I++)
    {
        F[I] = 0.0;
        M[I] = 0.0;
        for (int i = 0; i < numNeighbors[I]; i++)
        {
            int index = list[I][i];
            F[I] += W[I][index] / w_sum[index];
            //M[I] += C[I][index] / sumC[index];
        }
        M[I] = F[I];
    }

}

void Matrix_C::generateMatrixC()
{
    int     *numNeighbors = neighbor->numNeighbors;
    int     *numFgNeighbors = neighbor->numFgNeighbors;
    int    **list = neighbor->list;
    int    **fgList = neighbor->fgList;

    //Generating the C Matrix
    for (int i = 0; i < cg_num; i++)
    {
        double norm = 0.0;
        for (int j = 0; j < numNeighbors[i]; j++)
        {
            int index = list[i][j];
            C[i][index] = W[i][index] / w_sum[index];
            norm += C[i][index];
        }

        for (int j = 0; j < fg_num; j++)
        {
            C[i][j] /= norm;
        }
    }

    for (int i = 0; i < fg_num; i++)
    {
        sumC[i] = 0.0;
        for (int j = 0; j < numFgNeighbors[i]; j++)
        {
            int index = fgList[i][j];
            sumC[i] += C[index][i];
        }
    }
}

void Matrix_C::compute()
{
    //iterateSolver();
    matrixGenerator();
}
