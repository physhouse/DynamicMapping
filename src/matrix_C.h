#ifndef _MATRIX_H
#define _MATRIX_H

#include "pointers.h"

class Matrix_C : protected Pointers
{
public:
    Matrix_C(class Mapping *);
    ~Matrix_C();

    void init();
    double w_Ij(double *R, double *r);
    double calc_proximity(const double r) const;
    void calc_proximity_deriv(const double * const R, const double * const r, double * const dwvec) const;
    void calc_proximity_and_deriv(const double * const R, const double * const r, double &w, double * const dw_vec) const;

    void matrixGenerator();
    double generateWeights();
    void generateMass();
    void generateMatrixC();
    void iterateSolver();
    void cleanup();
    void compute();

    double **W;  // cg_num by fg_num weights Matrix W
    double  *w_sum;  // fg vector, W_sum[j] = sum over J W[J][j]
    double **    *dw; // dw/dx matrix, dimension fg by cg by 3
    double **dw_sum; // dw/dx sum matrix, dimension fg by 3
    double **C;  // C matrix, dimension cg*fg
    double      *sumC;   // sumC[j] = sum over i, C[i][j]
    double      *F;   // F[I]
    int      cg_num;
    int      fg_num;
    double   L;
    double   rcut;   // Cutoff distance for the tanh(r-rc/sigma) function
    double   sigma;  // value of sigma of the tanh(r-rc/sigma) function

    double   width;  // Width of Gaussian Distribution
    double   M0; // Constraint Mass

};

#endif
