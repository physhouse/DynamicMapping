#ifndef _MATRIX_M_H
#define _MATRIX_M_H

#include "pointers.h"

class Matrix_M : protected Pointers
{
public:
    Matrix_M(class Mapping *);
    ~Matrix_M();
    void init();
    void generate_M();
    void compute();
    void cleanup();
    double distance(double *, double *);

    double **   **M; // M matrix, dim * dim * cg * cg
    int      cg_num;
    int      fg_num;
};

#endif
