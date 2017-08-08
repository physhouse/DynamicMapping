#ifndef _MATRIX_N_H
#define _MATRIX_N_H

#include "pointers.h"

class Matrix_N : protected Pointers
{
public:
    Matrix_N(class Mapping *);
    ~Matrix_N();
    void init();
    void generate_N();
    void compute();
    void cleanup();

    double **   **N; // N matrix, dim * dim * cg * fg
    int      cg_num;
    int      fg_num;
};

#endif
