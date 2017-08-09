#ifndef _VECTOR_CNV_H
#define _VECTOR_CNV_H

#include "pointers.h"

class Vector_CNv : protected Pointers
{
public:
    Vector_CNv(class Mapping *);
    ~Vector_CNv();
    void init();
    void generate_CNv();
    void compute();
    void cleanup();

    double   **CNv; // (C + N)v vector, dim * cg
    int      cg_num;
    int      fg_num;
};

#endif