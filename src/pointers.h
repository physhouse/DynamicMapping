#ifndef _POINTERS_H
#define _POINTERS_H

#include "mapping.h"

class Pointers
{
public:
    Pointers(Mapping *ptr) :
        map(ptr),
        cg_sites(ptr->cg_sites),
        fg_atoms(ptr->fg_atoms),
        engine(ptr->engine),
        matrix_C(ptr->matrix_C),
        matrix_M(ptr->matrix_M),
        vector_CNv(ptr->vector_CNv),
        neighbor(ptr->neighbor) {}

protected:
    Mapping   *map;
    Cg_sites*&cg_sites;
    Fg_atoms*&fg_atoms;
    Engine   *&engine;
    Matrix_C*&matrix_C;
    Matrix_M*&matrix_M;
    Vector_CNv*&vector_CNv;
    Neighbor*&neighbor;

};

#endif
