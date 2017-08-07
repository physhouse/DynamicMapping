#ifndef _CG_SITES_H
#define _CG_SITES_H

#include <fstream>
#include <cstring>
#include "pointers.h"

class Cg_sites : protected Pointers
{
public:
    Cg_sites(class Mapping *);
    ~Cg_sites();
    void init(int argc, char **argv);
    void firstMapping();
    void output();
    void cleanup();

    double   L;
    double   rcut;
    double       sigmaOverRcut;
    int      cg_num;
    double   timestep;
    double   tol;  //tolerance for the initial mapping
    int      niter;  //Rounds of Iteration
    int      freq; // frequency of re-mapping

    double **R;
    double  *V;
    double  *VMAP;
    double  *M;
    std::ofstream cgtrj;
};

#endif
