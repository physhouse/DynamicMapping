#include "cg_sites.h"
#include "fg_atoms.h"
#include "engine.h"
#include "mapping.h"
#include "matrix_C.h"
#include "matrix_M.h"
#include "vector_CNv.h"
#include "neighbor.h"
#include <cstring>

Mapping::Mapping(int narg, char **arg)
{
    cg_sites = new Cg_sites(this);
    fg_atoms = new Fg_atoms(this);
    engine = new Engine(this);
    matrix_C = new Matrix_C(this);
    matrix_M = new Matrix_M(this);
    vector_CNv = new Vector_CNv(this);
    neighbor = new Neighbor(this);

    init(narg, arg);
}

Mapping::~Mapping()
{
    destroy();
}

void Mapping::init(int narg, char **arg)
{
    fg_atoms->init(narg, arg);
    cg_sites->init(narg, arg);
    engine->init(narg);
    neighbor->init();
    matrix_C->init();
    matrix_M->init();
    vector_CNv->init();
}

void Mapping::exec()
{
    engine->exec();
}

void Mapping::destroy()
{
    delete engine;
    delete fg_atoms;
    delete cg_sites;
    delete matrix_M;
    delete vector_CNv;
    delete matrix_C;
    delete neighbor;
}
