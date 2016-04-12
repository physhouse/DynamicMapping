#include "cg_sites.h"
#include "fg_atoms.h"
#include "engine.h"
#include "mapping.h"
#include "matrix_C.h"
#include "matrix_M.h"
#include "matrix_N.h"
#include <cstring>

Mapping::Mapping(int narg, char **arg)
{
  cg_sites = new Cg_sites(this);
  fg_atoms = new Fg_atoms(this);
  engine = new Engine(this);
  matrix_C = new Matrix_C(this);
  matrix_M = new Matrix_M(this);
  matrix_N = new Matrix_N(this);

  init(narg, arg);
}

Mapping::~Mapping()
{
  destroy();
}

void Mapping::init(int narg, char **arg)
{
  fg_atoms->init(narg, arg);
  cg_sites->init();
  engine->init();
  matrix_C->init();
  matrix_M->init();
  matrix_N->init();
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
  delete matrix_N;
  delete matrix_C;
}
