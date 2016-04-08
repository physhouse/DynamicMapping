#include "mapping.h"
#include "engine.h"
#include "fg_atoms.h"
#include "cg_sites.h"
#include "matrix_C.h"
#include "matrix_M.h"
#include "matrix_N.h"

int main(int narg, char** args)
{
  Mapping *mapping = new Mapping(narg, args);
  mapping->exec();
  delete mapping;
  return 0;
}
