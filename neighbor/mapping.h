#ifndef _MAPPING_H
#define _MAPPING_H

#include <cstdio>

class Mapping
{
public:
  class Fg_atoms *fg_atoms;
  class Cg_sites *cg_sites;
  class Engine	 *engine;
  class Matrix_C *matrix_C;
  class Matrix_M *matrix_M;
  class Matrix_N *matrix_N;
  class Neighbor *neighbor;

  Mapping(int, char **);
  ~Mapping();

  void init(int, char **);
  void exec();
  void destroy();
};

#endif
