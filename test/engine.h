#ifndef ENGINE_H
#define ENGINE_H

#include "pointers.h"

class Engine : protected Pointers
{
public:
  Engine(class Mapping *);
  ~Engine();
  void init();

  void update();
  void exec();

  void initFrame();
  void matrixSolver();
  void integrate();
  void endOfFrame();
  void cleanup();

private:
  double** MCG;
  double** MFG;

  int     cg_num;
  int	  fg_num;
  int	  nsteps;
};

#endif
