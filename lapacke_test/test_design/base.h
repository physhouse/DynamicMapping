#ifndef base
#define base
#include "machine.h"

class Base
{
public:
  Base(Machine *ptr) :
    m(ptr),
    ma(ptr->ma),
    mb(ptr->mb) {}

protected:
   Machine *m;
   A	   *&ma;
   B	   *&mb;
};

#endif
