#include "classa.h"
#include "classb.h"
#include <cstdlib>
#include <cstdio>

Machine::Machine()
{
  ma = new A(this);
  mb = new B(this);
}

void Machine::init()
{
  ma->init();
  mb->init();
}

void Machine::exec()
{
  ma->exec();
  mb->exec();
}
