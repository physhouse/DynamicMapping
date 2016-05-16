#ifndef a
#define a
#include "base.h"

class A : public Base
{
public:
    A(class Machine *m);
    void  init();
    void  exec();
    int entry;
};

#endif
