#ifndef b
#define b

#include "base.h"

class B : public Base
{
public:
    B(class Machine *m);
    void init();
    void exec();
    int entry;
};

#endif
