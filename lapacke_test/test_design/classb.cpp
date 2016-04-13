#include <cstdio>
#include "classb.h"
#include "classa.h"
#include "base.h"
#include <cstdlib>

B::B(class Machine *m) : Base(m) {}
void B::init()
{
    printf("B was initialized\n");
    entry = 30;
}

void B::exec()
{
    printf("From B, A.entry = %d\n", ma->entry);
}
