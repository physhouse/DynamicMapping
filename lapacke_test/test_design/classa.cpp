#include <cstdio>
#include <cstdlib>
#include "classa.h"
#include "classb.h"

A::A(class Machine *m) : Base(m) {}

void A::init()
{
    printf("A was initialized\n");
    entry = 10;
}

void A::exec()
{
    printf("From A B.entry = %d\n", mb->entry);
}
