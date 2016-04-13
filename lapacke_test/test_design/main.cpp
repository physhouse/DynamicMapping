#include "machine.h"

int main()
{
    Machine *m = new Machine();
    m->init();
    m->exec();
    return 0;
}
