#include "mapping.h"

int main(int narg, char** args)
{
   Mapping* test = new Mapping(narg, args);
   test->exec();
   delete test;
   return 0;
}
