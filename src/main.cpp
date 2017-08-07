#include "mapping.h"

int main(int narg, char **args)
{
    Mapping *mapping_obj = new Mapping(narg, args);
    mapping_obj->exec();
    delete mapping_obj;
    return 0;
}
