# DynamicMapping

The code in src files should not be used yet, still being tested

build:

cd test/

make

./test.x trj.in 20 // arg[1]->filename of FG trajectory, arg[2]->Number of Frames in FG trajectories

This code seems to be numerically stable for the LJ system, still need to be verified.

TODO: Implement neighbor list to accelerate, now the algorithm scales as O(N^2*M^2)
