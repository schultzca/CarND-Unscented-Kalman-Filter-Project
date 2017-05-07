#!usr/bin/sh

# create build directory
mkdir -p build

# move into build directory
cd build

# compile project
cmake .. && make

# move to base directory
cd ..

# Create output directory
mkdir -p output

# run against third dataset
./build/UnscentedKF ./data/obj_pose-laser-radar-synthetic-input.txt ./output/output.txt