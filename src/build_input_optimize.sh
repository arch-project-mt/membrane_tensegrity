mkdir -p build
g++ geometry.h rmsd.h data_io.hpp data_io.cpp optimize_structure_from_input.cpp -O3 -o ./build/io_optimization -std=c++14 --verbose
echo "[Build finished!]"
#./build/io_optimization
