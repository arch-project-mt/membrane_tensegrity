mkdir -p build && cd build
cmake ..
make
if [ $? -ne 0 ]; then
    echo "Build failed"
    exit 1
fi
./membrane_tensegrity
