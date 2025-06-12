mkdir seallibs
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=./seallibs
cmake --build build
sudo cmake --install build