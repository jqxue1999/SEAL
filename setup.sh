mkdir seallibs
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=./seallibs
cmake --build build
sudo cmake --install build

git clone https://github.com/google/benchmark.git
cd benchmark
git clone https://github.com/google/googletest.git
cmake -E make_directory build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/jiaq/Research/SEAL/benchmarklibs -DBENCHMARK_DOWNLOAD_DEPENDENCIES=ON
cmake --build build --config Release
cmake --install build --config Release