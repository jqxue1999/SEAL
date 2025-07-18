## Library Setup
```bash
bash setup.sh
```

## Compile
```bash
cd mytest
mkdir build
cd build
cmake ..
make
```
## Run the Benchmarks
```bash
# Test Baseline
./bin/baseline_bench

# Test Coefficients with BLAS
./bin/coefficients_blas_bench

# Test Coefficients with SEAL Function
# (negacyclic_multiply_poly_mono_coeffmod)
./bin/coefficients_seal_bench

# Test Digits
./bin/digits_bench

# Test Plaintext with BLAS
./bin/plain_blas_bench

# Test Plaintext with naive OMP
./bin/plain_omp_bench
```