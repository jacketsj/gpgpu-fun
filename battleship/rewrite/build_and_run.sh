/opt/hipSYCL/bin/syclcc-clang --hipsycl-platform=cpu  main_gpu.cpp -O2 --std=c++17 -o main_cpu
echo "Compilation done, running HipSycl (CPU)..."
time ./main_cpu
g++ -O3 -g main.cpp -std=c++17 -o main
echo "Compilation done, running basic (CPU)..."
time ./main
# non-gpu version
#g++ -O3 precomputed_bitmask_unroll_custom.cpp -o run_precomputed_bitmask_unroll_custom
#time ./run_precomputed_bitmask_unroll_custom
