#/opt/hipSYCL/bin/syclcc-clang main_gpu.cpp --hipsycl-targets=cuda:sm_75 --cuda-path=/usr/local/cuda -O2 -o main_gpu -funroll-loops
/opt/hipSYCL/bin/syclcc-clang main_gpu.cpp --hipsycl-targets=cuda:sm_75 --cuda-path=/usr/lib/cuda -O2 -o main_gpu -funroll-loops
echo "Compilation done, running..."
time ./main_gpu
# non-gpu version
#g++ -O3 precomputed_bitmask_unroll_custom.cpp -o run_precomputed_bitmask_unroll_custom
#time ./run_precomputed_bitmask_unroll_custom
