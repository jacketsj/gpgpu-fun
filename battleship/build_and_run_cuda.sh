/opt/hipSYCL/bin/syclcc-clang precomputed_bitmask_unroll_custom_gpu_parallelized_mem_opt.cpp --hipsycl-targets=cuda:sm_75 --cuda-path=/usr/local/cuda -O3 -o precomputed_bitmask_unroll_custom_gpu_parallelized_mem_opt_cuda
time ./precomputed_bitmask_unroll_custom_gpu_parallelized_mem_opt_cuda
# non-gpu version
#g++ -O3 precomputed_bitmask_unroll_custom.cpp -o run_precomputed_bitmask_unroll_custom
#time ./run_precomputed_bitmask_unroll_custom
