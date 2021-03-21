#/opt/hipSYCL/bin/syclcc-clang --hipsycl-platform=rocm --hipsycl-gpu-arch=gfx802 precomputed_bitmask_unroll_custom_gpu.cpp -O2 --std=c++17 -o run_precomputed_bitmask_unroll_custom_gpu
/opt/hipSYCL/bin/syclcc-clang --hipsycl-platform=rocm --hipsycl-gpu-arch=gfx802 precomputed_bitmask_unroll_custom_gpu_parallelized.cpp -O2 --std=c++17 -o run_precomputed_bitmask_unroll_custom_gpu_parallelized
#/opt/hipSYCL/bin/syclcc-clang --hipsycl-platform=cpu  precomputed_bitmask_unroll_custom_gpu_parallelized.cpp -O2 --std=c++17 -o run_precomputed_bitmask_unroll_custom_gpu_parallelized
