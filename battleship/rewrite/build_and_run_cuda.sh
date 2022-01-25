# old version
# /opt/hipSYCL/bin/syclcc-clang main_gpu.cpp --hipsycl-targets=cuda:sm_75 --cuda-path=/usr/lib/cuda -O2 -o main_gpu -funroll-loops
# echo "Compilation done, running..."
# time ./main_gpu
#
#/opt/hipSYCL/bin/syclcc-clang main_gpu.cpp --hipsycl-targets=cuda:sm_75 --cuda-path=/usr/local/cuda -O2 -o main_gpu -funroll-loops
# non-gpu version
#g++ -O3 precomputed_bitmask_unroll_custom.cpp -o run_precomputed_bitmask_unroll_custom
#time ./run_precomputed_bitmask_unroll_custom
# api version
/opt/hipSYCL/bin/syclcc-clang main_gpu_api.cpp --hipsycl-targets=cuda:sm_75 --cuda-path=/usr/lib/cuda -O2 -o main_gpu_api -funroll-loops -g
echo "GPU compilation done, compiling front-end..."
g++ -O3 api_user.cpp -o api_user
echo "Compilation done, running..."
./api_user
