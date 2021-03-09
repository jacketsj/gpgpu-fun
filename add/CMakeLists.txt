add_executable(add add.cpp)
add_sycl_to_target(TARGET add SOURCES add.cpp)
#install(TARGETS add RUNTIME DESTINATION share/hipSYCL/examples/)
