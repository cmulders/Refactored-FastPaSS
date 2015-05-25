cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
find_package(CUDA REQUIRED)

set(CPP_SOURCES
	Fastpass_main.cpp
	Fastpass_main.h	
 	fastPassLibrary.cpp
	fastPassLibrary.h
	fastPassOutput.cpp
	fastPassOutput.h
	fastPassParams.cpp
	fastPassParams.h
	fastPassQuery.cpp
	fastPassQuery.h
)

set(CUDA_SOURCES
	gpuDevice.cu
	gpuDevice.cuh
)

CUDA_ADD_EXECUTABLE(FastPaSS ${CPP_SOURCES} ${CUDA_SOURCES})

cuda_add_cublas_to_target(FastPaSS)
target_link_libraries(FastPaSS cuda)
target_link_libraries(FastPaSS PugiXML StringEncoders)
install(TARGETS FastPaSS DESTINATION ${PROJECT_SOURCE_DIR}/bin)
