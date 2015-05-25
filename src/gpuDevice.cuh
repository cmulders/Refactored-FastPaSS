/*
 * cudaInitialize.h
 *
 *  Created on: Jan 8, 2015
 *      Author: geforce
 */

#ifndef CUDAINITIALIZE_H_
#define CUDAINITIALIZE_H_

#define MZ_START 10
#define MZ_END 2010
#define MZ_RANGE (MZ_END - MZ_START + 1)

#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <iosfwd>
#include <iostream>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline bool gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		std::cerr << "GPUassert:" << cudaGetErrorString(code) << " in " << file << ":" << line << std::endl;
		if (abort)
		{
			exit(code);
		}
		return false;
	}
	return true;
}

#define gpuErrchkCublas(ans) { gpuAssertCublas((ans), __FILE__, __LINE__); }
inline bool gpuAssertCublas(cublasStatus_t code, const char *file, int line, bool abort = true)
{
	if (code != CUBLAS_STATUS_SUCCESS)
	{
		std::cerr << "GPUassert CUBLAS:" << code << " in " << file << ":" << line << std::endl;
		if (abort)
		{
			exit(code);
		}
		return false;
	}
	return true;
}

namespace GPU {
class Device
{
		//Flag if there is currently memroy allocated
		bool allocatedMemory;

		//Keep track of the real matrix sizes
		int numberOfQueries, numberOfLibraryEntries;

		//Device memory pointers
		float *deviceQueryMatrix, *deviceLibraryMatrix, *deviceDotMatrix, *deviceBiasMatrix;

		bool matrixOnCPU();

	public:
		size_t sizeQueries, sizeLibrary, sizeResults;

		int paddedQueryRange, paddedLibraryRange;

		//Host memory pointers
		float *hostQueryMatrix, *hostLibraryMatrix, *hostDotMatrix, *hostBiasMatrix;

		int deviceIndex;
		cudaDeviceProp *deviceProp;
		cublasHandle_t cublasHandle;
		Device();
		~Device();

		size_t getFreeMemory();

		bool setMatrixSizes(int numberOfQueries, int libraryRange);

		bool freeMatrix();

		bool addQueriesToMemory(float* binnedSpectra, int index);
		bool performDotProduct();

		friend std::ostream& operator<<(std::ostream& os, const Device* gpuDevice);

};

} // namespace GPU

#endif /* CUDAINITIALIZE_H_ */
