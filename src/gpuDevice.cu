/*
 * cudaInitialize.cu
 *
 *  Created on: Jan 8, 2015
 *      Author: geforce
 */

#include "gpuDevice.cuh"
#include <ostream>
#include <cuda_occupancy.h>

#define BLOCK_SIZE 32
#define TILES (MZ_RANGE/BLOCK_SIZE-1)
__global__ void MatrixDotDotBias(
		float * deviceQueryMatrix, float * deviceLibraryMatrix, float * deviceDotMatrix, float * deviceBiasMatrix)
{

	__shared__ float A_tile[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float B_tile[BLOCK_SIZE][BLOCK_SIZE];

	float dotSub = 0, dotBias = 0, tempValue;

	// Row i of matrix deviceDotMatrix
	int iRow = blockIdx.y * blockDim.y + threadIdx.y;
	// Column j of matrix deviceDotMatrix
	int jCol = blockIdx.x * blockDim.x + threadIdx.x;

	// Accumulate deviceDotMatrix tile by tile.
#pragma unroll
	for (int tileIdx = 0; tileIdx < TILES; tileIdx++)
	{
		// Load one tile of A and one tile of B into shared mem

		int queryIndex = iRow * MZ_RANGE + tileIdx * blockDim.x + threadIdx.x;
		int libIndex = jCol * MZ_RANGE + tileIdx * blockDim.y + threadIdx.y;

		// Load A(i,j) to shared mem
		A_tile[threadIdx.y][threadIdx.x] = deviceQueryMatrix[queryIndex];
		// Load B(i,j) to shared mem
		B_tile[threadIdx.y][threadIdx.x] = deviceLibraryMatrix[libIndex];
		// Synchronize before computation
		__syncthreads();

		// Accumulate one tile of C from tiles of A and B in shared mem

		for (int k = 0; k < blockDim.x; k++)
		{
			tempValue = A_tile[threadIdx.y][k] * B_tile[k][threadIdx.x];
			// Accumulate for matrix C    // No Bank conflict
			dotSub += tempValue;
			dotBias += tempValue * tempValue;
		}
		// Synchronize
		__syncthreads();
	}

	// Store accumulated value to deviceDotMatrix(i,j) and deviceBiasMatrix(i,j)
	int outputIndex = iRow + jCol * (blockDim.y * gridDim.y);

	if ((dotBias < 0.0001) || (dotSub < 0.01))
		dotBias = 0.0;
	else
		dotBias = (sqrt(dotBias) / dotSub);

	deviceDotMatrix[outputIndex] = dotSub;
	deviceBiasMatrix[outputIndex] = dotBias;

}

namespace GPU {
/**
 * initialize cuda hardware
 */
Device::Device() :
		allocatedMemory(false), deviceLibraryMatrix(0), deviceQueryMatrix(0), deviceDotMatrix(0)
{

	int iDeviceCount;
	cudaDeviceProp *currentDeviceProp = new cudaDeviceProp();
	this->deviceProp = new cudaDeviceProp();

	double dUseDeviceScore = 0.0; // device score = clockRate * multiProcessorCount

	cuInit(0);
	cudaGetDeviceCount(&iDeviceCount);
	for (int curDevice = 0; curDevice < iDeviceCount; curDevice++)
	{
		memset(currentDeviceProp, 0, sizeof(*currentDeviceProp));
		if (cudaSuccess == cudaGetDeviceProperties(currentDeviceProp, curDevice))
		{
			if (currentDeviceProp->multiProcessorCount * currentDeviceProp->clockRate > dUseDeviceScore)
			{
				dUseDeviceScore = currentDeviceProp->multiProcessorCount * currentDeviceProp->clockRate;
				this->deviceIndex = curDevice;
			}
		}
		else
		{
			std::cerr << "error cannot read device properties for " << curDevice << std::endl;
		}
	}

	delete currentDeviceProp;

	memset(this->deviceProp, 0, sizeof(*this->deviceProp));
	if (cudaSuccess != cudaGetDeviceProperties(this->deviceProp, this->deviceIndex))
	{
		std::cerr << "error cannot read device properties for " << this->deviceIndex << std::endl;
		exit(EXIT_FAILURE);
	}

	gpuErrchk(cudaSetDevice(this->deviceIndex));

	gpuErrchkCublas(cublasCreate(&this->cublasHandle));
}

std::ostream& operator<<(std::ostream& out, const Device* gpuDevice)
{
	int version;
	cublasGetVersion(gpuDevice->cublasHandle, &version);
	out << "Using device " << gpuDevice->deviceIndex << ": '";
	out << gpuDevice->deviceProp->name << "', ";
	out << gpuDevice->deviceProp->totalGlobalMem / 1024 / 1024 << " MB, ";
	out << gpuDevice->deviceProp->multiProcessorCount << " multiprocessors, ";
	out << gpuDevice->deviceProp->clockRate / 1000 << " MHz" << std::endl;
	out << "CUBLAS version: " << version;
	return out;
}

Device::~Device()
{
	cublasDestroy(this->cublasHandle);
	this->freeMatrix();
}

size_t Device::getFreeMemory()
{
	size_t freeMemory, totalMemory;
	gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));

	return freeMemory;
}

bool Device::setMatrixSizes(int numberOfQueries, int numberOfLibraryEntries)
{
	//If there is no memory to be freed, exit
	if (this->allocatedMemory)
	{
		this->freeMatrix();
#ifdef DEBUG
		std::cout << "Auto freeing memory!" << std::endl;
#endif
	}
	//Keep track of real matrix dimensions
	this->numberOfQueries = numberOfQueries;
	this->numberOfLibraryEntries = numberOfLibraryEntries;

	//Calculate the padding
	this->paddedQueryRange = this->numberOfQueries + (BLOCK_SIZE - this->numberOfQueries % BLOCK_SIZE);
	this->paddedLibraryRange = this->numberOfLibraryEntries + (BLOCK_SIZE - this->numberOfLibraryEntries % BLOCK_SIZE);

	//Calculate the size of the float matrix
	this->sizeQueries = paddedQueryRange * MZ_RANGE * sizeof(float);
	this->sizeLibrary = paddedLibraryRange * MZ_RANGE * sizeof(float);
	this->sizeResults = this->numberOfQueries * this->numberOfLibraryEntries * sizeof(float);

	int deviceSizeResults = paddedQueryRange * paddedLibraryRange * sizeof(float);

#ifdef DEBUG_
	std::cout << "Query size:" << sizeQueries << std::endl;
#endif

	//Allocate the memory on the host
	gpuErrchk(cudaMallocHost((void ** ) &this->hostQueryMatrix, this->sizeQueries));
	gpuErrchk(cudaMallocHost((void ** ) &this->hostLibraryMatrix, this->sizeLibrary));
	gpuErrchk(cudaMallocHost((void ** ) &this->hostDotMatrix, this->sizeResults)); //not padded
	gpuErrchk(cudaMallocHost((void ** ) &this->hostBiasMatrix, this->sizeResults)); //not padded

	//Allocate the memory on the device
	gpuErrchk(cudaMalloc((void ** ) &this->deviceQueryMatrix, this->sizeQueries));
	gpuErrchk(cudaMalloc((void ** ) &this->deviceLibraryMatrix, this->sizeLibrary));
	gpuErrchk(cudaMalloc((void ** ) &this->deviceDotMatrix, deviceSizeResults)); //padded
	gpuErrchk(cudaMalloc((void ** ) &this->deviceBiasMatrix, deviceSizeResults)); //padded

	//Reset result and bias matrices
	cudaMemset(this->deviceDotMatrix, 0, deviceSizeResults);
	cudaMemset(this->deviceBiasMatrix, 0, deviceSizeResults);

	this->allocatedMemory = true;

	return true;
}

bool Device::freeMatrix()
{
	//If there is no memory to be freed, exit
	if (!this->allocatedMemory) return false;

	//Free the memory on the device
	gpuErrchk(cudaFree(this->deviceQueryMatrix));
	gpuErrchk(cudaFree(this->deviceLibraryMatrix));
	gpuErrchk(cudaFree(this->deviceDotMatrix));
	gpuErrchk(cudaFree(this->deviceBiasMatrix));

	//Free the memory on the host
	gpuErrchk(cudaFreeHost(this->hostQueryMatrix));
	gpuErrchk(cudaFreeHost(this->hostLibraryMatrix));
	gpuErrchk(cudaFreeHost(this->hostDotMatrix));
	gpuErrchk(cudaFreeHost(this->hostBiasMatrix));

	this->allocatedMemory = false;
	return true;
}

/*
 * Adds a query to the gpu memory and calculates the
 */
bool Device::addQueriesToMemory(float* binnedSpectra, int index)
{
	//If there is no memory set yet
	if (!this->allocatedMemory)
	{
		std::cerr << "No memory set yet with Device::setMatrixSizes" << std::endl;
		return false;
	}

	size_t offestMemory = index * MZ_RANGE;
	if (index > this->paddedQueryRange)
	{
		std::cerr << "Index out of bounds" << index << " = " << offestMemory << " Max index:" << this->paddedQueryRange
				<< std::endl;
		return false;
	}

	gpuErrchkCublas(
			cublasSetVector(MZ_RANGE, sizeof(binnedSpectra[0]), binnedSpectra, 1, (this->hostQueryMatrix + offestMemory), 1));

	return true;
}

bool Device::performDotProduct()
{
	if ((this->numberOfQueries*this->numberOfLibraryEntries) <= 10000 && false)
	{
		//Perform on CPU
		return matrixOnCPU();
	}
	else
	{
		//If there is no memory set yet
		if (!this->allocatedMemory)
		{
			std::cerr << "No memory set yet with Device::setMatrixSizes" << std::endl;
			return false;
		}

		//Copy over Host memory to device memory
		gpuErrchkCublas(
				cublasSetVector(this->paddedQueryRange * MZ_RANGE, sizeof(this->hostQueryMatrix[0]), this->hostQueryMatrix, 1, this->deviceQueryMatrix, 1));
		gpuErrchkCublas(
				cublasSetVector(this->paddedLibraryRange * MZ_RANGE, sizeof(this->hostLibraryMatrix[0]), this->hostLibraryMatrix, 1, this->deviceLibraryMatrix, 1));

		//Wait to finish kernel (possible query Euclidean norm calculation ongoing)
		gpuErrchk(cudaDeviceSynchronize());

		dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE, 1);
		dim3 dimGridResults(this->paddedLibraryRange / dimBlock.x, this->paddedQueryRange / dimBlock.y, 1);

#ifdef DEBUG_
		std::cout << "Grid (x:" << dimGridResults.x << ",y:" << dimGridResults.y << ") \t";
		std::cout << "Block (x:" << dimBlock.x << ",y:" << dimBlock.y << ")" << std::endl;
#endif

		MatrixDotDotBias<<<dimGridResults, dimBlock>>>(this->deviceQueryMatrix, this->deviceLibraryMatrix,
				this->deviceDotMatrix, this->deviceBiasMatrix);

		gpuErrchk(cudaDeviceSynchronize());

		//Copy over Device memory to host memory
		gpuErrchkCublas(
				cublasGetMatrix(this->numberOfQueries, this->numberOfLibraryEntries, sizeof(this->hostDotMatrix[0]),
						this->deviceDotMatrix, this->paddedQueryRange, this->hostDotMatrix, this->numberOfQueries));
		gpuErrchkCublas(
				cublasGetMatrix(this->numberOfQueries, this->numberOfLibraryEntries, sizeof(this->hostBiasMatrix[0]),
						this->deviceBiasMatrix, this->paddedQueryRange, this->hostBiasMatrix, this->numberOfQueries));

		return true;
	}
}

bool Device::matrixOnCPU()
{

	int outIndex;
	float dot = 0;
	float bias = 0;
	float temp = 0;

	for(int iQueryIndex=0;iQueryIndex<this->numberOfQueries;iQueryIndex++)
	{
		for(int iLibIndex=0;iLibIndex<this->numberOfLibraryEntries;iLibIndex++)
		{
			dot = 0;
			bias = 0;

			for(int iBinIndex = 0; iBinIndex < MZ_RANGE; iBinIndex++)
			{
				temp = this->hostQueryMatrix[iBinIndex+iQueryIndex*MZ_RANGE] * this->hostLibraryMatrix[iBinIndex+iLibIndex*MZ_RANGE];
				dot += temp;
				bias += temp*temp;
			}

			if ((bias < 0.0001) || (dot < 0.01))
				bias = 0.0;
			else
				bias = (sqrt(bias) / dot);

			outIndex = iLibIndex+iQueryIndex*this->numberOfLibraryEntries;

			this->hostDotMatrix[outIndex] = dot;
			this->hostBiasMatrix[outIndex] = bias;
		}
	}
	return true;
}

} //namespace GPU
