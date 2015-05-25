/*
 * Fastpass_main.cpp
 *
 *  Created on: Jan 6, 2015
 */

#include "Fastpass_main.h"
#include "fastPassParams.h"
#include "fastPassLibrary.h"
#include "fastPassQuery.h"
#include "fastPassOutput.h"
#include "gpuDevice.cuh"

#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ctime> //For timing
#ifdef DEBUG
#define __GLIBCXX_DEBUG
#define __GLIBCXX_FORCE_NEW
#endif

#define DEF_STARTTIMER  startTime = std::clock(); endTime = 0;

#define DEF_TIME(MES) endTime = std::clock(); totalTime += (endTime-startTime); \
				std::cout << MES << " - " \
				<< std::fixed << std::setprecision(3) << (endTime - startTime) / (double) (CLOCKS_PER_SEC) << " s " \
				<< std::fixed << std::setprecision(2) << "(Total time: " << totalTime/ (double) (CLOCKS_PER_SEC) << " s)" \
				<< std::endl << std::endl;

// display a matrix
//											  rows    columns
void display(const char* str, const float* a, int nr, int nc)
{
//	int temp = nr;
//	nr = nc;
//	nc = temp;

	std::cout << str << std::endl;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
			std::cout << std::setw(6) << std::fixed << std::setprecision(3) << a[j + i * nc]; // a[j + i * nc] * 100;
		std::cout << "\n";
	}
	std::cout << std::endl;
}

int main(int argc, char *argv[])
{
	//Timers for user information
	std::clock_t startTime, endTime; //Block timers
	std::clock_t totalTime = 0; //Total block time
	std::clock_t wallTime = std::clock(); //Total time for the while program

	//Batch progress bar
	size_t totalNumberOfBatches = 0;
	size_t currentBatch = 0;
#define PROGRESSBAR_MAX 20
	int lastProgress = 0;

	//The device we are going to use
	GPU::Device *gpuDevice;
	//The parameters for this programs
	Parameters::PStruct *fastpassParameters;
	//The library used for searching the query spectra
	FastPassLibrary::Info *fastpassLibraryInfo;
	//The current query we are processing
	FastPassQuery::Query *fastPassCurrentQuery;
	//Output each query into their respective file
	FastPassOutput::PepXML *fastPassOutput;

	//Introduction to FastPass on the commandline
	std::cout << "FastPaSS" << std::endl;
	std::cout << "Created by Lydia Ashleigh Baumgardner et al. in  ";
	std::cout << "'Fast parallel tandem mass spectral library searching using GPU hardware acceleration'" << std::endl;
	std::cout << "Adapted by Coen Mulders at the CPM department of the LUMC in Leiden (jan 2015)" << std::endl;
	std::cout << std::endl;
	std::cout << "mzXML parser depends on pugixml." << std::endl;
	std::cout << "This software is based on pugixml library (http://pugixml.org)." << std::endl;
	std::cout << "pugixml is Copyright (C) 2006-2014 Arseny Kapoulkine." << std::endl;
	std::cout << "--------------------------------------------------" << std::endl << std::endl;

	//Initialize parameter list with default values
	fastpassParameters = getParameters(argc, argv);

	//No parameters set, show help
	if (argc == 1)
	{
		showHelpParameters();
		return 0;
	}

	if (fastpassParameters->showParameters)
	{
		//Print default parameters for the user
		printParameters(fastpassParameters);

		//If this was the only argument, then exit
		if (argc == 2) return 0;
	}

	//Check if there are input files, otherwise show help
	if (argc > 1 && fastpassParameters->queryFiles.size() == 0)
	{
		fastpassParameters->flagShowHelp = true;
		std::cerr << "No input files were specified." << std::endl;
	}

	//Check if the user wants to convert a library ad an .splib is supplied
	if (fastpassParameters->convertSpLib && fastpassParameters->queryFiles.size() > 0
			&& fastpassParameters->queryFiles.at(0).substr(fastpassParameters->queryFiles.at(0).rfind(".")) == ".splib")
	{
		//Convert the splib to fpbin, with all the set parameters (so that the created library matches SpectraST as good as possible)
		FastPassLibrary::Converter *fastpassLibraryConverter = new FastPassLibrary::Converter(); //The converter for the splib
		fastpassLibraryConverter->setParameters(fastpassParameters);

		//Process each library
		for (auto & libFile : fastpassParameters->queryFiles)
		{
			fastpassLibraryConverter->convert(libFile);
		}
		delete fastpassLibraryConverter;
		return 0;
	}

	//Check if there is a library specified, if there is more than 1 argument and we are not converting a library
	if (argc > 1 && !fastpassParameters->libraryFile.empty()
			&& fastpassParameters->libraryFile.substr(fastpassParameters->libraryFile.rfind(".")) == ".splib")
	{
		//Library file is a .splib, try to find the .fpbin or convert the .splib
		std::cout << std::endl << "Splib supplied, check if " << Parameters::libraryExtension
				<< " lib exists in the same directory." << std::endl;

		std::string fplib(fastpassParameters->libraryFile);
		fplib.replace(fastpassParameters->libraryFile.rfind("."), std::string::npos, Parameters::libraryExtension);

		//We have an splib as library file, check if the corresponding fpbin exists
		std::ifstream f(fplib);
		bool fpLibExists = f.good();
		f.close();

		//Check if we do a forced convert
		if (!fpLibExists || fastpassParameters->convertSpLib)
		{
			//Notifiy the user about the converison that takes place
			if (fastpassParameters->convertSpLib)
				std::cout << "Performing forced conversion" << std::endl;
			else
				std::cout << "Fpbin did not exist, starting the conversion of the library" << std::endl;

			DEF_STARTTIMER
			{
				//Convert the splib to fpbin, with all the set parameters (so that the created library matches SpectraST as good as possible)
				FastPassLibrary::Converter *fastpassLibraryConverter = new FastPassLibrary::Converter(); //The converter for the splib
				fastpassLibraryConverter->setParameters(fastpassParameters);
				fastpassLibraryConverter->convert(fastpassParameters->libraryFile);
				delete fastpassLibraryConverter;
			}
			DEF_TIME("Library converted")
		}
		else
		{
			//fplib exists and we are going to use that one
			f.close();
			std::cout << "Fplib exists (" << fplib << "). No forced conversion required (add -C)" << std::endl;

		}
		//Assign our new library to the library parameter
		std::cout << "Replacing the libraryFile parameter '" << fastpassParameters->libraryFile << "' > '" << fplib
				<< "'" << std::endl << std::endl;
		fastpassParameters->libraryFile.assign(fplib);

	}

	//Check if there is a library specified, if there is more than 1 argument and we are not converting a library
	if (argc > 1 && fastpassParameters->libraryFile.empty())
	{
		fastpassParameters->flagShowHelp = true;
		std::cerr << "No library was specified." << std::endl;
	}

	/*
	 * Display help
	 * 	1. when asked with -H or flagShowHelp
	 * 	2. when an value was invalid.
	 * 	3. when there are no input files
	 *
	 * 	will exit the program
	 */
	if (fastpassParameters->flagShowHelp)
	{
		showHelpParameters();
		destroyParameters(fastpassParameters);
		return 0;
	}

	{ //Output the input files to the user, put into a block for convenience
		std::cout << "Input file" << (fastpassParameters->queryFiles.size() == 1 ? " is" : "s are") << ":";
		int index = 0;
		for (auto & queryFile : fastpassParameters->queryFiles)
		{
			//Add some nice spacing
			if ((index % 3) == 0)
				std::cout << std::endl << "	";
			else if (index > 0) std::cout << ", ";

			std::cout << queryFile;
			index++;
		}
		std::cout << std::endl << std::endl;
	}

	/*
	 * Done parameter checking
	 * 1. Loading the library
	 * 2. Initializing GPU
	 * 3. Running each query file
	 */

	//Load the library information
	DEF_STARTTIMER
	{
		fastpassLibraryInfo = new FastPassLibrary::Info();
		if (!fastpassLibraryInfo->readFPLibrary(fastpassParameters->libraryFile))
		{
			//Failed to open or read library: abort
			std::cerr << "Failed to open or read library." << std::endl;
			return 0;
		}
		fastpassLibraryInfo->printLibraryInfo();
	}
	DEF_TIME("Read library into memory")

	std::cout << "Initializing GPU device" << std::endl;
//Check if where on CUDA capable hardware and initialize it, getting the properties of this device
	gpuDevice = new GPU::Device();
	std::cout << gpuDevice << std::endl;
	std::cout << std::endl;

	uint64_t totalNumberOfSpectra = 0;
	/*
	 * Main loop to process each query file
	 */
	for (auto & queryFile : fastpassParameters->queryFiles)
	{
		fastPassCurrentQuery = new FastPassQuery::Query();
		fastPassCurrentQuery->setParameters(fastpassParameters);
		fastPassCurrentQuery->setLibrary(fastpassLibraryInfo);

		/*
		 * Open the Query file
		 */DEF_STARTTIMER
		{
			if (!fastPassCurrentQuery->load(queryFile))
			{
				std::cerr << "Could not open '" << queryFile << "', skipping this file" << std::endl;
				delete fastPassCurrentQuery;
				continue;
			}
			std::cout << "Opened '" << queryFile << "', processed " << fastPassCurrentQuery->queryPeptides.size()
					<< " queries" << std::endl;
		}
		DEF_TIME("File read CPU time")

		/*
		 * Filter the queries to the user specifications
		 */DEF_STARTTIMER
		{
			fastPassCurrentQuery->filterQueries();
		}
		DEF_TIME("Filtering Queries (Good spectra left:" << fastPassCurrentQuery->queryPeptides.size() << ")")

		if (fastPassCurrentQuery->queryPeptides.size() == 0)
		{
			std::cerr << "No queries left, skipping this file. " << std::endl;
		}
		totalNumberOfSpectra += fastPassCurrentQuery->queryPeptides.size();

		/*
		 * Bin the query spectra for multiplication
		 */DEF_STARTTIMER
		{
			fastPassCurrentQuery->binSpectraQueries();
		}
		DEF_TIME("Binning Queries")

		/*
		 * Assign the library ranges for each query
		 */DEF_STARTTIMER
		{
			fastPassCurrentQuery->assignLibraryrange();
		}
		DEF_TIME("Assigning library ranges to queries")

//		if (fastPassCurrentQuery->queryType != FastPassQuery::QueryType::mzxml)
//		{
//mzXML reader sorts already on mass, much faster
		/*
		 * Sort the queries on peptide mass
		 */DEF_STARTTIMER
		{
			fastPassCurrentQuery->sortOnPeptideMass();
		}
		DEF_TIME("Sorting on peptide mass")
//		}

		/*
		 * Create batches
		 */DEF_STARTTIMER
		{
			//Create the batches
			fastPassCurrentQuery->createBatches();

			std::cout << "Total number of batches: " << fastPassCurrentQuery->batches.size() << std::endl;
#ifdef DEBUG //Show amount of GPU memory available in debug mode
			size_t freeMemoryDevice = gpuDevice->getFreeMemory();
			std::cout << "GPU free memory:" << freeMemoryDevice / 1024 / 1024 << "MB" << std::endl;
#endif
		}
		DEF_TIME("Creating batches for GPU")
		totalNumberOfBatches = fastPassCurrentQuery->batches.size();
		currentBatch = 0;

		/*
		 * Process each batch on the GPU
		 */DEF_STARTTIMER
		std::cout << "Starting batch processing on the GPU" << std::endl;
//		int index = 0;
		for (auto & batch : fastPassCurrentQuery->batches)
		{
#ifdef DEBUG
			std::clock_t startBatch = std::clock();
			std::cout << "Batch: Queries: " << batch.numberOfQueries << " x Library entries: "
			<< batch.numberOfLibraryEntries << std::endl;
#endif
			//Initialize next dot product
			gpuDevice->setMatrixSizes(batch.numberOfQueries, batch.numberOfLibraryEntries);

			//read the library peaks into the matrix pointer
			fastpassLibraryInfo->readPeaksForLibraryRangeIntoPointer(batch.startQueryBatch->libraryMatchStart,
					batch.numberOfLibraryEntries, gpuDevice->hostLibraryMatrix);

			//Assign the queries to the Query memory
			for (auto currentQuery = batch.startQueryBatch; currentQuery != batch.endQueryBatch; currentQuery++)
			{
				gpuDevice->addQueriesToMemory(&(*currentQuery).binnedSpectra[0],
						std::distance(batch.startQueryBatch, currentQuery));
			}

#ifdef DEBUG
			double time = ((std::clock() - startBatch) / (double) (CLOCKS_PER_SEC));
			std::cout << "Memory copy time: " << std::fixed << std::setprecision(2) << time * 1000 << "ms"<< std::endl;

			startBatch = std::clock();
#endif

			//Copies over the matrices and performs dot product
			gpuDevice->performDotProduct();

//			std::cout << "---------------" << std::endl;
//			std::cout << batch.numberOfLibraryEntries << " x " << batch.numberOfQueries<< std::endl;
//			display("Host dot", gpuDevice->hostDotMatrix, batch.numberOfLibraryEntries, batch.numberOfQueries);

#ifdef DEBUG
			time = ((std::clock() - startBatch) / (double) (CLOCKS_PER_SEC));
			std::cout << std::fixed << std::setprecision(2) << time * 1000 << "ms -> " << std::setprecision(1)
			<< std::fixed << std::setw(7) << (time * 1000000000) / (gpuDevice->sizeResults / sizeof(float))
			<< "ns/cell (total:" << gpuDevice->sizeResults / sizeof(float) << ")" << std::endl;

			startBatch = std::clock();
#endif
//			std::cout << index++ << "/" << fastPassCurrentQuery->batches.size() << std::endl<< std::endl;
			//Assign the scores to the queries
			fastPassCurrentQuery->assignDotAndBiasFromDevice(&batch, gpuDevice->hostDotMatrix,
					gpuDevice->hostBiasMatrix);

#ifdef DEBUG
			std::cout << "Copying data from device output: ";
			std::cout << std::fixed << std::setprecision(1) << std::fixed << std::setw(7)
			<< ((std::clock() - startBatch) / (double) (CLOCKS_PER_SEC)) * 1000000 << "μs" << std::endl
			<< std::endl;
#endif

			gpuDevice->freeMatrix();

#ifndef DEBUG
			int cur(std::ceil((++currentBatch) / (float)totalNumberOfBatches * PROGRESSBAR_MAX));
			if(cur != lastProgress)
			{
				std::cerr << std::fixed << std::setprecision(2) //
						<< "\r\t[" << std::string(cur, '#') << std::string(PROGRESSBAR_MAX - cur, ' ') << "] " //
						<< currentBatch << "/" << totalNumberOfBatches; //
			}
#endif
		}
#ifndef DEBUG
		//Update progress bar to 100%
		std::cout << std::fixed << std::setprecision(2) //
				<< "\r\t[" << std::string(PROGRESSBAR_MAX, '#') << "] " << totalNumberOfBatches << "/" << totalNumberOfBatches << std::endl;
#endif
		DEF_TIME("Finished all batches on the GPU")

//If we are not testing a variable, proceed with normal output
		if (fastpassParameters->normalOutput)
		{
			/*
			 * Score queries
			 */DEF_STARTTIMER
			{
				fastPassCurrentQuery->scoreQueries();
			}
			DEF_TIME("Scoring queries")

			/*
			 * Open the output file and write the header
			 */DEF_STARTTIMER
			{
				fastPassOutput = new FastPassOutput::PepXML();
				fastPassOutput->setParameters(fastpassParameters);
				fastPassOutput->setOutputFile(queryFile);
				fastPassOutput->addLibraryInfo(fastpassLibraryInfo);
				fastPassOutput->printHeader(fastPassCurrentQuery->getMetaData());
			}
			DEF_TIME("Starting output to " << fastPassOutput->outputFullPathFileName << "")

			/*
			 * Add queries to output
			 */DEF_STARTTIMER
			{
				fastPassOutput->addQueries(&fastPassCurrentQuery->queryPeptides);
			}
			DEF_TIME("Adding queries to output")

			/*
			 * Output queries
			 */DEF_STARTTIMER
			{
				fastPassOutput->printFooter();

			}
			DEF_TIME("Output written to " << fastPassOutput->outputFullPathFileName << "")
			delete fastPassOutput;
		}
		else //We are testing a variable, so no normal output
		{

			if (fastpassParameters->testMinValueVariable == fastpassParameters->testMaxValueVariable)
			{
				std::cout << "Cannot loop over " << fastpassParameters->testVariable
						<< " because the min and max value are equal " << "("
						<< fastpassParameters->testMinValueVariable << " - " << fastpassParameters->testMaxValueVariable
						<< ")" << std::endl;
			}

			//Swap the values if there are not in order
			if (fastpassParameters->testMinValueVariable > fastpassParameters->testMaxValueVariable)
			{
				int temp = fastpassParameters->testMinValueVariable;
				fastpassParameters->testMinValueVariable = fastpassParameters->testMaxValueVariable;
				fastpassParameters->testMaxValueVariable = temp;
			}

			int iterations = (fastpassParameters->testMaxValueVariable - fastpassParameters->testMinValueVariable)
					/ fastpassParameters->testValueVariableStepSize;

			std::string outFile(
					queryFile.substr(0, queryFile.rfind(".")) + ".test." + fastpassParameters->testVariable + ".xls");
			std::ofstream outputFile;

			/*
			 * Time the loop that we are going to run
			 */DEF_STARTTIMER
			{
				std::cout << "Going to loop over '" << fastpassParameters->testVariable << "' " << std::endl;
				std::cout << fastpassParameters->testMinValueVariable << "->"
						<< fastpassParameters->testMaxValueVariable;
				std::cout << " in steps of " << fastpassParameters->testValueVariableStepSize;
				std::cout << " and confidence: fVal >= " << fastpassParameters->testFvalConfidence << std::endl;
				std::cout << "Starting " << iterations << " iterations" << std::endl;
				std::cout << "Outputting in " << outFile << std::endl;

				outputFile.open(outFile, std::ofstream::out | std::ofstream::trunc);

				outputFile << std::endl << "New-value	hitsAboveConfident" << std::endl;
				int hitsAboveConfident;

				for (float newValue = fastpassParameters->testMinValueVariable;
						newValue <= fastpassParameters->testMaxValueVariable;
						newValue += fastpassParameters->testValueVariableStepSize)
				{
					//Reset the counter
					hitsAboveConfident = 0;

					//Update the parameter
					setParameter(fastpassParameters, fastpassParameters->testVariable, std::to_string(newValue));

					//Score the queries again with the changed parameter
					fastPassCurrentQuery->scoreQueries();

					//Retrieve the number of hits above the threshold
					for (auto &query : fastPassCurrentQuery->queryPeptides)
					{
						if (query.matches.size() > 0
								&& query.matches.at(0).fVal >= fastpassParameters->testFvalConfidence)
						{
							hitsAboveConfident++;
						}
					}

					//Output the hits
					outputFile << newValue << "	" << hitsAboveConfident << std::endl;

					{ //Block for the normal pep.xml file output
						/*
						 * Open the outputfile and write the header
						 */DEF_STARTTIMER
						{
							fastPassOutput = new FastPassOutput::PepXML();
							fastPassOutput->setParameters(fastpassParameters);
							fastPassOutput->setOutputFile(
									queryFile + ".test." + fastpassParameters->testVariable + "."
											+ std::to_string(newValue));
							fastPassOutput->addLibraryInfo(fastpassLibraryInfo);
							fastPassOutput->printHeader(fastPassCurrentQuery->getMetaData());
						}
						DEF_TIME("Starting output to " << fastPassOutput->outputFullPathFileName << "")

						/*
						 * Add queries to output
						 */DEF_STARTTIMER
						{
							fastPassOutput->addQueries(&fastPassCurrentQuery->queryPeptides);
						}
						DEF_TIME("Adding queries to output")

						/*
						 * Output queries
						 */DEF_STARTTIMER
						{
							fastPassOutput->printFooter();

						}
						DEF_TIME("Output written to " << fastPassOutput->outputFullPathFileName << "")
						delete fastPassOutput;
					}

				}
				//Close the general outputfile with the threshold hits
				outputFile.close();
			}
			DEF_TIME("Done looping (" << iterations << ")")
		}
		//Delete current  query memory
		delete fastPassCurrentQuery;

		std::cout << "Running at " << totalNumberOfSpectra / ((std::clock() - wallTime) / (double) (CLOCKS_PER_SEC))
				<< " spectra/second" << std::endl;
		std::cout << std::endl << std::endl << std::endl;
	}

	std::cout << "Total filtered spectra processed: " << totalNumberOfSpectra << " in "

	<< fastpassParameters->queryFiles.size() << " files. Taking "
			<< ((std::clock() - wallTime) / (double) (CLOCKS_PER_SEC)) << " seconds -> "
			<< totalNumberOfSpectra / ((std::clock() - wallTime) / (double) (CLOCKS_PER_SEC)) << " spectra/second"
			<< std::endl;

	//Free the memory
	delete gpuDevice;
	delete fastpassLibraryInfo;
	destroyParameters(fastpassParameters);
	std::cout << "Done" << std::endl << std::endl;
	return 0;
}
