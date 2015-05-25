/*
 * fastPassQuery.cpp
 *
 *  Created on: Jan 9, 2015
 *      Author: geforce
 */

#include "fastPassQuery.h"

#include <pugixml/pugixml.hpp>

#include <stringencoders/modp_b64.h>

#include <iostream>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cctype>

namespace FastPassQuery {

Query::Query() :
		numberOfQueries(0)
{
	this->metaData = new QueryMetaData();
}

Query::~Query()
{
	delete this->metaData;
}

bool Query::load(std::string queryFile)
{
	std::string extension(queryFile.substr(queryFile.rfind("."), std::string::npos));
	if (extension == ".mgf")
	{
		this->queryType = mgf;
		return this->loadMgf(queryFile);
	}
	else if (extension == ".mzXML")
	{
		this->queryType = mzxml;
		return this->loadMzXML(queryFile);
	}
	std::cerr << "File extension " << extension << " not recognized, please use .mgf or .mzXML files." << std::endl;
	return false;
}

/**
 * This sorts the queries already before the next sort later on. Reduces a lot of time in the second sorting
 */
int sortOnMass(pugi::xml_node firstNode, pugi::xml_node secondNode)
{
	return atof(firstNode.child_value("precursorMz")) < atof(secondNode.child_value("precursorMz"));
}

typedef union
{
		uint32_t uint;
		float decimal;
} U32;

typedef union
{
		uint64_t uint;
		double decimal;
} U64;

bool Query::loadMzXML(std::string queryFile)
{
	pugi::xml_document doc;

	pugi::xml_parse_result result = doc.load_file(queryFile.c_str());

	if (result)
	{
		std::cout << "pugixml: [" << queryFile << "] parsed without errors, total scans: "
				<< doc.child("mzXML").child("msRun").attribute("scanCount").value() << "\n\n";
	}
	else
	{
		std::cout << "pugixml: [" << queryFile << "] parsed with errors\n";
		std::cout << "Error description: " << result.description() << "\n";
		std::cout << "Error offset: " << result.offset;
		return false;
	}

	pugi::xml_node msRun = doc.child("mzXML").child("msRun");

	if (msRun.child("msInstrument"))
	{
		//mzXML contains some information about the machine used
		for (pugi::xml_node instrumentNode : msRun.child("msInstrument").children())
		{
			if (instrumentNode.attribute("category"))
			{
				//read instrument info in a map for easy acces in the output functions
				this->metaData->instrumentInfoMap.insert(
						std::make_pair(instrumentNode.attribute("category").value(),
								instrumentNode.attribute("value").value()));
			}
		}
		std::cout << "Read instrument info: ";
		for (auto pair : this->metaData->instrumentInfoMap)
		{
			std::cout << pair.first << ":'" << pair.second << "' ";
		}
		std::cout << std::endl;
	}

	/*
	 * First retrieve all the scans
	 */
	std::vector<pugi::xml_node> msLevel2ScanNodes; //Our buffer for the scan pointers
	for (pugi::xml_node ms1scanNode : msRun.children("scan"))
	{
		//Two layered if, could happen that the file includes only msLevel 2 or msLevel 1 with msLevel 2 inside
		if (ms1scanNode.attribute("msLevel").as_int() == 1)
		{
			//Ms1 scan
			for (pugi::xml_node ms2scanNode : ms1scanNode.children("scan"))
			{
				if (ms2scanNode.attribute("msLevel").as_int() == 2)
				{
					//Ms2 scan, the ones we want
					msLevel2ScanNodes.push_back(ms2scanNode);
				}
				else
				{
					//MsX scan (not 1 or 2)
					std::cerr << "Found scan with unkown msLevel: " << ms1scanNode.attribute("msLevel").value()
							<< std::endl;
					continue;
				}

			}
		}
		else if (ms1scanNode.attribute("msLevel").as_int() == 2)
		{
			//Ms2 scan
			msLevel2ScanNodes.push_back(ms1scanNode);
		}
		else
		{
			//MsX scan (not 1 or 2)
			std::cerr << "Found scan with unkown msLevel: " << ms1scanNode.attribute("msLevel").value() << std::endl;
			continue;
		}
	}

	std::cout << "Found " << msLevel2ScanNodes.size() << " msLevel=2 scans" << std::endl;

	/*
	 * Now read them
	 */
	this->queryPeptides.reserve(msLevel2ScanNodes.size()); //reserve the space to prevent allocation delay

	for (pugi::xml_node scan : msLevel2ScanNodes)
	{
		this->queryPeptides.emplace_back();

		this->queryPeptides.back().fileLocation = scan.offset_debug() - 1; //adjusted so it matches the original offset table

		//Check if there is a retention time included, otherwise an empty string is returned
		const char* retentionTime = scan.attribute("retentionTime").as_string("");
		//PTx.xxS format
		if (retentionTime[0] == 'P' && retentionTime[1] == 'T' && retentionTime[strlen(retentionTime) - 1] == 'S')
		{
			//Skip the PT prefix of the seconds by doing pointer arithmetic
			this->queryPeptides.back().retentionTime = atof(retentionTime + 2);
		}

		//Set the precursor mass
		float precursorMass = atof(scan.child_value("precursorMz"));
		this->queryPeptides.back().peptideMass = precursorMass;
		this->queryPeptides.back().intPeptideMass = (int) precursorMass;

		//Strip the path from the query
		std::string::size_type startOffset = queryFile.rfind("/");
		if (startOffset == std::string::npos)
		{
			startOffset = 0;
		}
		else
		{
			startOffset += 1;
		}
		//Construct a TPP style query name QUERY.SCAN.SCAN.CHARGE
		std::stringstream ss;
		ss << queryFile.substr(startOffset, queryFile.rfind("."));
		ss << '.' << std::setw(5) << std::setfill('0');
		ss << std::right << scan.attribute("num").as_int(0);
		ss << '.' << std::setw(5) << std::setfill('0');
		ss << std::right << scan.attribute("num").as_int(0);
		ss << '.' << std::right;
		ss << scan.child("precursorMz").attribute("precursorCharge").as_int(0);

		this->queryPeptides.back().title = ss.str();

		//Read the number of peaks
		this->queryPeptides.back().numberOfPeaks = scan.attribute("peaksCount").as_int(0);

		//Decode the base64 encoded peak data
		std::string base64peaks(scan.child_value("peaks"));
		std::string rawPeaks(modp::b64_decode(base64peaks));

		//Get the precision of the encoded data
		unsigned precision = scan.child("peaks").attribute("precision").as_int(0);

		if (precision != 32 && precision != 64)
		{
			std::cerr << "Unknown precision of the peaks! " << precision << std::endl;
			continue;
		}

		//Be certain that the conventinal network type is used (big endian)
		if (strcmp(scan.child("peaks").attribute("byteOrder").as_string("network"), "network") != 0)
		{
			std::cerr << "byte order of the input mzXML is unkown, "
					<< scan.child("peaks").attribute("byteOrder").as_string("") << ". Skipping this file" << std::endl;
			return false;
		}

		//Check if we are dealing with normal mz-intensity pairs, attribute name was changed between mzXML 2.2 mzXML 3.x
		const char* pairOrder = NULL;
		if (scan.child("peaks").attribute("contentType"))
		{
			pairOrder = scan.child("peaks").attribute("contentType").as_string();
		}
		else if (scan.child("peaks").attribute("pairOrder"))
		{
			pairOrder = scan.child("peaks").attribute("pairOrder").as_string();
		}
		else
		{
			std::cerr << "Unknown pairorder! " << std::endl;
			return false;
		}
		if (pairOrder != NULL && strcmp(pairOrder, "m/z-int") != 0)
		{
			std::cerr << "Content type of the input mzXML is unkown, " << pairOrder << ". Only m/z-int supported"
					<< std::endl;
			return false;
		}

		this->queryPeptides.back().fragments.resize(this->queryPeptides.back().numberOfPeaks);

		//Generate a byte array of the raw data to read the floats from
		const unsigned char* rawData = (unsigned char*) rawPeaks.c_str();
		float floatPrec = 0, floatInten = 0;
		for (int i = 0; i < this->queryPeptides.back().numberOfPeaks; i++)
		{
			//Get the offset of the current peak-intensity pair
			size_t offset = i * 2 * (precision / CHAR_BIT);

			//Read floats
			if (precision == 32)
			{
				U32 mz, intensity;
				/*
				 * Read the mass and intensity from the byte array in big endian format, system endian does not matter
				 * Because we are adding numbers together, the compiler should take care of system endianness
				 */
				mz.uint = ((uint32_t) rawData[3 + offset] << 0) | ((uint32_t) rawData[2 + offset] << 8)
						| ((uint32_t) rawData[1 + offset] << 16) | (rawData[offset] << 24);
				intensity.uint = ((uint32_t) rawData[7 + offset] << 0) | ((uint32_t) rawData[6 + offset] << 8)
						| ((uint32_t) rawData[5 + offset] << 16) | ((uint32_t) rawData[4 + offset] << 24);

				floatPrec = mz.decimal;
				floatInten = intensity.decimal;

			}
			//Read doubles precision == 64
			else
			{
				U64 mz, intensity;
				/*
				 * Read the mass and intensity from the byte array in big endian format, system endian does not matter
				 * Because we are adding numbers together, the compiler should take care of system endianness
				 */
				mz.uint = ((uint64_t) rawData[7 + offset] << 0) | ((uint64_t) rawData[6 + offset] << 8)
						| ((uint64_t) rawData[5 + offset] << 16) | ((uint64_t) rawData[4 + offset] << 24)
						| ((uint64_t) rawData[3 + offset] << 32) | ((uint64_t) rawData[2 + offset] << 40)
						| ((uint64_t) rawData[1 + offset] << 48) | ((uint64_t) rawData[offset] << 56);
				intensity.uint = ((uint64_t) rawData[15 + offset] << 0) | ((uint64_t) rawData[14 + offset] << 8)
						| ((uint64_t) rawData[13 + offset] << 16) | ((uint64_t) rawData[12 + offset] << 24)
						| ((uint64_t) rawData[11 + offset] << 32) | ((uint64_t) rawData[10 + offset] << 40)
						| ((uint64_t) rawData[9 + offset] << 48) | ((uint64_t) rawData[8 + offset] << 56);
				//Convert double to float
				floatPrec = mz.decimal;
				floatInten = intensity.decimal;
			}

			this->queryPeptides.back().fragments.at(i).first = floatPrec;
			this->queryPeptides.back().fragments.at(i).second = floatInten;
		} //peak reading loop

	} //msLevel2nodes loop
	  //Done reading
	return true;
}

bool Query::loadMgf(std::string queryFile)
{
	this->queryFileStream.open(queryFile, std::ifstream::in);

	if (this->queryFileStream.fail())
	{
		this->queryFileStream.close();
		std::cerr << "Could not open the query file: '" << queryFile << "'" << std::endl;
		return false;
	}

	//Our buffer for each line
	std::string line;

	float peakMZ = 0, peakInten = 10; // buffer for peaks
	float tempPeptideMass = 0; //buffer value for type conversion (string to int)
	//Pointer for fast peak reading
	char* ptrLinePtr;
	char* ptrLineCopy = new char[50];

	while (std::getline(this->queryFileStream, line))
	{
		//This 'if' is the slowest of them all..
		if (isdigit(line.c_str()[0])) //start of peaks
		{
			//Faster atof than the C++ or C standard, optimizing for MGF specific
			strcpy(ptrLineCopy, line.c_str());
			ptrLinePtr = ptrLineCopy;

			//Get the precursor mass first
			peakMZ = atof(ptrLineCopy);

			//Skip the whitespace between the two peaks
			while (!std::isspace(*ptrLinePtr) && ptrLinePtr++)
				;
			ptrLinePtr++; //skip the last space
			//Get the intensity of the peak
			peakInten = atoi(ptrLinePtr);

			//Put the peak on the vector
			this->queryPeptides.back().fragments.emplace_back(peakMZ, peakInten);

			this->queryPeptides.back().numberOfPeaks++;
		}
		else if (line.compare(0, 10, "BEGIN IONS") == 0)
		{
			//increase our amount of queries
			this->numberOfQueries++;

			//Put a new QueryPeptide at the end of the vector
			this->queryPeptides.emplace_back();

			//Set the location of the file for quick lookup
			this->queryPeptides.back().fileLocation = this->queryFileStream.tellg();

			//Reserve already some peaks
			this->queryPeptides.back().fragments.reserve(300);
			this->queryPeptides.back().binnedSpectra.reserve(MZ_RANGE);
		}
		else if (line.compare(0, 12, "RTINSECONDS=") == 0)
		{
			this->queryPeptides.back().retentionTime = atof(line.substr(12).c_str());
		}

		else if (line.compare(0, 6, "TITLE=") == 0)
		{
			this->queryPeptides.back().title = line.substr(6);
		}
		else if (line.compare(0, 8, "PEPMASS=") == 0)
		{
			tempPeptideMass = atof(line.c_str() + 8);
			this->queryPeptides.back().peptideMass = tempPeptideMass;
			this->queryPeptides.back().intPeptideMass = (int) tempPeptideMass;
			tempPeptideMass = 0;
		}

	}

	delete[] ptrLineCopy;
	return true;
}

/**
 * peptideMassAsc is a function to sort the queries on peptide mass
 */
int peptideMassAsc(QueryPeptide firstQuery, QueryPeptide secondQuery)
{
	return firstQuery.peptideMass < secondQuery.peptideMass;
}

void Query::sortOnPeptideMass()
{
	std::sort(this->queryPeptides.begin(), this->queryPeptides.end(), peptideMassAsc);
}

/**
 * peaksCompareDesc is a function to sort floats descending order
 */
int peaksCompareDesc(const std::pair<float, float> firstPeak, const std::pair<float, float> secondPeak)
{
	return firstPeak.second > secondPeak.second;
}

/*
 * filterQueries filters the spectra for quality set in the params file
 */

void Query::filterQueries()
{
	float totalIntensity, totalIntensityBelowMZ, maxMz = 0, minMz = 999999;
	unsigned goodPeaks;

	for (auto queryPeptide = this->queryPeptides.begin(); queryPeptide != this->queryPeptides.end(); ++queryPeptide)
	{

		totalIntensity = 0;
		totalIntensityBelowMZ = 0;
		goodPeaks = 0;

		for (auto & fragment : queryPeptide->fragments)
		{
			//Remove peaks with to low intensity or negative values
			if (fragment.second < (float) this->parameters->filterRemovePeakIntensityThreshold)
			{
				fragment.second = 0.0;
				continue;
			}

			if (fragment.second < this->parameters->filterCountPeakIntensityThreshold
					|| fragment.first < this->parameters->filterLightIonsMzThreshold) continue; //Do not count this small peak

			//Increase number of good peaks
			goodPeaks++;
			//Compute total intensity
			totalIntensity += fragment.second;

			if (fragment.first < parameters->filterAllPeaksBelowMz)
			{
				//Compute total intensity below user specified mz
				totalIntensityBelowMZ += fragment.second;
			}
			maxMz = std::max(maxMz, fragment.first);
			minMz = std::min(minMz, fragment.first);
		}

		/*
		 * Filter spectra out that have a low MZ range
		 */
		if (maxMz - minMz < parameters->filterMinMzRange)
		{
			queryPeptide->errorCode = lowMZrange;
			queryPeptide->validSpectra = false;
		}

		/*
		 * Filter spectra out that have not enough peaks
		 */
		if (goodPeaks < this->parameters->filterMinPeakCount)
		{
			queryPeptide->errorCode = notEnoughPeaks;
			queryPeptide->validSpectra = false;
		}

		/*
		 * Check if 95% of the intensity is below the specified m/z
		 */

		if (totalIntensityBelowMZ >= (0.95 * totalIntensity))
		{
			queryPeptide->errorCode = intensityBelowMZ;
			queryPeptide->validSpectra = false;
		}

	}
	for (size_t index = 0; index < this->queryPeptides.size();)
	{
		//Remove any unvalid spectra with efficient swap and pop_back
		if (!this->queryPeptides.at(index).validSpectra)
		{
//std::cerr << "Spectra removed: " << this->queryPeptides.at(index).errorCode << std::endl;
			std::swap(this->queryPeptides.at(index), this->queryPeptides.back()); //swap with the back
			this->queryPeptides.pop_back();     //erase the element
			//No need to go to the next spectra, the swap gives us a new one on the current position
			continue;
		}
		index++;
	}
}

/*
 * Assign library range sets the range of the density struct for each query
 */
void Query::assignLibraryrange()
{
	for (auto & queryPeptide : this->queryPeptides)
	{

		int precursorStart, precursorEnd;

		//assign the precursor mass and add or subtract the tolerance
		precursorStart = queryPeptide.intPeptideMass - (int) (this->parameters->precursorMzTolerance);
		precursorEnd = queryPeptide.intPeptideMass + (int) (this->parameters->precursorMzTolerance);
		//Be sure that the precursor range is within library range
		if (precursorStart < 0 || precursorStart > (int) this->libraryInfo->precursorRange)
		{
			precursorStart = this->libraryInfo->precursorRange - 1;
		}
		if (precursorEnd < 0 || precursorEnd > (int) this->libraryInfo->precursorRange)
		{
			//0 indexing
			precursorEnd = this->libraryInfo->precursorRange - 1;
		}

		queryPeptide.libraryMatchStart = this->libraryInfo->getDensityForPrecursorMass(precursorStart)->firstLocation;
		queryPeptide.libraryMatchEnd = this->libraryInfo->getDensityForPrecursorMass(precursorEnd)->lastLocation;
	}
	for (auto & queryPeptide : this->queryPeptides)
	{
		queryPeptide.matches.reserve(queryPeptide.libraryMatchEnd - queryPeptide.libraryMatchStart + 1);
	}
}

void Query::createBatches()
{
//Some local values
	int paddedLibEntries = 0;
	int paddedQueryEntries = 0;
	uint32_t calculatedMemoryUse = 0;

//Helper variables for easy accessing the query iterator
	Batch *currentBatch;

	this->batches.emplace_back();
	currentBatch = &this->batches.back();

	currentBatch->startQueryBatch = this->queryPeptides.begin();
	bool bValid = true;
	for (auto it = this->queryPeptides.begin(); it != this->queryPeptides.end(); it++)
	{
		int numberOfQueries = std::distance(currentBatch->startQueryBatch, it) + 1;
		int numberOfLibraryEntries = it->libraryMatchEnd - currentBatch->startQueryBatch->libraryMatchStart + 1;

		//The queries need to be padded to multiplies of 16 for the GPU (these will be zeroed)
		paddedLibEntries = numberOfLibraryEntries + (16 - numberOfLibraryEntries % 16);
		paddedQueryEntries = numberOfQueries + (16 - numberOfQueries % 16);

		//Batch queries memory usage on the GPU
		calculatedMemoryUse = ((paddedQueryEntries + paddedLibEntries) * MZ_RANGE
				+ paddedLibEntries * paddedQueryEntries);
		calculatedMemoryUse *= sizeof(float);

		if (this->parameters->precursorMzTolerance > 4
				&& (it->peptideMass - currentBatch->startQueryBatch->peptideMass)
						> 3 * this->parameters->precursorMzTolerance)
		{
			bValid = false;
		}
		else if (calculatedMemoryUse > this->parameters->maxGPUMemoryUsage)
		{
			bValid = false;
		}

		currentBatch->endQueryBatch = it;

		if (!bValid)
		{

			currentBatch->numberOfQueries = std::distance(currentBatch->startQueryBatch, currentBatch->endQueryBatch);
			currentBatch->numberOfLibraryEntries = (currentBatch->endQueryBatch - 1)->libraryMatchEnd
					- currentBatch->startQueryBatch->libraryMatchStart + 1;

			this->batches.emplace_back();
			currentBatch = &this->batches.back();

			currentBatch->startQueryBatch = it;
			bValid = true;
		}
	}

	currentBatch->endQueryBatch = this->queryPeptides.end();

	currentBatch->numberOfQueries = std::distance(currentBatch->startQueryBatch, currentBatch->endQueryBatch);
	currentBatch->numberOfLibraryEntries = (currentBatch->endQueryBatch - 1)->libraryMatchEnd
			- currentBatch->startQueryBatch->libraryMatchStart + 1;

}

void Query::binSpectraQueries()
{
	int intMZ = 0;
	float floatMZ = 0, floatInten = 0;
	float basePeakIntensity;

	for (auto & query : this->queryPeptides)
	{
		//Reserve the space for the peaks, ranging from MZ_START - MZ_END (=MZ_RANGE) 0 - (MZ_END-MZ_START)
		query.binnedSpectra.assign(MZ_RANGE, 0);

		basePeakIntensity = 0;

		//First calculate some values such as the base peak, that is necessary for the peak settings
		for (auto & fragment : query.fragments)
		{
			floatMZ = fragment.first;

			//Remove light peaks
			if (floatMZ < (float) this->parameters->filterLightIonsMzThreshold)
			{
				fragment.second = 0.0;
			}

			//Remove the peaks around the precursor
			if (floatMZ > (query.peptideMass - 60) && floatMZ < (query.peptideMass + 20))
			{
				fragment.second = 0.0;
			}

			//Initialize scaledIntensity
			if (this->parameters->peakScalingMzPower == 0.0 && this->parameters->peakScalingIntensityPower == 0.5)
			{
				fragment.second = sqrt(fragment.second); //Shortcut to sqrt Intensity
			}
			else
			{
				fragment.second = pow(floatMZ, this->parameters->peakScalingMzPower)
						* pow(fragment.second, this->parameters->peakScalingIntensityPower);
			}
		}

		/*
		 * Sort the peaks of a query and constrain the vector to the specified number of peaks by the user
		 */
		std::sort(query.fragments.begin(), query.fragments.end(), peaksCompareDesc);
		query.fragments.resize(parameters->filterMaxPeaksUsed);

		basePeakIntensity = query.fragments.at(0).second;
		//Now set the values of the peaks in the corresponding bins
		for (auto & fragment : query.fragments)
		{
			//Adjust the start so the vector start at 0 instead of MZ_START (default:10)
			intMZ = (int) fragment.first - MZ_START;
			floatMZ = fragment.first;
			floatInten = fragment.second;

			//Remove peak if it does not fit the set dynamic range
			if (floatInten < (basePeakIntensity / this->parameters->filterMaxDynamicRange))
			{
				floatInten = 0;
				continue;
			}

			//Check if the peak is in bounds
			if (intMZ < 0)
			{
				//Out of bound, set lowest peaks
				query.binnedSpectra[0] += floatInten;
				query.binnedSpectra[1] += floatInten * this->parameters->peakBinningFractionToNeighbor;
			}
			else if (intMZ > (int) (query.binnedSpectra.size() - 2))
			{
				//Out of bound, set highest two peaks
				query.binnedSpectra[query.binnedSpectra.size() - 2] += floatInten
						* this->parameters->peakBinningFractionToNeighbor;
				query.binnedSpectra[query.binnedSpectra.size() - 1] += floatInten;
			}
			else
			{
				//Set the intensity
				query.binnedSpectra[intMZ - 1] += floatInten * this->parameters->peakBinningFractionToNeighbor;
				query.binnedSpectra[intMZ] += floatInten;
				query.binnedSpectra[intMZ + 1] += floatInten * this->parameters->peakBinningFractionToNeighbor;
			}
		}

		/*
		 * In SpectraST this step happens after the dot product, but can be performed now already?
		 *
		 */
		float squareMagnitude = 0, magnitude = 0;
		for (auto & peakIntensity : query.binnedSpectra)
		{
			squareMagnitude += peakIntensity * peakIntensity;
		}

		magnitude = sqrt(squareMagnitude);

		if (magnitude > 0)
		{
			for (auto & peakIntensity : query.binnedSpectra)
				peakIntensity /= magnitude;
		}

	}

}

/*
 * Assign the results from the GPU to the queries, according to the library entries
 */
void Query::assignDotAndBiasFromDevice(Batch *currentBatch, float* dotMatrix, float* biasMatrix)
{
	int queryIndex = 0;

	for (auto query = currentBatch->startQueryBatch; query != currentBatch->endQueryBatch; query++, queryIndex++)
	{

		//std::cout << query->libraryMatchStart << " - " << query->libraryMatchEnd << std::endl;
		for (int libindex = query->libraryMatchStart; libindex <= query->libraryMatchEnd; libindex++)
		{
			int index = (libindex - currentBatch->startQueryBatch->libraryMatchStart) * currentBatch->numberOfQueries
					+ queryIndex;

			/*
			 * Double check if our mass is correct (we do not do this earlier when loading the library,
			 * because then the advantage of reading the library sequentially is gone.
			 * This check is based on the floating point mass and not on integer mass
			 */
			if (fabs(query->peptideMass - this->libraryInfo->getLibraryPeptide(libindex)->PrecursorMZ)
					<= this->parameters->precursorMzTolerance /*&& dotMatrix[index] > 0.001f*/ ) //Filter out small dot values
			{

				//		std::cout << "(" << libindex << " - " << currentBatch->startQueryBatch->libraryMatchStart << ") " << (libindex - currentBatch->startQueryBatch->libraryMatchStart)<< "x"
				//				<< queryIndex << "	" << currentBatch->numberOfQueries << " -> " << index << "/"
				//				<< currentBatch->numberOfQueries * currentBatch->numberOfLibraryEntries << " \n";

				query->matches.emplace_back(libindex, dotMatrix[index], biasMatrix[index]);

				//Count the number of hits (dot > 0.01) for later score calculation
				if (dotMatrix[index] > 0.01)
				{
					query->numHits++;
				}
			}

		}
//		query->binnedSpectra.clear();
//		query->binnedSpectra.shrink_to_fit();
//		query->matches.shrink_to_fit();
	}
}

/**
 * peptideMassDesc is a function to sort the queries on dot prod
 */
int dotProdDesc(Candidate firstCandidate, Candidate secondCandidate)
{
	return firstCandidate.dotProd > secondCandidate.dotProd;
}
/**
 * fValDesc is a function to sort the queries on fValue
 */

int fValDesc(Candidate firstCandidate, Candidate secondCandidate)
{
	return firstCandidate.fVal > secondCandidate.fVal;
}

/*
 * Calculates the score of each query based on the gathered dot and dotBias
 */
void Query::scoreQueries()
{
	for (auto & query : this->queryPeptides)
	{

		std::sort(query.matches.begin(), query.matches.end(), dotProdDesc);

		if (this->parameters->detectHomologs > 1)
		{
			this->detectHomologs(&query);
		}
		else
		{
			if (query.matches.size() > 1)
			{
				query.matches[0].firstNonHomolog = 2;
			}
		}

		double totalDot = 0, totalSqDot = 0;
		for (auto & candidate : query.matches)
		{
			if (candidate.firstNonHomolog > 0)
			{
				candidate.deltaDot = candidate.dotProd - query.matches.at(candidate.firstNonHomolog - 1).dotProd;
			}
			this->calcFval(&candidate, query.numHits);
			totalDot += candidate.dotProd;
			totalSqDot += candidate.dotProd * candidate.dotProd;
		}

		query.hitsMean = totalDot / (double) query.numHits;
		query.hitsDotStdev = totalSqDot / (double) query.numHits - query.hitsMean * query.hitsMean;
		if (query.hitsDotStdev > 0.000001)
		{
			query.hitsDotStdev = sqrt(query.hitsDotStdev);
		}
		else if (query.hitsDotStdev < 0) //rare case of -0 value
		{
			query.hitsDotStdev = 0;
		}

		std::sort(query.matches.begin(), query.matches.end(), fValDesc);

	}
}

/*
 * Calcuates the f value
 */

void Query::calcFval(Candidate *candidate, int totalHits)
{

	if (candidate->dotProd < 0.00001)
	{
		candidate->fVal = -0.00001;
		return;
	}

	float normDelta = candidate->deltaDot / candidate->dotProd;

	if (normDelta > 2.0 * candidate->dotProd)
	{
		normDelta = 2.0 * candidate->dotProd;
	}

	float fractionDelta = this->parameters->fvalFractionDelta;

	candidate->fVal = (1 - fractionDelta) * candidate->dotProd + fractionDelta * normDelta;

	if (candidate->fVal > 0.1 && totalHits < 20)
	{
		candidate->fVal = (1 - 0.5 * fractionDelta) * candidate->dotProd;
	}

	if (this->parameters->fvalUseDotBias && candidate->dotBias >= 0.000001)
	{
// impose dot bias penalty
		if (candidate->dotBias < 0.09)
		{
			candidate->fVal -= 0.12;
		}
		else if (candidate->dotBias > 0.32 && candidate->dotBias <= 0.35)
		{
			candidate->fVal -= (candidate->dotBias - 0.32) * 4.0;
		}
		else if (candidate->dotBias > 0.35 && candidate->dotBias <= 0.45)
		{
			candidate->fVal -= (0.12 + (candidate->dotBias - 0.35) * 1.2);
		}
		else if (candidate->dotBias > 0.45)
		{
			candidate->fVal -= 0.24;
		}

		if (candidate->fVal <= 0.0)
		{
			candidate->fVal = -0.00001;
		}
	}
}

bool equalPeptides(FastPassLibrary::Library_peptide* first, FastPassLibrary::Library_peptide* second)
{
	std::string stripedPeptideFirst(first->FullName);
	std::map<int, std::string> modsFirst;
	unsigned int chargeFirst = atoi(
			stripedPeptideFirst.substr(stripedPeptideFirst.length() - 1, std::string::npos).c_str());

	std::string stripedPeptideSecond(second->FullName);
	std::map<int, std::string> modsSecond;
	unsigned int chargeSecond = atoi(
			stripedPeptideSecond.substr(stripedPeptideSecond.length() - 1, std::string::npos).c_str());

	std::string::size_type currentPos;
	while ((currentPos = stripedPeptideFirst.find('[')) != std::string::npos)
	{
		//Keep track of the modifications
		modsFirst.insert(std::make_pair((int) currentPos, stripedPeptideFirst.substr(currentPos - 1, 1)));

		//N or C terminal modification, delete the marker also
		if (currentPos > 0
				&& (stripedPeptideFirst.at(currentPos - 1) == 'n' || stripedPeptideFirst.at(currentPos - 1) == 'c'))
		{
			currentPos -= 1;
		}
		stripedPeptideFirst.erase(currentPos, stripedPeptideFirst.find(']', currentPos) - currentPos + 1);
	}

	while ((currentPos = stripedPeptideSecond.find('[')) != std::string::npos)
	{
		//Keep track of the modifications
		modsFirst.insert(std::make_pair((int) currentPos, stripedPeptideSecond.substr(currentPos - 1, 1)));

		//N or C terminal modification, delete the marker also
		if (currentPos > 0
				&& (stripedPeptideSecond.at(currentPos - 1) == 'n' || stripedPeptideSecond.at(currentPos - 1) == 'c'))
		{
			currentPos -= 1;
		}
		stripedPeptideSecond.erase(currentPos, stripedPeptideSecond.find(']', currentPos) - currentPos + 1);
	}

	if (stripedPeptideFirst != stripedPeptideSecond)
	{
		/*  if (considerIsobaricAASame) {
		 // if considerIsobaricAASame is set, this function will consider K and Q to be the same,
		 // and I and L to be the same. The degenerate symbols B (for N and D) and Z (for Q and E)
		 // will be considered the same as their nondegenerate ones.

		 // check the length of the peptides
		 if (this->stripped.length() != p.stripped.length()) {
		 return false;
		 }

		 // they are the same length, now go through each AA to compare.
		 string::size_type i, j;

		 for (i = 0, j = 0; i < this->stripped.length() && j < p.stripped.length(); i++, j++) {
		 if (!((this->stripped[i] == p.stripped[j]) ||
		 (this->stripped[i] == 'K' && p.stripped[j] == 'Q') ||
		 (this->stripped[i] == 'Q' && p.stripped[j] == 'K') ||
		 (this->stripped[i] == 'I' && p.stripped[j] == 'L') ||
		 (this->stripped[i] == 'L' && p.stripped[j] == 'I') ||
		 (this->stripped[i] == 'B' && p.stripped[j] == 'N') ||
		 (this->stripped[i] == 'B' && p.stripped[j] == 'D') ||
		 (this->stripped[i] == 'N' && p.stripped[j] == 'B') ||
		 (this->stripped[i] == 'D' && p.stripped[j] == 'B') ||
		 (this->stripped[i] == 'Z' && p.stripped[j] == 'Q') ||
		 (this->stripped[i] == 'Z' && p.stripped[j] == 'E') ||
		 (this->stripped[i] == 'Q' && p.stripped[j] == 'Z') ||
		 (this->stripped[i] == 'E' && p.stripped[j] == 'Z'))) {
		 return false;
		 }
		 }

		 } else {*/
		return false;
//}
	}
	if (chargeFirst != 0 && chargeSecond != 0 && chargeFirst != chargeSecond)
	{
		return false;
	}
//	if (this->hasUnknownMod || p.hasUnknownMod)
//	{
//		return false;
//	}
	if (modsFirst != modsSecond)
	{
		return false;
	}
//	if (this->isModsSet && p.isModsSet && ((this->nTermMod != p.nTermMod) || (this->cTermMod != p.cTermMod)))
//	{
//		return false;
//	}
	return true;
}

bool isSubsequence(FastPassLibrary::Library_peptide* first, FastPassLibrary::Library_peptide* second, bool ignoreMods)
{
	std::string peptideFirst(first->FullName);
	std::string peptideSecond(second->FullName);
	if (ignoreMods)
	{
		peptideFirst.assign(first->StripedName);
		peptideSecond.assign(second->StripedName);
	}

	std::string::size_type found;
	if (peptideFirst.length() < peptideSecond.length())
	{
		found = peptideSecond.find(peptideFirst, 0);
	}
	else if (peptideFirst.length() > peptideSecond.length())
	{
		found = peptideFirst.find(peptideSecond, 0);
	}
	else
	{
		//Equal length
		return peptideFirst == peptideSecond;
	}
	return (found != std::string::npos);

}

//From SpectraST
// isHomolog - see if the two peptide sequences are homologous.
// use a dynamic programming algorithm to align the two sequences
// At the end, T[m][n] will contain a measure of sequence identity
// Basically, an aligned AA will get 1 point, and an aligned AA following another aligned AA will get 2 points.
// Examples:
//    DEFGHI-
// vs DEFGHIK  (1+2+2+2+2+2 = 11 points)
//    DEFGGHI-
// vs DEFG-HIK (1+2+2+2+1+2 = 10 points)
//    DEFFG-HI-
// vs YEF-GGHIK (1+2+1+1+2 = 7 points)
bool isHomolog(FastPassLibrary::Library_peptide* first, FastPassLibrary::Library_peptide* second, double threshold)
{
	std::string peptideFirst(first->StripedName);
	std::string peptideSecond(second->StripedName);

	unsigned int m = (unsigned int) (peptideFirst.length());
	unsigned int n = (unsigned int) (peptideSecond.length());

	std::vector<std::string> aa1(m + 1);
	std::vector<std::string> aa2(n + 1);

	std::string::size_type shorter = m;
	if (m > n) shorter = n;

	// if the shorter sequence is a perfect sub-sequence of the longer one, its score is (shorter * 2 - 1).
	// so we are calculating our threshold identity score by multiplying the specified threshold by this "perfect" score
	double minIdentity = threshold * (shorter * 2 - 1);

	int** alignedLength = new int*[m + 1];
	int** identical = new int*[m + 1];

	unsigned int i = 0;
	unsigned int j = 0;

	for (i = 0; i <= m; i++)
	{
		alignedLength[i] = new int[n + 1];
		identical[i] = new int[n + 1];

		alignedLength[i][0] = 0;
		if (i > 0) aa1[i] = peptideFirst[i - 1];

		for (j = 0; j <= n; j++)
		{
			identical[i][j] = 0;
		}
	}

	for (j = 0; j <= n; j++)
	{
		alignedLength[0][j] = 0;
		if (j > 0) aa2[j] = peptideSecond[j - 1];
	}

	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
		{
			if (aa1[i] == aa2[j] || (aa1[i] == "L" && aa2[j] == "I")
					|| // we consider I and L the same, K and Q the same
					(aa1[i] == "I" && aa2[j] == "L") || (aa1[i] == "K" && aa2[j] == "Q")
					|| (aa1[i] == "Q" && aa2[j] == "K"))
			{
				alignedLength[i][j] = alignedLength[i - 1][j - 1] + 1 + identical[i - 1][j - 1];
				identical[i][j] = 1;
			}
			else
			{
				if (alignedLength[i - 1][j] > alignedLength[i][j - 1])
				{
					alignedLength[i][j] = alignedLength[i - 1][j];
				}
				else
				{
					alignedLength[i][j] = alignedLength[i][j - 1];
				}

			}
		}
	}

//identity = alignedLength[m][n];
	bool result = ((double) (alignedLength[m][n]) >= minIdentity);

	for (unsigned int r = 0; r <= m; r++)
	{
		delete[] alignedLength[r];
		delete[] identical[r];
	}
	delete[] alignedLength;
	delete[] identical;

	return (result);
}

//Detects homologs for a specific query
void Query::detectHomologs(QueryPeptide *query)
{

	if (query->matches.empty() || query->matches.size() == 1) return;

	FastPassLibrary::Library_peptide* topHit = this->libraryInfo->getLibraryPeptide(query->matches[0].libId);

	bool homologFound = false;
	unsigned int curRank = 0;

// go down the hit lists until hitting a nonhomologous hit
	do
	{
		homologFound = false;
		curRank++;

		FastPassLibrary::Library_peptide* thisHit = this->libraryInfo->getLibraryPeptide(query->matches[curRank].libId);

		//std::cout << topHit->FullName << " =?= " << thisHit->FullName;

		if (equalPeptides(topHit, thisHit))
		{
			//std::cout << " Equal peptides!";
			// identical!
			homologFound = true;
			//}else if (fabs(topHit->PrecursorMZ - thisHit->PrecursorMZ) < 0.01)
			//{
			//  // exactly same mass, consider as homolog (possibly target/decoy pair)
			//  homologFound = true;
		}
		else if (isSubsequence(topHit, thisHit, true))
		{
			//std::cout << " Subsequence!";
			// one is subsequence of the other!
			homologFound = true;
		}
		else if (isHomolog(topHit, thisHit, 0.7))
		{
			//std::cout << " Homolog!";
			homologFound = true;

		}
		else
		{
			//std::cout << "No homolog";
		}
		//std::cout << std::endl;
	} while (homologFound && curRank < this->parameters->detectHomologs - 1
			&& curRank < (unsigned int) (query->matches.size()) - 1);
//std::cout << "Setting firstNonHomolog to: "<< this->libraryInfo->getLibraryPeptide(query->matches[curRank].libId)->FullName;
//std::cout << std::endl << std::endl;

// setting the field firstNonHomolog for all the homologs found
	for (unsigned int rank = 0; rank < curRank; rank++)
	{
		query->matches[rank].firstNonHomolog = curRank + 1;
	}

	return;
}

} //namespace FastPassQuery

