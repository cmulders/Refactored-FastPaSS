/*
 * fastPassLibrary.c++
 *
 *  Created on: Jan 9, 2015
 *      Author: geforce
 */

#include "fastPassLibrary.h"
#include <iostream>
#include <cmath>
#include <string>
#include "assert.h"
#include <algorithm>

namespace FastPassLibrary {

Info::Info() :
		numberOfEntries(0), precursorRange(0)
{
}

Info::~Info()
{
	this->libraryFileStream.close();
}
/*
 * Read the library into the internal variables and read the density struct into a vector
 */
bool Info::readFPLibrary(std::string libraryFile)
{
	this->libraryFile = libraryFile;

	this->libraryFileStream.open(this->libraryFile, std::ifstream::binary);
	if (this->libraryFileStream.fail())
	{
		this->libraryFileStream.close();
		std::cerr << "Could not open the library file: '" << this->libraryFile << "'" << std::endl;
		return false;
	}

	this->libraryFileStream.read(reinterpret_cast<char*>(&this->numberOfEntries), sizeof(this->numberOfEntries));
	this->libraryFileStream.read(reinterpret_cast<char*>(&this->precursorRange), sizeof(this->precursorRange));

	this->densityList.reserve(this->precursorRange);

	//Read the density structs into our library info
	Density *densityStruct = new Density();
	for (uint32_t i = 0; i < this->precursorRange; i++)
	{
		this->libraryFileStream.read(reinterpret_cast<char*>(densityStruct), sizeof(Density));
		this->densityList.push_back(*densityStruct);
	}
	delete densityStruct;

	//Skip the library spectra and start reading after that again
	this->libraryFileStream.seekg(this->numberOfEntries * MZ_RANGE * sizeof(float), std::ios_base::cur);

	this->libraryPeptides.reserve(this->numberOfEntries);
	//Read the density structs into our library info
	Library_peptide *libraryPeptidesStruct = new Library_peptide();
	for (uint32_t i = 0; i < this->numberOfEntries; i++)
	{
		this->libraryFileStream.read(reinterpret_cast<char*>(libraryPeptidesStruct), sizeof(Library_peptide));
		this->libraryPeptides.push_back(*libraryPeptidesStruct);
	}
	delete libraryPeptidesStruct;

	return true;
}

/*
 * Print general info about the library
 */
void Info::printLibraryInfo()
{
	std::cout << "Library info" << std::endl;
	std::cout << "	Library file: " << this->libraryFile;
	std::cout << " Number of entries:" << this->numberOfEntries;
	std::cout << ", ";
	std::cout << "Precursor range:" << this->precursorRange;
	std::cout << std::endl;
}

/*
 * Retrieves density information for a specific mass
 * Because the density vector contains the density's according to their mass, the indexing is easy
 */
Density *Info::getDensityForPrecursorMass(int mass)
{
	if (mass < 0 || (size_t) mass > (this->densityList.size() - 1))
	{
		std::cout << "Reading out of bounds! " << 0 << "<=" << mass << " < " << this->densityList.size() << std::endl;
		return new Density(0, 0, 0);
	}
	return &this->densityList.at(mass);
}

/*
 * Retrieves the library peptide information for a specific Lib index
 * Because the Library_peptide vector contains the peptides according to their lib index, the indexing is easy
 */
Library_peptide *Info::getLibraryPeptide(int libIndex)
{
	if (libIndex > (int) (this->libraryPeptides.size() - 1))
	{
		std::cout << "Out of bounds!: " << libIndex << std::endl;
		abort();
	}
	return &this->libraryPeptides.at(libIndex);
}

/*
 * Gets the peaks of the library entries for a specific range
 */
bool Info::readPeaksForLibraryRangeIntoPointer(int libStart, int numberOfSpectra, float *pointer)
{
	//Reposition the seek pointer to the start of the library spectra
	this->libraryFileStream.seekg(
			sizeof(this->numberOfEntries) + sizeof(this->precursorRange) + this->precursorRange * sizeof(Density),
			std::ios_base::beg);

	//Skip the other library spectra
	this->libraryFileStream.seekg(libStart * MZ_RANGE * sizeof(float), std::ios_base::cur);

	//Read the current library spectra (multiple at once)
	this->libraryFileStream.read((char*) pointer, numberOfSpectra * MZ_RANGE * sizeof(float));
	return true;
}

Converter::Converter()
{

}

Converter::~Converter()
{
	entries.clear();
	libraryPeptides.clear();
	densityList.clear();
}

/*
 * From SpectraST codebase: FileUtils.cpp::289
 */
bool nextLine(std::istream& fileStream, std::string &line)
{
	if (fileStream.eof())
	{
		line = "_EOF_";
		return (false);
	}

	std::getline(fileStream, line);
	if (line.empty()) return (true);

	std::string::size_type last = line.length() - 1;
	if (line[last] == '\r')
	{
		line.erase(last);
	}
	return (true);
}

/**
 * peaksCompareDesc is a function to sort floats descending order
 */
int peaksCompareDesc(const Peak firstPeak, const Peak secondPeak)
{
	return firstPeak.intensity > secondPeak.intensity;
}

/**
 * peaksCompareDesc is a function to sort floats descending order
 */
int entriesCompareDesc(const SPLibEntry firstEntry, const SPLibEntry secondEntry)
{
	return firstEntry.precurosrMz < secondEntry.precurosrMz;
}

void Converter::convert(std::string splibfile)
{
	std::string::size_type extension = splibfile.rfind(".");
	if (splibfile.substr(extension) != ".splib")
	{
		std::cerr << "Library file is not an .splib file " << splibfile << " skipping this file" << std::endl;
		return;
	}
	std::string outLibFile(splibfile);
	outLibFile.replace(extension, std::string::npos, Parameters::libraryExtension);

	std::ifstream libraryFileStream(splibfile, std::ifstream::binary | std::ifstream::in);
	if (libraryFileStream.fail())
	{
		libraryFileStream.close();
		std::cerr << "Could not open the library file: '" << splibfile << "'" << std::endl;
		return;
	}

	std::ofstream libraryOutStream(outLibFile, std::ofstream::binary | std::ofstream::out | std::ofstream::trunc);
	if (libraryOutStream.fail())
	{
		libraryOutStream.close();
		std::cerr << "Could not open the output library file: '" << outLibFile << "'" << std::endl;
		return;
	}

	std::cout << "Converting " << splibfile << " to " << outLibFile << std::endl;

	char firstChar = libraryFileStream.peek();
	if (firstChar == '#' || firstChar == 'N')
	{
		//sptxt file not binary
		libraryFileStream.close();
		std::cerr << "This library file is not the binary .splib but the plaintext file .sptxt" << std::endl;
		return;
	}

	//Retrieve the file size
	libraryFileStream.seekg(0, libraryFileStream.end);
	std::istream::pos_type fileSize = libraryFileStream.tellg();
	libraryFileStream.seekg(0, libraryFileStream.beg);

	/*
	 * Preamble
	 */
	int spectrastVersion = 0;
	int spectrastSubVersion = 0;
	unsigned int numLines = 0;
	std::string bufferLine;

	libraryFileStream.read((char*) (&spectrastVersion), sizeof(int));
	libraryFileStream.read((char*) (&spectrastSubVersion), sizeof(int));

	if (!nextLine(libraryFileStream, bufferLine))
	{
		//splib is not formatted properly
		libraryFileStream.close();
		std::cerr << "Corrupt .splib." << std::endl;
		return;
	}

	std::string fileName(bufferLine);

	std::string firstLine("");

	libraryFileStream.read((char*) (&numLines), sizeof(unsigned int));

	for (unsigned int i = 0; i < numLines; i++)
	{
		if (!nextLine(libraryFileStream, bufferLine))
		{
			//splib is not formatted properly
			libraryFileStream.close();
			std::cerr << "Corrupt .splib." << std::endl;
			return;
		}
		if (firstLine.empty())
		{
			firstLine = fileName + " : " + bufferLine;
			splibPreamble.push_back("> " + firstLine);
		}
		else
		{
			splibPreamble.push_back((bufferLine[0] == '>' ? "" : "> ") + bufferLine);
		}
	}

	//Print some information for the user
	std::cout << "> Library created by SpectraST " << spectrastVersion << "." << spectrastSubVersion << std::endl;
	for (auto &line : splibPreamble)
	{
		std::cout << line << std::endl;
	}

	std::cout << std::endl << "Start reading splib into memory" << std::endl;
	/*
	 * Entries
	 */
	std::istream::pos_type lastProgressOffset = 0;
	do
	{
		entries.emplace_back();

		entries.back().offset = libraryFileStream.tellg();

		//Read general info about the peptide entry
		libraryFileStream.read((char*) (&entries.back().libId), sizeof(entries.back().libId));
		nextLine(libraryFileStream, entries.back().fullName);
		libraryFileStream.read((char*) (&entries.back().precurosrMz), sizeof(entries.back().precurosrMz));
		nextLine(libraryFileStream, entries.back().status);

		//Read the peaks
		libraryFileStream.read((char*) (&entries.back().numPeaks), sizeof(entries.back().numPeaks));
		for (unsigned int i = 0; i < entries.back().numPeaks; i++)
		{
			//add new peak
			entries.back().peaks.emplace_back();

			//add the peak information
			libraryFileStream.read((char*) (&entries.back().peaks.back().mz), sizeof(Peak::mz));
			libraryFileStream.read((char*) (&entries.back().peaks.back().intensity), sizeof(Peak::intensity));
			nextLine(libraryFileStream, entries.back().peaks.back().annotation);
			nextLine(libraryFileStream, entries.back().peaks.back().info);
		}
		//Reserve the space for the binned spectra
		entries.back().libSpectra.assign(MZ_RANGE, 0);

		//Process the command line to a map for easy acces
		std::string tempComments;
		nextLine(libraryFileStream, tempComments);

		bool inQuotes = false;
		std::string::size_type startOffset = 0, curOffset;
		std::string oneComment;
		for (std::string::iterator it = tempComments.begin(); it != tempComments.end(); it++)
		{
			//Spaces in quotes we skip
			if ((*it) == '"') inQuotes = !inQuotes;

			if ((*it) == ' ' && !inQuotes)
			{
				curOffset = std::distance(tempComments.begin(), it);
				oneComment.assign(tempComments.substr(startOffset, curOffset - startOffset));

				entries.back().comments[oneComment.substr(0, oneComment.find("="))] = oneComment.substr(
						oneComment.find("=") + 1, std::string::npos);

				startOffset = curOffset + 1;
			}
		}

		//Progress bar for the user
		if ((libraryFileStream.tellg() - lastProgressOffset) > fileSize / 10)
		{
			lastProgressOffset = libraryFileStream.tellg();
			std::cout << (int) (((float) lastProgressOffset / fileSize) * 100) << "%.." << std::flush;
		}
	} while (libraryFileStream.peek() >= 0 && !libraryFileStream.eof());
	std::cout << "100%" << std::endl << "Done reading splib file" << std::endl << std::endl;

	//Process the entries and filter them
	processLibEntries();

	//Bin the spectra in 1Da bins
	binLibEntries();

	std::cout << "Start sorting library entries (may take a while, without visual indication)" << std::endl;
	std::sort(entries.begin(), entries.end(), entriesCompareDesc);
	std::cout << "Done sorting" << std::endl << std::endl;

	int processed = 0;
	std::cout << "Creating density array" << std::endl;

	densityList.clear();
	densityList.reserve(MZ_END);
	//Start of density list (mass 0, so all zero), because we use zero-indexing
	densityList.emplace_back(0, 0, 0);

	//Create all the densities
	for (int i = 1; i <= MZ_END; i++)
	{
		densityList.emplace_back(0, UINT32_MAX, 0);
	}

	//populate the density structs
	for (size_t index = 0; index < entries.size(); index++)
	{
		//std::cout << entries.back().comments.begin()->second << std::endl;
		if ((int) (index - processed) > (int) entries.size() / 10)
		{
			processed = index;
			std::cout << (int) (((float) processed / entries.size()) * 100) << "%.." << std::flush;
		}

		//Restrict mass to (0)-(MZ_END-1)
		int mass = std::min(MZ_END - 1, std::max(0, (int) entries.at(index).precurosrMz));

		densityList.at(mass).numberOfEntries++;
		densityList.at(mass).firstLocation = std::min(densityList.at(mass).firstLocation, (uint32_t) index);
		densityList.at(mass).lastLocation = std::max(densityList.at(mass).lastLocation, (uint32_t) index);
	}
	//Relink the density structs so it is one big continuous list
	for (int i = 1; i <= MZ_END; i++)
	{
		//If the previous entry and this entry differ by more than one, reset this entry
		if (std::abs(densityList.at(i).firstLocation - densityList.at(i - 1).lastLocation) > 1)
		{
			densityList.at(i).firstLocation = densityList.at(i - 1).lastLocation;
		}

		//If the current first location is higher than the last location, set the alst lcoation to the first location
		//Resulting in a single entry
		if (densityList.at(i).firstLocation > densityList.at(i).lastLocation)
		{
			densityList.at(i).lastLocation = densityList.at(i).firstLocation;
		}
	}
	std::cout << "100%" << std::endl << "Done creating density array" << std::endl << std::endl;

	/*
	 * Final write to file
	 */
	std::cout << "Writing peptides" << std::endl;
	uint32_t numEntries = entries.size();
	uint32_t precursorRange = densityList.size();
	libraryOutStream.write((char*) (&numEntries), sizeof(numEntries));
	libraryOutStream.write((char*) (&precursorRange), sizeof(precursorRange));
	libraryOutStream.write((char*) (densityList.data()), sizeof(densityList[0]) * densityList.size());

	//Creating the real lib entries for the binary file with fixed lengths (for fast retrieval)
	std::size_t copiedCharsLength; //Stores the char pointer for the 0 byte of the strings
	for (auto &entry : entries)
	{
		libraryOutStream.write((char*) (entry.libSpectra.data()), sizeof(entry.libSpectra[0]) * MZ_RANGE);

		//Get a new peptide
		libraryPeptides.emplace_back();

		libraryPeptides.back().LibID = entry.libId;

		std::string stripedName(entry.fullName.substr(2,entry.fullName.rfind(".")-2));
		std::string::size_type currentPos;
		while ((currentPos = stripedName.find('[')) != std::string::npos)
		{
			//N or C terminal modification, delete the marker also
			if (currentPos > 0 && (stripedName.at(currentPos - 1) == 'n' || stripedName.at(currentPos - 1) == 'c'))
			{
				currentPos -= 1;
			}
			stripedName.erase(currentPos, stripedName.find(']', currentPos) - currentPos + 1);
		}
		copiedCharsLength = stripedName.copy(libraryPeptides.back().StripedName, 255);
		libraryPeptides.back().StripedName[copiedCharsLength] = '\0';

		copiedCharsLength = entry.fullName.copy(libraryPeptides.back().FullName, 255);
		libraryPeptides.back().FullName[copiedCharsLength] = '\0';

		copiedCharsLength = entry.comments.at("Protein").copy(libraryPeptides.back().protein, 255);
		libraryPeptides.back().protein[copiedCharsLength] = '\0';

		//Check if the remark comment exists
		if (entry.comments.count("Remark") == 0) entry.comments["Remark"] = "_NONE_";
		copiedCharsLength = entry.comments.at("Remark").copy(libraryPeptides.back().libRemark, 99);
		libraryPeptides.back().libRemark[copiedCharsLength] = '\0';

		//Copy over the status line
		copiedCharsLength = entry.status.copy(libraryPeptides.back().libStatus, 99);
		libraryPeptides.back().libStatus[copiedCharsLength] = '\0';

		//Precursor MZ and number of epaks
		libraryPeptides.back().PrecursorMZ = entry.precurosrMz;
		libraryPeptides.back().NumPeaks = entry.numPeaks;

		//The bianry position in the original library
		libraryPeptides.back().splibPos = entry.offset;

		//Library probility and num of replicates
		libraryPeptides.back().libProbability = atof(entry.comments["Prob"].c_str());
		libraryPeptides.back().libReplicates = atoi(entry.comments["Nreps"].c_str());
	}
	//Write all the library peptides in one chunk
	libraryOutStream.write((char*) (libraryPeptides.data()), sizeof(libraryPeptides[0]) * libraryPeptides.size());
	libraryOutStream.close();
	std::cout << "Done writing peptides" << std::endl << std::endl;

	//Cleanup the mess
	entries.clear();
	splibPreamble.clear();
	libraryPeptides.clear();
	densityList.clear();
}

void Converter::processLibEntries()
{
	std::cout << "Start converting " << entries.size() << " entries to fpbin" << std::endl;
	int processed = 0;

	float totalIntensity, totalIntensityBelowMZ;
	double maxMz = 0, minMz = 999999;
	unsigned goodPeaks;

	for (auto entry = entries.begin(); entry != entries.end(); entry++)
	{
		if ((int) (std::distance(entries.begin(), entry) - processed) > (int) entries.size() / 10)
		{
			processed = std::distance(entries.begin(), entry);
			std::cout << (int) (((float) processed / entries.size()) * 100) << "%.." << std::flush;
		}

		totalIntensity = 0;
		totalIntensityBelowMZ = 0;
		goodPeaks = 0;

		for (auto & peak : entry->peaks)
		{
			if ((peak.mz > (entry->precurosrMz - 60) && peak.mz < (entry->precurosrMz + 20))
					|| peak.mz < parameters->filterLightIonsMzThreshold)
			{
				peak.intensity = 0;
			}

			if (parameters->peakScalingMzPower == 0.0 && parameters->peakScalingIntensityPower == 0.5)
			{
				peak.intensity = sqrt(peak.intensity); //Shortcut to sqrt Intensity
			}
			else
			{
				peak.intensity = pow(peak.mz, parameters->peakScalingMzPower)
						* pow(peak.intensity, parameters->peakScalingIntensityPower);
			}

			if (!peak.annotation.empty() && peak.annotation[0] == '?')
			{
				peak.intensity *= parameters->peakScalingUnassignedPeaks;
			}
		}

		/*
		 * Sort the peaks of library entry and constrain by the amount of library peaks used
		 */
		std::sort(entry->peaks.begin(), entry->peaks.end(), peaksCompareDesc);
		if (parameters->filterLibMaxPeaksUsed < entry->peaks.size())
		{
			entry->peaks.resize(parameters->filterLibMaxPeaksUsed);
		}
	}

	std::cout << "100%" << std::endl;
	std::cout << "Done filtering, " << entries.size() << " entries left" << std::endl << std::endl;
}

void Converter::binLibEntries()
{
	std::cout << "Start binning " << entries.size() << " entries" << std::endl;
	int intMZ = 0;

//Reset progress counter
	int processed = 0;
	for (auto entry = entries.begin(); entry != entries.end(); entry++)
	{
		if ((int) (std::distance(entries.begin(), entry) - processed) > (int) entries.size() / 10)
		{
			processed = std::distance(entries.begin(), entry);
			std::cout << (int) (((float) processed / entries.size()) * 100) << "%.." << std::flush;
		}

		//Now set the values of the peaks in the corresponding bins
		for (auto & peak : entry->peaks)
		{
			//Adjust the start so the vector start at 0 instead of MZ_START (default:10)
			intMZ = (int) peak.mz - MZ_START;

			//Check if the peak is in bounds
			if (intMZ < 0)
			{
				//Out of bound, set lowest peaks
				entry->libSpectra[0] += peak.intensity;
				entry->libSpectra[1] += peak.intensity * parameters->peakBinningFractionToNeighbor;
			}
			else if (intMZ > (int) (entry->libSpectra.size() - 2))
			{
				//Out of bound, set highest two peaks
				entry->libSpectra[entry->libSpectra.size() - 2] += peak.intensity
						* parameters->peakBinningFractionToNeighbor;
				entry->libSpectra[entry->libSpectra.size() - 1] += peak.intensity;
			}
			else
			{
				//Set the intensity
				entry->libSpectra[intMZ - 1] += peak.intensity * parameters->peakBinningFractionToNeighbor;
				entry->libSpectra[intMZ] += peak.intensity;
				entry->libSpectra[intMZ + 1] += peak.intensity * parameters->peakBinningFractionToNeighbor;
			}
		}

		/*
		 * In SpectraST this step happens after the dot product, but can be performed already
		 * slight change in the results will happen
		 */
		float squareMagnitude = 0, magnitude = 0;
		for (auto & peakIntensity : entry->libSpectra)
		{
			squareMagnitude += peakIntensity * peakIntensity;
		}

		magnitude = sqrt(squareMagnitude);

		if (magnitude > 0)
		{
			for (auto & peakIntensity : entry->libSpectra)
				peakIntensity /= magnitude;
		}
	}
	std::cout << "100%" << std::endl << "Done binning" << std::endl << std::endl;

}
} //namespace FastPassLibrary
