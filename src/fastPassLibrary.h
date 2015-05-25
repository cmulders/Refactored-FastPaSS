/*
 * fastPassLibrary.h
 *
 *  Created on: Jan 9, 2015
 *      Author: geforce
 */

#ifndef FASTPASSLIBRARY_H_
#define FASTPASSLIBRARY_H_

#include "fastPassParams.h"

#include <string>
#include <vector>
#include <fstream>

namespace FastPassLibrary {

struct Density
{
		uint32_t numberOfEntries;
		uint32_t firstLocation;
		uint32_t lastLocation;
		Density()
		{
		}
		Density(int numberOfEntries, int firstLocation, int lastLocation) :
				numberOfEntries(numberOfEntries), firstLocation(firstLocation), lastLocation(lastLocation)
		{
		}
};

struct Library_peptide
{
		char FullName[256];
		char StripedName[256];
		char protein[256];
		char libStatus[100];
		char libRemark[100];
		int LibID;
		float PrecursorMZ;
		int NumPeaks;
		long int splibPos;
		float libProbability;
		unsigned int libReplicates;
};

class Info
{

		std::ifstream libraryFileStream;
	public:
		std::string libraryFile;
		uint32_t numberOfEntries, precursorRange;
		std::vector<Density> densityList; //indexed with precursor mass as index

		std::vector<Library_peptide> libraryPeptides; //indexed with libid index

		Info(); //Constructor
		~Info(); //Destructor

		/*
		 * Read the library into the internal variables and read the density struct into a vector
		 */
		bool readFPLibrary(std::string libraryFile);
		void printLibraryInfo();

		/*
		 * Retrieves density information for a specific mass
		 */
		Density *getDensityForPrecursorMass(int mass);

		/*
		 * Retrieves library peptide information for a specific libid
		 */
		Library_peptide *getLibraryPeptide(int libIndex);
		/*
		 * Gets the peaks of the library entries for a specific range
		 */
		bool readPeaksForLibraryRangeIntoPointer(int libStart, int libEnd, float *pointer);
};



struct Peak
{
		double mz;
		double intensity;
		std::string annotation;
		std::string info;
};

struct SPLibEntry
{
		long int offset;
		bool valid = true;
		int libId;
		std::string fullName;
		double precurosrMz;
		std::string status;
		unsigned int numPeaks;
		std::vector<Peak> peaks;
		std::map<std::string,std::string> comments;
		std::vector<float> libSpectra;
};

class Converter
{
		std::vector<std::string> splibPreamble; // the preamble of this lib file
		std::vector<SPLibEntry> entries; // all the entries
		Parameters::PStruct *parameters;

		void processLibEntries();
		void binLibEntries();
		void writeFplib();

		//Output structs for .fpbin
		std::vector<Library_peptide> libraryPeptides; //indexed with libid index
		std::vector<Density> densityList; //indexed with precursor mass as index

	public:
		Converter();
		~Converter();

		void convert(std::string splib);

		void setParameters(Parameters::PStruct *parameters)
		{
			this->parameters = parameters;
		}
};

} //namespace FastPassLibrary

#endif /* FASTPASSLIBRARY_H_ */
