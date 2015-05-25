/*
 * fastPassQuery.h
 *
 *  Created on: Jan 9, 2015
 *      Author: geforce
 */

#ifndef FASTPASSQUERY_H_
#define FASTPASSQUERY_H_

#include "fastPassParams.h"
#include "fastPassLibrary.h"
#include "gpuDevice.cuh"

#include <string>
#include <map>
#include <vector>
//#include <memory>
#include <fstream>
#include <utility>
#include <stdint.h>

namespace FastPassQuery {

enum InvalidQueryError
{
	none, notEnoughPeaks, intensityBelowMZ, lowMZrange, maxIntensityBelow

};

enum QueryType
{
	mgf, mzxml

};

struct QueryMetaData
{
	std::map<std::string,std::string> instrumentInfoMap;
};

struct Candidate
{
		int libId;
		float dotProd;
		float dotBias;
		float deltaDot;
		float fVal;
		unsigned int firstNonHomolog;
		Candidate(int libId, float dotProd, float dotBias) :
				libId(libId), dotProd(dotProd), dotBias(dotBias), deltaDot(0), fVal(0), firstNonHomolog(0)
		{
		}
};

struct QueryPeptide
{
		InvalidQueryError errorCode = none;
		int validSpectra = true;
		std::string title;
		int intPeptideMass;
		float retentionTime = -1;
		float peptideMass;
		int numberOfPeaks = 0;
		int fileLocation;
		int libraryMatchStart;
		int libraryMatchEnd;
		long int position;
		std::vector<std::pair<float, float> > fragments;
		std::vector<float> binnedSpectra;
		float totalIntensity;
		//Match; libid, {dot, dotbias, delta}
		std::vector<Candidate> matches;

		double hitsMean;
		double hitsDotStdev;

		int numHits = 0; //Number of hits (dotProd > 0.1)
};

struct Batch
{
		int batchNumber = 0;

		int numberOfQueries;
		std::vector<QueryPeptide>::iterator startQueryBatch;
		std::vector<QueryPeptide>::iterator endQueryBatch;

		int numberOfLibraryEntries;
};

class Query
{
		QueryMetaData *metaData;
		Parameters::PStruct *parameters;
		FastPassLibrary::Info *libraryInfo;
		std::ifstream queryFileStream;

		void calcFval(Candidate *candidate, int totalHits);
		void detectHomologs(QueryPeptide *query);

	public:
		QueryType queryType;

		int numberOfQueries;
		std::vector<Batch> batches;
		std::vector<QueryPeptide> queryPeptides;

		Query();
		~Query();

		bool load(std::string queryFile);
		bool loadMzXML(std::string queryFile);
		bool loadMgf(std::string queryFile);

		void setParameters(Parameters::PStruct *parameters)
		{
			this->parameters = parameters;
		}
		void setLibrary(FastPassLibrary::Info *libraryInfo)
		{
			this->libraryInfo = libraryInfo;
		}

		void sortOnPeptideMass();

		void filterQueries();
		void binSpectraQueries();

		void assignLibraryrange();

		void createBatches();

		void assignDotAndBiasFromDevice(Batch *currentBatch, float* dotMatrix, float* biasMatrix);
		void scoreQueries();

		QueryMetaData * getMetaData()
		{
			return this->metaData;
		}
};

} //namespace FastPassQuery
#endif /* FASTPASSQUERY_H_ */
