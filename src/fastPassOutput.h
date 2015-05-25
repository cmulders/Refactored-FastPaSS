/*
 * fastPassOutput.h
 *
 *  Created on: Jan 21, 2015
 *      Author: geforce
 */

#ifndef FASTPASSOUTPUT_H_
#define FASTPASSOUTPUT_H_

#include "fastPassParams.h"
#include "fastPassQuery.h"
#include "fastPassLibrary.h"
#include <iosfwd>
#include <string>

//File location of pepXML_std.xsl
#define PEPXML_STD_XSL "/usr/local/tpp/schema/"
#define PEPXML_NAMESPACE "http://regis-web.systemsbiology.net/pepXML"
#define PEPXML_SCHEMA "pepXML_v18.xsd"

#define PROTON_MASS 1.00739

namespace FastPassOutput {

class PepXML
{
		u_int64_t queryHitCount = 0;

		FastPassLibrary::Info* fastPassLibary;
		Parameters::PStruct *parameters;
		std::ofstream* outputStream;
		std::string originalQueryFile;

		void printStartQuery(FastPassQuery::QueryPeptide *query);
		void printEndQuery(FastPassQuery::QueryPeptide *query);
		void printQueryHit(FastPassQuery::QueryPeptide *query, FastPassQuery::Candidate *queryHit, unsigned int hitRank, FastPassLibrary::Library_peptide *libraryPeptide);

	public:
		std::string outputFilePath;
		std::string outputFileName;
		std::string outputFileExtension;
		std::string outputFullPathFileName;

		PepXML();
		~PepXML();

		void setParameters(Parameters::PStruct *parameters)
		{
			this->parameters = parameters;
		}
		void printHeader(FastPassQuery::QueryMetaData*);
		void printFooter();

		bool setOutputFile(std::string outputFile);

		void addQueries(std::vector<FastPassQuery::QueryPeptide> *queries);
		void addLibraryInfo(FastPassLibrary::Info *fastPassLibrary);


};

} //namespace FastPassOutput
#endif /* FASTPASSOUTPUT_H_ */
