/*
 * fastPassOutput.cpp
 *
 *  Created on: Jan 21, 2015
 *      Author: geforce
 */

#include "fastPassOutput.h"
#include <fstream>
#include "time.h"
#include <iostream>
#include <iomanip>

namespace FastPassOutput {

PepXML::PepXML()
{

}

PepXML::~PepXML()
{
	if (this->outputStream->is_open())
	{
		this->outputStream->close();
	}
}

bool PepXML::setOutputFile(std::string queryFileName)
{
	//The query file with the data
	this->originalQueryFile = queryFileName;

	//Initialize the filename and path
	this->outputFilePath = "";
	this->outputFileName = queryFileName.substr(0, queryFileName.rfind(".")); //"output"; //queryFileName;
	this->outputFileExtension = ".fastpass." + Parameters::Output::OutMap[this->parameters->outputExtension];
	this->outputFullPathFileName = this->outputFilePath + this->outputFileName + this->outputFileExtension;

	//Open the file
	this->outputStream = new std::ofstream();
	this->outputStream->open(this->outputFullPathFileName, std::ofstream::out | std::ofstream::trunc);
	return true;
}

void PepXML::printHeader(FastPassQuery::QueryMetaData* queryMetaData)
{
	(*this->outputStream) << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n"
			<< "<?xml-stylesheet type=\"text/xsl\" href=\"" << PEPXML_STD_XSL << "pepXML_std.xsl\"?>" << "\n";

	time_t now;
	time(&now);
	char dateTime[sizeof "2015-01-21T07:07:09"];
	strftime(dateTime, sizeof dateTime, "%FT%T", gmtime(&now));

	(*this->outputStream) << "<msms_pipeline_analysis date=\"" << dateTime << "\" xmlns=\"" << PEPXML_NAMESPACE << "\" "
			<< "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"" << PEPXML_NAMESPACE
			<< " " << PEPXML_STD_XSL << PEPXML_SCHEMA << "\" summary_xml=\""
			<< this->outputFileName + this->outputFileExtension << "\">" << "\n";

	(*this->outputStream) << "<msms_run_summary base_name=\"" << this->outputFileName << "\" ";

	//TODO: intstrumetn info could be implemented when MZXML was read

	// if instrument info is available (which is only the case for .mzXML input), also print them
	if (queryMetaData->instrumentInfoMap.size() > 0)
	{
		std::map<std::string, std::string>::iterator it;
		if ((it = queryMetaData->instrumentInfoMap.find("msManufacturer")) != queryMetaData->instrumentInfoMap.end())
			(*this->outputStream) << "msManufacturer=\"" << it->second << "\" ";

		if ((it = queryMetaData->instrumentInfoMap.find("msModel")) != queryMetaData->instrumentInfoMap.end())
			(*this->outputStream) << "msModel=\"" << it->second << "\" ";

		if ((it = queryMetaData->instrumentInfoMap.find("msIonization")) != queryMetaData->instrumentInfoMap.end())
			(*this->outputStream) << "msIonization=\"" << it->second << "\" ";

		if ((it = queryMetaData->instrumentInfoMap.find("msMassAnalyzer")) != queryMetaData->instrumentInfoMap.end())
			(*this->outputStream) << "msMassAnalyzer=\"" << it->second << "\" ";

		if ((it = queryMetaData->instrumentInfoMap.find("msDetector")) != queryMetaData->instrumentInfoMap.end())
			(*this->outputStream) << "msDetector=\"" << it->second << "\" ";
	}

	std::string queryFileExtension = this->originalQueryFile.substr(this->originalQueryFile.find_last_of("."));
	(*this->outputStream) << "raw_data_type=\"" << queryFileExtension << "\" ";
	(*this->outputStream) << "raw_data=\"" << queryFileExtension << "\">" << "\n";

	//TODO: hardcoded enzym for now
	(*this->outputStream) << "<sample_enzyme name=\"trypsin\">" << "\n";
	(*this->outputStream) << "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>" << "\n";
	(*this->outputStream) << "</sample_enzyme>" << "\n";

	/*	if (m_searchParams.enzymeForPepXMLOutput.empty())
	 {
	 // default to trypsin if -s_ENZ option is not set.
	 // hard-coded sample_enzyme element for trypsin
	 (*m_fout) << "<sample_enzyme name=\"trypsin\">" << endl;
	 (*m_fout) << "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>" << endl;
	 (*m_fout) << "</sample_enzyme>" << endl;

	 }
	 else
	 {

	 #ifndef STANDALONE_LINUX
	 ProteolyticEnzymeFactory* enzFactory = new ProteolyticEnzymeFactory();
	 ProteolyticEnzyme* enz = enzFactory->getProteolyticEnzyme(m_searchParams.enzymeForPepXMLOutput.c_str());

	 if (enz)
	 {
	 enz->writePepXMLTags(*m_fout);
	 m_proteolyticEnzyme = enz;
	 }
	 delete (enzFactory);
	 #else
	 (*m_fout) << "<sample_enzyme name=\"trypsin\">" << endl;
	 (*m_fout) << "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>" << endl;
	 (*m_fout) << "</sample_enzyme>" << endl;
	 #endif

	 }
	 */

	// search_summary element (not quite applicable to SpectraST), but have to be hard-coded
	(*this->outputStream) << "<search_summary base_name=\"" << this->outputFileName << "\" ";
	(*this->outputStream) << "search_engine=\"" << "SpectraST" << "\" ";
	(*this->outputStream) << "search_engine_real=\"" << "FastPaSS" << "\" ";
	(*this->outputStream) << "precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" ";
	(*this->outputStream) << "out_data_type=\"out\" out_data=\".tgz\" search_id=\"1\">" << "\n";

	/* TODO: implement sequence database reference
	 if (!this->parameters->databaseFile.empty())
	 {
	 (*this->outputStream) << "<search_database local_path=\"" << this->parameters->databaseFile << "\" type=\""
	 << this->parameters->databaseType << "\"/>" << std::endl;
	 }
	 //*/

	//Add FastPass parameters to this output
	printPepXMLSearchParams(this->parameters, this->outputStream);

	/*
	 string icatType("");
	 if (m_searchParams.expectedCysteineMod == "ICAT_cl") {
	 icatType = "cl";
	 } else if (m_searchParams.expectedCysteineMod == "ICAT_uc") {
	 icatType = "uc";
	 }

	 (*m_fout) << "<parameter name=\"icat_type\" value=\"" << icatType << "\"/>" << endl;
	 */

	// etc
	(*this->outputStream) << "</search_summary>" << "\n";
}

void PepXML::addLibraryInfo(FastPassLibrary::Info *fastPassLibrary)
{
	this->fastPassLibary = fastPassLibrary;
}

/*
 * PepXML::printStartQuery print the opening tags and query information
 */
void PepXML::printStartQuery(FastPassQuery::QueryPeptide *query)
{

	/*
	 *  parse the query name to find the startScan, endScan etc.
	 *  the query name is of the format: <baseName>.<startScan>.<endScan>.<charge>
	 *  										dot1Pos		dot2Pos	  dot3Pos
	 *  note that in a typical SpectraST search, the charge is not specified, so it will be set to zero
	 */
	std::string::size_type dot1Pos, dot2Pos, dot3Pos;
	std::string startScan("0"), endScan("0");
	int assumedCharge = 0;

	dot3Pos = query->title.rfind('.');
	dot2Pos = query->title.rfind('.', dot3Pos - 1);
	dot1Pos = query->title.rfind('.', dot2Pos - 1);

	//Check if we have a compatible format
	if (dot1Pos != std::string::npos && dot2Pos != std::string::npos && dot3Pos != std::string::npos)
	{
		startScan = query->title.substr(dot1Pos + 1, dot2Pos - dot1Pos - 1);
		endScan = query->title.substr(dot2Pos + 1, dot3Pos - dot2Pos - 1);
	}
	else
	{
		//Construct a fake tpp_compatible spectrum
		startScan.assign(std::to_string(query->fileLocation).substr(0, 4));
		endScan.assign(std::to_string(query->fileLocation).substr(0, 4));
		query->title.assign(query->title.substr(0, query->title.find("/")) + "." + startScan + "." + endScan + ".0");
		dot3Pos = query->title.length() - 2;
	}

	//Get the charge from the top match
	std::string topMatchPeptideName(this->fastPassLibary->getLibraryPeptide(query->matches[0].libId)->FullName);
	assumedCharge = atoi(topMatchPeptideName.substr(topMatchPeptideName.rfind('/') + 1).c_str());

	(*this->outputStream) << "<spectrum_query spectrum=\"";
	(*this->outputStream) << query->title.substr(0, dot3Pos + 1);
	(*this->outputStream) << assumedCharge << "\" start_scan=\"" << startScan << "\" end_scan=\"" << endScan << "\" ";

	double precursorNeutralMass = assumedCharge * (query->peptideMass - PROTON_MASS);

	(*this->outputStream) << "precursor_neutral_mass=\"" << std::setprecision(4) << std::fixed << precursorNeutralMass
			<< "\" ";
	(*this->outputStream) << "assumed_charge=\"" << assumedCharge << "\" ";
	(*this->outputStream) << "index=\"" << ++this->queryHitCount << "\"";

	if (query->retentionTime >= 0.0)
	{
		this->outputStream->precision(2);
		(*this->outputStream) << " retention_time_sec=\"" << std::fixed << query->retentionTime << "\"";
	}

	(*this->outputStream) << ">" << "\n";
	(*this->outputStream) << "<search_result>" << "\n";

}

/*
 * PepXML::::printEndQuery print the ending tags
 */
void PepXML::printEndQuery(FastPassQuery::QueryPeptide *query)
{
	(*this->outputStream) << "</search_result>" << "\n";
	(*this->outputStream) << "</spectrum_query>" << "\n";
}

/*
 * PepXML::::printQueryHit print a search hit
 */
//Macro for easy search score addition
#define SEARCH_SCORE(PRECISION,KEY,VALUE) { (*this->outputStream).precision(PRECISION); \
	(*this->outputStream) << "<search_score name=\"" << KEY << "\" value=\"" << VALUE << "\"/>" << "\n"; }

void PepXML::printQueryHit(
		FastPassQuery::QueryPeptide *query, FastPassQuery::Candidate *queryHit, unsigned int hitRank,
		FastPassLibrary::Library_peptide *libraryPeptide)
{

	std::string fullPeptide(libraryPeptide->FullName);
	//Striped peptide from X.AAA.X/C -> AAA
	std::string peptide(
			fullPeptide.substr(fullPeptide.find('.') + 1, fullPeptide.rfind('.') - fullPeptide.find('.') - 1));

	int assumedCharge = atof(fullPeptide.substr(fullPeptide.rfind("/") + 1).c_str());

	double libraryPeptideMass = assumedCharge * (libraryPeptide->PrecursorMZ - PROTON_MASS);
	double precursorNeutralMass = assumedCharge * (query->peptideMass - PROTON_MASS);

	//Strip any modifactions from the peptide
	size_t currentPos = 0;
	while ((currentPos = peptide.find('[')) != std::string::npos)
	{
		//N or C terminal modification, delete the marker also
		if (currentPos > 0 && (peptide.at(currentPos - 1) == 'n' || peptide.at(currentPos - 1) == 'c'))
		{
			currentPos -= 1;
		}
		peptide.erase(currentPos, peptide.find(']', currentPos) - currentPos + 1);
	}

	//Calculates the number of missed cleavages and number of tryptic terminals
	unsigned int ntt = 0, nmc = 0;
	for (std::string::size_type i = 0; i < peptide.length() - 1; i++)
	{
		if ((peptide[i] == 'K' || peptide[i] == 'R') && (peptide[i + 1] != 'P'))
		{
			nmc++;
		}
	}
	std::string prevAA(fullPeptide.substr(0, 1)), nextAA(fullPeptide.substr(fullPeptide.rfind('.') + 1, 1));
	if (prevAA == "-" || prevAA == "[" || ((prevAA == "K" || prevAA == "R") && peptide.substr(0, 1) != "P"))
	{
		ntt++;
	}
	if (nextAA == "-" || nextAA == "]"
			|| ((peptide.substr(peptide.size() - 1, 1) == "K" || peptide.substr(peptide.size() - 1, 1) == "R")
					&& nextAA != "P"))
	{
		ntt++;
	}

	std::string stripedPeptide(peptide);

	//Parse out the proteins
	std::string fullprotein(libraryPeptide->protein);
	std::vector<std::string> proteins;
	std::string::size_type slashPos = fullprotein.find('/');
	if (slashPos == std::string::npos || slashPos >= fullprotein.length() - 1)
	{
		// single protein, old format
		proteins.push_back(fullprotein);
	}
	else
	{

		std::string::size_type lastSlash = slashPos;
		while ((slashPos = fullprotein.find('/', slashPos + 1)) != std::string::npos)
		{
			proteins.push_back(fullprotein.substr(lastSlash + 1, (slashPos - lastSlash) - 1));
			lastSlash = slashPos;
		}
		//Push on the last part
		proteins.push_back(fullprotein.substr(lastSlash + 1, std::string::npos));
	}

	(*this->outputStream) << "<search_hit hit_rank=\"" << hitRank << "\" peptide=\"" << peptide << "\" ";
	(*this->outputStream) << "peptide_prev_aa=\"" << prevAA << "\" peptide_next_aa=\"" << nextAA << "\" ";
	(*this->outputStream) << "protein=\"" << proteins[0] << "\" num_tot_proteins=\"" << proteins.size() << "\" ";
	(*this->outputStream) << "calc_neutral_pep_mass=\"";
	(*this->outputStream).precision(4);
	(*this->outputStream) << std::fixed << libraryPeptideMass;
	(*this->outputStream) << "\" massdiff=\"";
	(*this->outputStream).precision(4);
	(*this->outputStream) << std::fixed << std::showpoint << std::showpos << precursorNeutralMass - libraryPeptideMass;
	(*this->outputStream) << std::noshowpoint << std::noshowpos;
	(*this->outputStream) << "\" num_tol_term=\"" << ntt;
	(*this->outputStream) << "\" num_missed_cleavages=\"" << nmc << "\">" << "\n";

	if (proteins.size() > 1)
	{
		for (std::vector<std::string>::size_type index = 1; index < proteins.size(); index++)
		{
			(*this->outputStream) << "<alternative_protein protein=\"" << proteins[index] << "\"/>" << "\n";
		}
	}
	/*

	 if (!p->mods.empty() || !p->nTermMod.empty() || !p->cTermMod.empty()) {
	 // peptide is modified

	 (*m_fout) << "<modification_info ";
	 if (!p->nTermMod.empty()) {
	 if (m_searchParams.precursorMzUseAverage) {
	 (*m_fout) << "mod_nterm_mass=\"" << Peptide::getAAPlusModAverageMass('n', p->nTermMod) << "\" ";
	 } else {
	 (*m_fout) << "mod_nterm_mass=\"" << Peptide::getAAPlusModMonoisotopicMass('n', p->nTermMod) << "\" ";
	 }
	 }
	 if (!p->cTermMod.empty()) {
	 if (m_searchParams.precursorMzUseAverage) {
	 (*m_fout) << "mod_cterm_mass=\"" << Peptide::getAAPlusModAverageMass('c', p->cTermMod) << "\" ";
	 } else {
	 (*m_fout) << "mod_cterm_mass=\"" << Peptide::getAAPlusModMonoisotopicMass('c', p->cTermMod) << "\" ";
	 }
	 }

	 (*m_fout) << "modified_peptide=\"" << p->interactStyle() << "\">" << endl;

	 map<int, string>::iterator i;
	 for (i = p->mods.begin(); i != p->mods.end(); i++) {
	 // plus 1 because we want the mod position to be one-based; in the Peptide class it is stored as zero-based
	 (*m_fout) << "<mod_aminoacid_mass position=\"" << (*i).first + 1;

	 if (m_searchParams.precursorMzUseAverage) {
	 (*m_fout) << "\" mass=\"" << Peptide::getAAPlusModAverageMass(p->stripped[(*i).first], (*i).second) << "\"/>" << endl;
	 } else {
	 (*m_fout) << "\" mass=\"" << Peptide::getAAPlusModMonoisotopicMass(p->stripped[(*i).first], (*i).second) << "\"/>" << endl;
	 }
	 }
	 (*m_fout) << "</modification_info>" << endl;
	 }

	 }
	 */

	//SEARCH_SCORE(3, "LibID", queryHit->libId)
	//Output the search score
	SEARCH_SCORE(3, "dot", queryHit->dotProd)
	SEARCH_SCORE(3, "delta", queryHit->deltaDot)
	SEARCH_SCORE(3, "dot_bias", queryHit->dotBias)
	SEARCH_SCORE(3, "precursor_mz_diff", (precursorNeutralMass - libraryPeptideMass) / assumedCharge)
	SEARCH_SCORE(0, "hits_num", query->numHits)
	SEARCH_SCORE(3, "hits_mean", query->hitsMean)
	SEARCH_SCORE(3, "hits_stdev", query->hitsDotStdev)
	SEARCH_SCORE(3, "fval", queryHit->fVal)

	SEARCH_SCORE(3, "p_value", std::scientific << -1.0)
	//TODO: hard-coded
	SEARCH_SCORE(3, "KS_score", std::fixed << 0.0)
	//TODO: hard-coded

	// internal is zero-indexing, for output requires 1 indexing
	SEARCH_SCORE(0, "first_non_homolog", queryHit->firstNonHomolog)

	SEARCH_SCORE(3, "open_mod_mass", std::fixed << 0.0)
	//TODO: hard-coded
	SEARCH_SCORE(0, "open_mod_locations", "")
	//TODO: hard-coded

	SEARCH_SCORE(0, "charge", assumedCharge)
	SEARCH_SCORE(0, "lib_file_offset", libraryPeptide->splibPos)
	SEARCH_SCORE(4, "lib_probability", libraryPeptide->libProbability)
	SEARCH_SCORE(4, "lib_status", libraryPeptide->libStatus)
	SEARCH_SCORE(0, "lib_num_replicates", libraryPeptide->libReplicates)
	SEARCH_SCORE(4, "lib_remark", libraryPeptide->libRemark)

	(*this->outputStream) << "</search_hit>" << "\n";

}

void PepXML::addQueries(std::vector<FastPassQuery::QueryPeptide> *queries)
{
	FastPassQuery::Candidate *queryHit;
	for (auto query : *queries)
	{
/*		std::cout << query.title << std::endl;
		for (auto & match : query.matches)
		{
			std::cout << this->fastPassLibary->getLibraryPeptide(match.libId)->FullName
					<< "	" << match.dotProd
					<< "	" << match.dotBias
					<< "	" << match.dotProd
					<< "	" << match.fVal << std::endl;
		}
		std::cout << std::endl << std::endl << std::endl;
*/
		//If there are no matches above the threshold, skip this query
		if (query.matches.size() == 0 || query.matches.at(0).fVal < this->parameters->hitListTopHitFvalThreshold)
		{
			continue;
		}

		printStartQuery(&query);

		for (unsigned int hitRank = 0; hitRank < this->parameters->hitListShowMaxRank && hitRank < query.matches.size();
				hitRank++)
		{
			queryHit = &(query.matches.at(hitRank));

			//Filter out the low ranked hits
			if (hitRank == 0 && queryHit->fVal < this->parameters->hitListTopHitFvalThreshold) continue;
			if (hitRank > 0 && queryHit->fVal < this->parameters->hitListLowerHitsFvalThreshold) continue;

			printQueryHit(&query, queryHit, hitRank + 1, this->fastPassLibary->getLibraryPeptide(queryHit->libId));
		}

		printEndQuery(&query);
	}
}

void PepXML::printFooter()
{
	(*this->outputStream) << "</msms_run_summary>" << "\n";
	(*this->outputStream) << "</msms_pipeline_analysis>" << "\n";
}

} //namespace FastPassOutput

