/*
 * fastPassParams.cpp
 *
 *  Created on: Jan 6, 2015
 */
#include "fastPassParams.h"

#include <cstdlib>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <algorithm>

/*
 * Initialize the mapping for possible input to the specific Parameter in the struct.
 * This way new ways of initializing a parameter can be set here.
 */
void initializeEnumMapping()
{
	using namespace Parameters;

	PMap["queryFile"] = queryFiles;
	PMap["-Q"] = queryFiles;

	/*
	 * CONVERTING OPTIONS
	 */
	PMap["convertSpLib"] = convertSpLib;
	PMap["-C"] = convertSpLib;

	/*
	 * TESTING OPTIONS
	 */
	PMap["testFvalConfidence"] = testFvalConfidence;
	PMap["-TFvalConfi"] = testFvalConfidence;
	PMap["testVariable"] = testVariable;
	PMap["-TVariable"] = testVariable;
	PMap["testMinValueVariable"] = testMinValueVariable;
	PMap["-TMin"] = testMinValueVariable;
	PMap["testMaxValueVariable"] = testMaxValueVariable;
	PMap["-TMax"] = testMaxValueVariable;
	PMap["testValueVariableStepSize"] = testValueVariableStepSize;
	PMap["-TStep"] = testValueVariableStepSize;

	/*
	 * GENERAL OPTIONS
	 */
	PMap["flagShowHelp"] = flagShowHelp;
	PMap["-H"] = flagShowHelp;

	PMap["showParameters"] = showParameters;
	PMap["-DispP"] = showParameters;
	PMap["normalOutput"] = normalOutput;

	PMap["parameterFile"] = parameterFile;
	PMap["-F"] = parameterFile;

	PMap["libraryFile"] = libraryFile;
	PMap["-L"] = libraryFile;

	PMap["maxGPUMemoryUsage"] = maxGPUMemoryUsage;
	PMap["-GPUMem"] = maxGPUMemoryUsage;

	/*
	 * CANDIDATE SELECTION AND SCORING OPTIONS
	 */
	PMap["precursorMzTolerance"] = precursorMzTolerance;
	PMap["-sM"] = precursorMzTolerance;

	PMap["detectHomologs"] = detectHomologs;
	PMap["fvalFractionDelta"] = fvalFractionDelta;
	PMap["fvalUseDotBias"] = fvalUseDotBias;

	/*
	 * OUTPUT AND DISPLAY OPTIONS
	 */
	PMap["hitListTopHitFvalThreshold"] = hitListTopHitFvalThreshold;
	PMap["hitListLowerHitsFvalThreshold"] = hitListLowerHitsFvalThreshold;
	PMap["hitListShowMaxRank"] = hitListShowMaxRank;
	PMap["outputExtension"] = outputExtension;
	PMap["-OExt"] = outputExtension;

	/*
	 * SPECTRUM FILTERING OPTIONS
	 */
	PMap["filterMinPeakCount"] = filterMinPeakCount;
	PMap["filterAllPeaksBelowMz"] = filterAllPeaksBelowMz;
	PMap["filterMaxIntensityBelow"] = filterMaxIntensityBelow;
	PMap["filterMinMzRange"] = filterMinMzRange;
	PMap["filterCountPeakIntensityThreshold"] = filterCountPeakIntensityThreshold;

	/*
	 * SPECTRUM PROCESSING OPTIONS
	 */
	PMap["filterRemovePeakIntensityThreshold"] = filterRemovePeakIntensityThreshold;
	PMap["filterMaxPeaksUsed"] = filterMaxPeaksUsed;
	PMap["filterMaxDynamicRange"] = filterMaxDynamicRange;
	PMap["peakScalingMzPower"] = peakScalingMzPower;
	PMap["peakScalingIntensityPower"] = peakScalingIntensityPower;
	PMap["peakScalingUnassignedPeaks"] = peakScalingUnassignedPeaks;
	PMap["peakBinningFractionToNeighbor"] = peakBinningFractionToNeighbor;
	PMap["filterLibMaxPeaksUsed"] = filterLibMaxPeaksUsed;
	PMap["filterLightIonsMzThreshold"] = filterLightIonsMzThreshold;

}

std::string getParameterSyntax(std::map<int, std::vector<std::string> > Options, int Parameter)
{
	//Stream to append our syntax to
	std::stringstream syntax;

	//First option in the list flag
	bool firstOption = true;

	//Go thorugh all the syntaxes of this argument
	for (std::vector<std::string>::iterator optionSyntax = Options.find(Parameter)->second.begin();
			optionSyntax != Options.find(Parameter)->second.end(); ++optionSyntax)
	{
		syntax << (firstOption ? "" : ", ") << (*optionSyntax).c_str();
		if (firstOption) firstOption = false;
	}
	return syntax.str();
}

inline void printParameterHelp(std::map<int, std::vector<std::string> > Options, int Parameter, std::string information)
{
	std::cout << "\t" << getParameterSyntax(Options, Parameter) << "\n\t\t" << information << std::endl;
}

/*
 * Shows command line help when options for the possible options
 */
void showHelpParameters()
{
	using namespace Parameters;
	std::cout << "\nFastPaSS Help\n\n";

	//General introduction
	std::cout << "FastPaSS is a spectral library search program that is based on SpectraST 4 but uses also the GPU. ";
	std::cout << "It uses a matrix multiplication to calculate the dot product of a query spectrum ";
	std::cout << "and library spectrum. ";
	std::cout << "This matrix multiplication is done on the GPU and is much faster than the ";
	std::cout << "original SpectraST calculations on the CPU. ";
	std::cout << std::endl;
	std::cout << "USAGE: FastPaSS [Options] queryFile..." << std::endl;

	std::map<int, std::vector<std::string> > Options;

	//Collect all the options and ways to write the options
	for (std::map<std::string, PEnum>::iterator PMapIterator = PMap.begin(); PMapIterator != PMap.end(); PMapIterator++)
	{
		if (PMapIterator->second == unknown) continue; //Skip the unknown dummy Enum type
		Options[PMapIterator->second].push_back(PMapIterator->first);
	}

	std::cout << std::endl;
	std::cout
			<< "Options starting with a dash '-' are command line options, the other could be put in a file 'parameter = value'";
	std::cout << std::endl << std::endl;

	std::cout << "CONVERTING OPTIONS" << std::endl;
	{
		printParameterHelp(Options, convertSpLib,
				"Converting of a .splib to .fpbin (fastpass compatible), use this option and add the splib's as queryfiles with -Q");
	}

	std::cout << "TESTING OPTIONS" << std::endl;
	{
		printParameterHelp(Options, testFvalConfidence,
				"Cutoff confidence interval used for testing purposes (default: 0.85)");
		printParameterHelp(Options, testVariable, "The variable which we are looping over (default: none)");
		printParameterHelp(Options, testMinValueVariable,
				"Minimum value for the variable were are looping (default: 0)");
		printParameterHelp(Options, testMaxValueVariable,
				"Maximum value for the variable were are looping  (default: 1)");
		printParameterHelp(Options, testValueVariableStepSize, "Step size for the loop (default: 1)");
	}

	std::cout << "GENERAL OPTIONS" << std::endl;
	{
		printParameterHelp(Options, queryFiles,
				"Query files can be added with these commands or added to the end of the command line");
		printParameterHelp(Options, flagShowHelp, "Show this help information.");
		printParameterHelp(Options, showParameters,
				"Show the used parameters. This output could be directly put into a parameter file.");
		printParameterHelp(Options, libraryFile,
				"The library that the program will use, must be already converted for FastPaSS.");
		printParameterHelp(Options, parameterFile,
				"The parameter file that the program will use to read it's parameters.");
		printParameterHelp(Options, maxGPUMemoryUsage, "The maximum amount of GPU memory usage.");

	}

	std::cout << std::endl << "CANDIDATE SELECTION AND SCORING OPTIONS" << std::endl;
	{
		printParameterHelp(Options, precursorMzTolerance,
				"Specify precursor m/z tolerance in Th. Monoisotopic mass is assumed, must be greater than 0.0. (default: 3.0)");
		printParameterHelp(Options, detectHomologs,
				"Detect homologous lower hits up to <rank>. Looks for lower hits homologous to the top hit and adjust delta accordingly. (Default: is 4)");
		printParameterHelp(Options, fvalFractionDelta,
				"Specify the fraction of the normalized delta score (delta/dot) in the F-value formula. (default: 0.4)");
		printParameterHelp(Options, fvalUseDotBias,
				"Use dot bias to penalize high-scoring matches with massive noise and/or dominant peak. (default: true)");
	}

	std::cout << std::endl << "OUTPUT AND DISPLAY OPTIONS" << std::endl;
	{
		printParameterHelp(Options, hitListTopHitFvalThreshold,
				"Minimum F value threshold for the top hit. Only top hits having F value greater than <thres> will be printed.	(Default = 0.03)");
		printParameterHelp(Options, hitListLowerHitsFvalThreshold,
				"Minimum F value threshold for the lower hits. Only lower hits having F value greater than <thres> will be printed. (Default = 0.45)");
		printParameterHelp(Options, hitListShowMaxRank,
				"Maximum rank for hits shown for each query. (default: 1, only show top rank)");
		printParameterHelp(Options, outputExtension,
				"Output format. The search result will be written to a file with the same base name as the search file, with extension <ext>.");
	}

	std::cout << std::endl << "SPECTRUM FILTERING OPTIONS" << std::endl;
	{
		printParameterHelp(Options, filterMinPeakCount,
				"Discard query spectra with fewer than <thres> peaks above threshold. (default: 10)");
		printParameterHelp(Options, filterAllPeaksBelowMz,
				"Discard query spectra with almost no peaks above a certain m/z value. All query spectra with 95%+ of the total intensity below <m/z> will be removed. (default: 520)");
		printParameterHelp(Options, filterMaxIntensityBelow,
				"Discard query spectra with no peaks with intensity above <inten>.	(Default is 0)");
		printParameterHelp(Options, filterMinMzRange,
				"Discard query spectra with m/z range narrower than <range>.	(Default is 350)");
		printParameterHelp(Options, filterCountPeakIntensityThreshold,
				"Minimum peak intensity for peaks to be counted. Only peaks with intensity above <thres> will be counted to meet the requirement for minimum number of peaks. (default: 2.01)");
	}

	std::cout << std::endl << "SPECTRUM PROCESSING OPTIONS" << std::endl;
	{
		printParameterHelp(Options, filterRemovePeakIntensityThreshold,
				"Noise peak threshold. All peaks with intensities below <thres> will be zeroed. (default: 2.01");
		printParameterHelp(Options, filterMaxPeaksUsed,
				"Remove all but the top <num> peaks in query spectra. (default: 150)");
		printParameterHelp(Options, filterMaxDynamicRange,
				"Remove all peaks smaller than 1/<num> of the base (highest) peak in query spectra. (default: 1000)");
		printParameterHelp(Options, peakScalingMzPower,
				"Intensity scaling power with respect to the m/z value and the raw intensity. The scaled intensity will be (m/z)^<mzpow> * (raw intensity)^<intpow> (default: 0.0)");
		printParameterHelp(Options, peakScalingIntensityPower,
				"Intensity scaling power with respect to the m/z value and the raw intensity. The scaled intensity will be (m/z)^<mzpow> * (raw intensity)^<intpow> (default: 0.5)");
		printParameterHelp(Options, peakBinningFractionToNeighbor,
				"Fraction of the scaled intensity assigned to neighboring bins. (default: 0.5)");
		printParameterHelp(Options, peakScalingUnassignedPeaks,
				"Scaling factor for unassigned peaks in library spectra. Unassigned peaks in the library spectra will be scaled by <factor>. (Default is 1.0)");
		printParameterHelp(Options, filterLibMaxPeaksUsed,
				"Remove all but the top <num> peaks in the LIBRARY spectra.	(Default is 50)");
		printParameterHelp(Options, filterLightIonsMzThreshold,
				"Remove all light ions with m/z lower than <thres> Th for both (library,not yet) and query spectra. (default: 180)");

	}
	std::cout << std::endl << std::endl;
}

void setDefaults(Parameters::PStruct *parameters)
{
	assert(parameters != NULL);

	parameters->queryFiles.reserve(1);

	/*
	 * CONVERTING OPTIONS
	 */
	setParameter(parameters, "convertSpLib", "false");

	/*
	 * TESTING OPTIONS
	 */
	setParameter(parameters, "testFvalConfidence", "0.85");
	setParameter(parameters, "testVariable", "none");
	setParameter(parameters, "testMinValueVariable", "0");
	setParameter(parameters, "testMaxValueVariable", "1");
	setParameter(parameters, "testValueVariableStepSize", "1");

	/*
	 * GENERAL OPTIONS
	 */
	setParameter(parameters, "flagShowHelp", "false");
	setParameter(parameters, "showParameters", "false");
	setParameter(parameters, "normalOutput", "true");
	setParameter(parameters, "parameterFile", "fastpass.params");
	setParameter(parameters, "maxGPUMemoryUsage", std::to_string(1024 * 1024 * 100)); //100MB

	/*
	 * CANDIDATE SELECTION AND SCORING OPTIONS
	 */
	setParameter(parameters, "precursorMzTolerance", "3.0");
	setParameter(parameters, "detectHomologs", "4");
	setParameter(parameters, "fvalFractionDelta", "0.4");
	setParameter(parameters, "fvalUseDotBias", "true");

	/*
	 * OUTPUT AND DISPLAY OPTIONS
	 */
	setParameter(parameters, "hitListTopHitFvalThreshold", "0.03");
	setParameter(parameters, "hitListLowerHitsFvalThreshold", "0.45");
	setParameter(parameters, "hitListShowMaxRank", "1");
	setParameter(parameters, "outputExtension", "xml");

	/*
	 * SPECTRUM FILTERING OPTIONS
	 */
	setParameter(parameters, "filterMinPeakCount", "10");
	setParameter(parameters, "filterAllPeaksBelowMz", "520");
	setParameter(parameters, "filterMaxIntensityBelow", "0");
	setParameter(parameters, "filterMinMzRange", "350");
	setParameter(parameters, "filterCountPeakIntensityThreshold", "2.01");

	/*
	 * SPECTRUM PROCESSING OPTIONS
	 */
	setParameter(parameters, "filterRemovePeakIntensityThreshold", "2.01");
	setParameter(parameters, "filterMaxPeaksUsed", "150");
	setParameter(parameters, "filterMaxDynamicRange", "1000");
	setParameter(parameters, "peakScalingMzPower", "0.0");
	setParameter(parameters, "peakScalingIntensityPower", "0.5");
	setParameter(parameters, "peakScalingUnassignedPeaks", "1.0");
	setParameter(parameters, "peakBinningFractionToNeighbor", "0.5");
	setParameter(parameters, "filterLibMaxPeaksUsed", "50");
	setParameter(parameters, "filterLightIonsMzThreshold", "180");

}

/*
 * Prints the parameter struct to stdout with printf
 * Uses local blocks to organize the code
 */
void printParameters(Parameters::PStruct *parameters)
{
	assert(parameters != NULL);

	std::cout << "----------------PARAMETERS---------------" << std::endl;

	std::cout << "# QUERY FILES" << std::endl;
	for (std::vector<std::string>::iterator itQueryFiles = parameters->queryFiles.begin();
			itQueryFiles != parameters->queryFiles.end(); ++itQueryFiles)
	{
		std::cout << "queryFile = " << (*itQueryFiles) << std::endl;
	}
	std::cout << std::endl;

	std::cout << "# TESTING OPTIONS" << std::endl;
	{
		std::cout << "testFvalConfidence = " << std::noboolalpha << parameters->testFvalConfidence << std::endl;
		std::cout << "testVariable = " << parameters->testVariable << std::endl;
		std::cout << "testMinValueVariable = " << parameters->testMinValueVariable << std::endl;
		std::cout << "testMaxValueVariable = " << parameters->testMaxValueVariable << std::endl;
		std::cout << "testValueVariableStepSize = " << parameters->testValueVariableStepSize << std::endl;
	}
	std::cout << std::endl;

	std::cout << "# GENERAL OPTIONS" << std::endl;
	{
		std::cout << "flagShowHelp = " << std::noboolalpha << parameters->flagShowHelp << std::endl;
		std::cout << "showParameters = " << std::noboolalpha << parameters->showParameters << std::endl;
		std::cout << "normalOutput = " << std::noboolalpha << parameters->normalOutput << std::endl;
		std::cout << "libraryFile = " << parameters->libraryFile << std::endl;
		std::cout << "parameterFile = " << parameters->parameterFile << std::endl;
		std::cout << "maxGPUMemoryUsage = " << parameters->maxGPUMemoryUsage << " ("
				<< parameters->maxGPUMemoryUsage / 1024 / 1024 << "MB)" << std::endl;
		std::cout << "useSp4Scoring = true" << std::endl;
	}
	std::cout << std::endl;

	std::cout << "# CANDIDATE SELECTION AND SCORING OPTIONS" << std::endl;
	{
		std::cout << "precursorMzTolerance = " << parameters->precursorMzTolerance << std::endl;
		std::cout << "detectHomologs = " << parameters->detectHomologs << std::endl;
		std::cout << "fvalFractionDelta = " << parameters->fvalFractionDelta << std::endl;
		std::cout << "fvalUseDotBias = " << std::noboolalpha << parameters->fvalUseDotBias << std::endl;
	}
	std::cout << "\n";

	std::cout << "# OUTPUT AND DISPLAY OPTIONS" << std::endl;
	{
		std::cout << "hitListTopHitFvalThreshold = " << parameters->hitListTopHitFvalThreshold << std::endl;
		std::cout << "hitListLowerHitsFvalThreshold = " << parameters->hitListLowerHitsFvalThreshold << std::endl;
		std::cout << "hitListShowMaxRank = " << parameters->hitListShowMaxRank << std::endl;
		std::cout << "outputExtension = " << Parameters::Output::OutMap.at(parameters->outputExtension) << std::endl;
	}
	std::cout << std::endl;

	std::cout << "# SPECTRUM FILTERING OPTIONS" << std::endl;
	{
		std::cout << "filterMinPeakCount = " << parameters->filterMinPeakCount << "\n";
		std::cout << "filterAllPeaksBelowMz = " << parameters->filterAllPeaksBelowMz << "\n";
		std::cout << "filterMaxIntensityBelow = " << parameters->filterMaxIntensityBelow << "\n";
		std::cout << "filterMinMzRange = " << parameters->filterMinMzRange << "\n";
		std::cout << "filterCountPeakIntensityThreshold = " << parameters->filterCountPeakIntensityThreshold
				<< std::endl;

	}
	std::cout << std::endl;

	std::cout << "# SPECTRUM PROCESSING OPTIONS" << std::endl;
	{
		std::cout << "filterRemovePeakIntensityThreshold = " << parameters->filterRemovePeakIntensityThreshold
				<< std::endl;
		std::cout << "filterMaxPeaksUsed = " << parameters->filterMaxPeaksUsed << std::endl;
		std::cout << "filterMaxDynamicRange = " << parameters->filterMaxDynamicRange << std::endl;

		std::cout << "#(m/z)^<peakScalingMzPower> * (raw intensity)^<peakScalingIntensityPower" << std::endl;
		std::cout << "peakScalingMzPower = " << parameters->peakScalingMzPower << std::endl;
		std::cout << "peakScalingIntensityPower = " << parameters->peakScalingIntensityPower << std::endl;

		std::cout << "peakScalingUnassignedPeaks = " << parameters->peakScalingUnassignedPeaks << std::endl;
		std::cout << "peakBinningFractionToNeighbor = " << parameters->peakBinningFractionToNeighbor << std::endl;
		std::cout << "filterLibMaxPeaksUsed = " << parameters->filterLibMaxPeaksUsed << std::endl;
		std::cout << "filterLightIonsMzThreshold = " << parameters->filterLightIonsMzThreshold << std::endl;

	}
	std::cout << "-------------END-PARAMETERS--------------\n" << std::endl;
}

#define DEF_PRINT_PARAM(KEY,VALUE) (*outputStream) << "<parameter name=\"" << KEY << "\" value=\"" << VALUE << "\"/>" << std::endl;

void printPepXMLSearchParams(Parameters::PStruct *parameters, std::ofstream *outputStream)
{
	assert(parameters != NULL);

	//(*outputStream) << "<!-- input parameters -->" << std::endl;
	// GENERAL OPTIONS
	{
		DEF_PRINT_PARAM("use_sp4_scoring", "true")
		//We rely on the old SpectraST algorithm
		DEF_PRINT_PARAM("spectral_library", parameters->libraryFile)
		DEF_PRINT_PARAM("params_file", parameters->parameterFile)
		DEF_PRINT_PARAM("max_gpu_memory_usage", parameters->maxGPUMemoryUsage)
	}

	//CANDIDATE SELECTION AND SCORING OPTIONS
	{
		DEF_PRINT_PARAM("precursor_mz_tolerance", parameters->precursorMzTolerance)
		DEF_PRINT_PARAM("detect_homologs_max_rank", parameters->detectHomologs)
		DEF_PRINT_PARAM("fval_fraction_delta", parameters->fvalFractionDelta)
		DEF_PRINT_PARAM("fval_use_dot_bias", std::boolalpha << parameters->fvalUseDotBias)
	}

	//OUTPUT AND DISPLAY OPTIONS
	{
//	 std::cout << "hitListShowMaxRank = " << parameters->hitListShowMaxRank << std::endl;
//	 std::cout << "outputExtension = " << Parameters::Output::OutMap.at(parameters->outputExtension) << std::endl;
	}

	//SPECTRUM FILTERING OPTIONS
	{
		DEF_PRINT_PARAM("filter_min_peak_count", parameters->filterMinPeakCount)
		DEF_PRINT_PARAM("filter_all_peaks_below_mz", parameters->filterAllPeaksBelowMz)
		DEF_PRINT_PARAM("filter_max_intensity_below", parameters->filterMaxIntensityBelow)
		DEF_PRINT_PARAM("filter_min_mz_range", parameters->filterMinMzRange)
		DEF_PRINT_PARAM("filter_count_peak_intensity_threshold", parameters->filterCountPeakIntensityThreshold)
	}

	//SPECTRUM PROCESSING OPTIONS
	{
//		std::cout << "filterRemovePeakIntensityThreshold = " << parameters->filterRemovePeakIntensityThreshold
//				<< std::endl;

		DEF_PRINT_PARAM("filter_max_peaks_used", parameters->filterMaxPeaksUsed)
		DEF_PRINT_PARAM("filter_max_dynamic_range", parameters->filterMaxDynamicRange)
		DEF_PRINT_PARAM("peak_scaling_mz_power", parameters->peakScalingMzPower)
		DEF_PRINT_PARAM("peak_scaling_intensity_power", parameters->peakScalingIntensityPower)
		DEF_PRINT_PARAM("peak_scaling_unassigned_peaks", parameters->peakScalingUnassignedPeaks)
		DEF_PRINT_PARAM("peak_binning_fraction_to_neighbor", parameters->peakBinningFractionToNeighbor)
		//TODO: parameterize
		DEF_PRINT_PARAM("peak_binning_num_bins_per_mz", 1)
		DEF_PRINT_PARAM("peak_no_binning", "false")

		DEF_PRINT_PARAM("filter_lib_max_peaks_used", parameters->filterLibMaxPeaksUsed)
		DEF_PRINT_PARAM("filter_light_ions_mz_threshold", parameters->filterLightIonsMzThreshold)

	}
	/*


	 fout << "<parameter name=\"precursor_mz_use_average\" value=\"" << (precursorMzUseAverage ? "true" : "false") << "\"/>" << endl;
	 //  fout << "<parameter name=\"precursor_mz_use_average\" value=\"" << (indexRetrievalUseAverage ? "true" : "false") << "\"/>" << endl;

	 fout << "<parameter name=\"use_P_value\" value=\"" << (usePValue ? "true" : "false") << "\"/>" << endl;

	 fout << "<parameter name=\"use_tierwise_open_modification_search\" value=\"" << (useTierwiseOpenModSearch ? "true" : "false") << "\"/>" << endl;

	 fout << "<parameter name=\"filter_iTRAQ_reporter_peaks\" value=\"" << (filterITRAQReporterPeaks ? "true" : "false") << "\"/>" << endl;

	 fout << "<parameter name=\"filter_TMT_reporter_peaks\" value=\"" << (filterTMTReporterPeaks ? "true" : "false") << "\"/>" << endl;

	 */
}

//helper function to set the help flag to true and give an error if an value was left empty while an value was expected
inline bool isValueEmpty(Parameters::PStruct *parameters, std::string key, std::string value)
{
	if (value.empty())
	{
		std::cerr << "No value is set for '" << key << "'" << std::endl;
		//If no value is specified
		parameters->flagShowHelp = true;
		return true;
	}
	return false;
}

/*
 * Sets a parameter according to a key value
 * Global function to handle the parameters, to keep it all together
 */
void setParameter(Parameters::PStruct *parameters, std::string key, std::string value)
{
	assert(parameters != NULL);

#ifdef DEBUG
	printf("Debug: setParameter('%s' = '%s')\n", key.c_str(), value.c_str());
#endif

	using namespace Parameters;

	//delete any whitespace in the input
	key.erase(0, key.find_first_not_of(" "));
	key.erase(key.find_last_not_of(" ") + 1);
	value.erase(0, value.find_first_not_of(" "));
	value.erase(value.find_last_not_of(" ") + 1);

	//Temp variables for getting the numbers from the
	float fTempFloat;
	int iTempInt;
	PEnum variable; // for assiging test variable

	//Map the arguments/keys according to the initialised enum map (made in initializeEnumMapping)
	switch (PMap.find(key)->second)
	{
	/*
	 * Add query files
	 */
	case queryFiles:
		if (isValueEmpty(parameters, key, value)) break;
		parameters->queryFiles.push_back(value);
		break;

		/*
		 * CONVERTING OPTIONS
		 */
	case convertSpLib:
		//Assign true if only -C is specified
		if (value.empty()) value.assign("true");
		parameters->convertSpLib = parseBoolFromString(value);
//		parameters->queryFiles.clear(); //clear any queryFiles that were specified in the param file
		break;

		/*
		 * TESTING OPTIONS
		 */
	case testFvalConfidence:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->testFvalConfidence = fTempFloat;
		else
			std::cerr << "'testFvalConfidence' should be equal or greater than 0.0" << std::endl;
		//Set the output to testing
		parameters->normalOutput = false;
		break;
	case testVariable:
		//Dummy value, so break
		if (value == "none") break;

		variable = PMap.find(value)->second;
		if (variable == unknown)
		{
			std::cout << "Unkown variable to loop over (" << value << ")" << std::endl;
			break;
		}
		parameters->testVariable = value;
		//Set the output to testing
		parameters->normalOutput = false;
		break;
	case testMinValueVariable:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->testMinValueVariable = fTempFloat;
		else
			std::cerr << "'testMinValueVariable' should be equal or greater than 0.0" << std::endl;
		//Set the output to testing
		parameters->normalOutput = false;
		break;
	case testMaxValueVariable:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat > 0.0)
			parameters->testMaxValueVariable = fTempFloat;
		else
			std::cerr << "'testMaxValueVariable' should be greater than 0.0" << std::endl;
		//Set the output to testing
		parameters->normalOutput = false;
		break;
	case testValueVariableStepSize:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat > 0.0)
			parameters->testValueVariableStepSize = fTempFloat;
		else
			std::cerr << "'testValueVariableStepSize' should be greater than 0.0" << std::endl;
		//Set the output to testing
		parameters->normalOutput = false;
		break;

		/*
		 * GENERAL OPTIONS
		 */
	case flagShowHelp:
		//Assign true if only -H is specified
		if (value.empty()) value.assign("true");
		parameters->flagShowHelp = parseBoolFromString(value);
		break;
	case showParameters:
		//Assign true if only -DispP is specified
		if (value.empty()) value.assign("true");
		parameters->showParameters = parseBoolFromString(value);
		break;
	case normalOutput:
		//Do not show normal output if one of the testing parameters is set
		if (value.empty()) value.assign("true");
		parameters->normalOutput = parseBoolFromString(value);
		break;
	case libraryFile:
		if (isValueEmpty(parameters, key, value)) break;
/*		if (value.compare(value.size() - libraryExtension.size(), libraryExtension.size(), libraryExtension) != 0)
		{
			std::cerr << "'libraryFile' must be a library file ending with " << libraryExtension << std::endl;
			break;
		}
*/		parameters->libraryFile = value;
		break;
	case parameterFile:
		parameters->parameterFile = value;
		break;
	case maxGPUMemoryUsage:
		if (isValueEmpty(parameters, key, value)) break;
		{
			uint32_t tempLong = atol(value.c_str());
			//1048576 = 1MB
			if (tempLong >= 1048576uL && tempLong <= 1073741824uL)
				parameters->maxGPUMemoryUsage = tempLong;
			else
				std::cerr
						<< "'maxGPUMemoryUsage' should be equal or greater than 1048576 (1MB) and smaller than 1073741824 (1GB)"
						<< std::endl;
		}
		break;

		/*
		 * CANDIDATE SELECTION AND SCORING OPTIONS
		 */
	case precursorMzTolerance:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->precursorMzTolerance = fTempFloat;
		else
			std::cerr << "'precursorMzTolerance' should be equal or greater than 0.0" << std::endl;
		break;
	case detectHomologs:
		if (isValueEmpty(parameters, key, value)) break;
		iTempInt = atoi(value.c_str());
		if (iTempInt >= 0)
			parameters->detectHomologs = iTempInt;
		else
			std::cerr << "'detectHomologs' should be equal or greater than 0" << std::endl;
		break;
	case fvalFractionDelta:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->fvalFractionDelta = fTempFloat;
		else
			std::cerr << "'fvalFractionDelta' should be equal or greater than 0.0" << std::endl;
		break;
	case fvalUseDotBias:
		//Assign true if only -H is specified
		if (value.empty()) value.assign("true");
		parameters->fvalUseDotBias = parseBoolFromString(value);
		break;

		/*
		 * OUTPUT AND DISPLAY OPTIONS
		 */
	case hitListTopHitFvalThreshold:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->hitListTopHitFvalThreshold = fTempFloat;
		else
			std::cerr << "'hitListTopHitFvalThreshold' should be equal or greater than 0.0" << std::endl;
		break;
	case hitListLowerHitsFvalThreshold:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->hitListLowerHitsFvalThreshold = fTempFloat;
		else
			std::cerr << "'hitListLowerHitsFvalThreshold' should be equal or greater than 0.0" << std::endl;
		break;
		//Maximum rank for hits shown for each query (number of hits to show (default: 1)
	case hitListShowMaxRank:
		if (isValueEmpty(parameters, key, value)) break;
		iTempInt = atoi(value.c_str());
		if (iTempInt >= 1)
			parameters->hitListShowMaxRank = iTempInt;
		else
			std::cerr << "'hitListShowMaxRank' should be equal or greater than 1" << std::endl;
		break;
	case outputExtension:
		if (isValueEmpty(parameters, key, value)) break;
		//Transform to lower case extension
		std::transform(value.begin(), value.end(), value.begin(), ::tolower);
		{ //limit the scope of iOutputExtension
			Output::Types iOutputExtension = Output::unknown;
			if (value.compare("xml") == 0 || value.compare("pepxml") == 0 || value.compare("pep.xml") == 0)
			{
				iOutputExtension = Output::xml;
			}
			else if (value.compare("xls") == 0)
			{
				iOutputExtension = Output::xls;
			}
			parameters->outputExtension = iOutputExtension;
		}
		break;
		/*
		 * SPECTRUM FILTERING OPTIONS
		 */
	case filterMinPeakCount:
		if (isValueEmpty(parameters, key, value)) break;
		iTempInt = atoi(value.c_str());
		if (iTempInt >= 0)
			parameters->filterMinPeakCount = iTempInt;
		else
			std::cerr << "'filterMinPeakCount' should be equal or greater than 0" << std::endl;
		break;
	case filterAllPeaksBelowMz:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->filterAllPeaksBelowMz = fTempFloat;
		else
			std::cerr << "'filterAllPeaksBelowMz' should be equal or greater than 0.0" << std::endl;
		break;
	case filterMinMzRange:
		if (isValueEmpty(parameters, key, value)) break;
		iTempInt = atoi(value.c_str());
		if (iTempInt >= 0)
			parameters->filterMinMzRange = iTempInt;
		else
			std::cerr << "'filterMinMzRange' should be equal or greater than 0" << std::endl;
		break;
	case filterMaxIntensityBelow:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->filterMaxIntensityBelow = fTempFloat;
		else
			std::cerr << "'filterMaxIntensityBelow' should be equal or greater than 0.0" << std::endl;
		break;
	case filterCountPeakIntensityThreshold:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->filterCountPeakIntensityThreshold = fTempFloat;
		else
			std::cerr << "'filterCountPeakIntensityThreshold' should be equal or greater than 0.0" << std::endl;
		break;

		/*
		 * SPECTRUM PROCESSING OPTIONS
		 */
	case filterRemovePeakIntensityThreshold:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->filterRemovePeakIntensityThreshold = fTempFloat;
		else
			std::cerr << "'filterRemovePeakIntensityThreshold' should be equal or greater than 0.0" << std::endl;
		break;
	case filterMaxPeaksUsed:
		if (isValueEmpty(parameters, key, value)) break;
		iTempInt = atoi(value.c_str());
		if (iTempInt >= 1)
			parameters->filterMaxPeaksUsed = iTempInt;
		else
			std::cerr << "'filterMaxPeaksUsed' should be equal or greater than 1" << std::endl;
		break;
	case filterMaxDynamicRange:
		if (isValueEmpty(parameters, key, value)) break;
		iTempInt = atoi(value.c_str());
		if (iTempInt >= 1)
			parameters->filterMaxDynamicRange = iTempInt;
		else
			std::cerr << "'filterMaxDynamicRange' should be equal or greater than 1" << std::endl;
		break;
	case peakScalingMzPower:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->peakScalingMzPower = fTempFloat;
		else
			std::cerr << "'peakScalingMzPower' should be equal or greater than 0.0" << std::endl;
		break;
	case peakScalingIntensityPower:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->peakScalingIntensityPower = fTempFloat;
		else
			std::cerr << "'peakScalingIntensityPower' should be equal or greater than 0.0" << std::endl;
		break;
	case peakBinningFractionToNeighbor:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->peakBinningFractionToNeighbor = fTempFloat;
		else
			std::cerr << "'peakBinningFractionToNeighbor' should be equal or greater than 0.0" << std::endl;
		break;
	case peakScalingUnassignedPeaks:
		if (isValueEmpty(parameters, key, value)) break;
		fTempFloat = atof(value.c_str());
		if (fTempFloat >= 0.0)
			parameters->peakScalingUnassignedPeaks = fTempFloat;
		else
			std::cerr << "'peakScalingUnassignedPeaks' should be equal or greater than 0.0" << std::endl;
		break;
	case filterLibMaxPeaksUsed:
		if (isValueEmpty(parameters, key, value)) break;
		iTempInt = atoi(value.c_str());
		if (iTempInt >= 1)
			parameters->filterLibMaxPeaksUsed = iTempInt;
		else
			std::cerr << "'filterLibMaxPeaksUsed' should be equal or greater than 1" << std::endl;
		break;
	case filterLightIonsMzThreshold:
		if (isValueEmpty(parameters, key, value)) break;
		iTempInt = atoi(value.c_str());
		if (iTempInt >= 0)
			parameters->filterLightIonsMzThreshold = iTempInt;
		else
			std::cerr << "'filterLightIonsMzThreshold' should be equal or greater than 0" << std::endl;
		break;

		/*
		 * Ignore unidentified arguments
		 */
	default:
		//Show help when an unidentified argument is passed.
//		parameters->flagShowHelp = true;
		std::cerr << "Unknown argument, ignored ('" << key.c_str() << "' = '" << value.c_str() << "'), for help -H"
				<< std::endl;
		break;
	}
}

/*
 * Reads parameters from a file
 */
void readParameterFile(Parameters::PStruct *parameters)
{
	assert(parameters != NULL);

	std::ifstream parameterFileStream;
	parameterFileStream.open(parameters->parameterFile.c_str());

//Check if file stream exists
	if (parameterFileStream.is_open())
	{
		std::string line;
		while (std::getline(parameterFileStream, line))
		{
			//delete starting whitespaces and tabs
			line.erase(0, line.find_first_not_of(" "));
			line.erase(0, line.find_first_not_of("	"));

			//Skip empty lines
			if (line.empty()) continue;

			//check if it is not a comment (# or //), skip these
			if (line.compare(0, 1, "#") == 0 || line.compare(0, 2, "//") == 0) continue;

			size_t equalPostion = line.find("=");
			if (equalPostion != std::string::npos)
			{
				//Retrieve the key and value of each parameter and add to the parameter struct
				setParameter(parameters, line.substr(0, equalPostion), line.substr(equalPostion + 1));
			}
			else
			{
				//Warn the user about bad formatted parameters
				std::cerr << "Warning: bad parameter line '" << line << "'" << std::endl;
			}

		}
	}
	else
	{
		//File stream could not be opened
		std::cerr << "Warning: Could not open params file (" << parameters->parameterFile
				<< "). Skipping reading params file" << std::endl;
	}

	parameterFileStream.close();
}

/*
 * Reads parameters from the cli arguments
 */
void readParamsFromArguments(Parameters::PStruct *parameters, std::vector<std::string> cliArguments)
{
	assert(parameters != NULL);

	std::string prevargument("");
	//Loop over all the arguments
	for (std::vector<std::string>::iterator argument = cliArguments.begin(); argument != cliArguments.end(); ++argument)
	{
		if ((*argument).compare(0, 1, "-") == 0)
		{
			std::string key((*argument));
			std::string value("");
			//If the next argument does not start with a - then it might possibly be a value for the option,
			//	otherwise it will be empty
			if ((argument + 1) != cliArguments.end() && (*(argument + 1)).compare(0, 1, "-") != 0)
			{
				value.assign(*(argument + 1));
			}
			//Set the parameter of the cli argument with the general function
			setParameter(parameters, key, value);
		}
		else if (!prevargument.empty() != 0 && prevargument.compare(0, 1, "-") != 0) //If this argument is no -Argument and the previous one neither, then this is a input file
		{
			//Add the file to the parameters
			setParameter(parameters, "queryFile", *argument);
		}

		prevargument.assign(*argument);
	}
}

/*
 * Puts the parameters of the default, file or command line into the parameter struct
 * and returns the pointer to this struct (must be deleted with destroyParameters!)
 */
Parameters::PStruct* getParameters(int argc, char *argv[])
{
	struct Parameters::PStruct *parameters = new Parameters::PStruct();
	assert(parameters != NULL);

	initializeEnumMapping();

	setDefaults(parameters);

	std::vector<std::string> cliArguments(argv, argv + argc);

	//Check if a parameter file is specified
	bool hasParameterFile = false;
	for (std::vector<std::string>::iterator argument = cliArguments.begin(); argument != cliArguments.end(); ++argument)
	{
		//Check if -F for parameter file and check if the argument thereafter is specified
		if (Parameters::PMap.find(*argument)->second == Parameters::parameterFile)
		{
			hasParameterFile = true;
			//If the next argument starts with - then assume the default parameter file, otherwise the one specified here
			if ((argument + 1) != cliArguments.end() && (*(argument + 1)).compare(0, 1, "-") != 0)
			{
				setParameter(parameters, "parameterFile", *(argument + 1));
			}
			else
			{
				//Give a warning for using the default parameter file
				std::cerr << "Warning: No file is specified for the parameters '" << parameters->parameterFile
						<< "' assumed" << std::endl;
			}
			//Check for only the first parameter argument and exclude/skip the rest
			break;
		}
	}

	if (hasParameterFile)
	{
		//Read the parameters from the file if necessary and overwrite any default values
		readParameterFile(parameters);
	}

	//Reading all the parameters from the command line arguments, overwriting the default values and specified in the params file
	readParamsFromArguments(parameters, cliArguments);

	return parameters;
}

/*
 * Free the parameter struct
 */
void destroyParameters(Parameters::PStruct *parameters)
{
	assert(parameters != NULL);
	delete parameters;
}

/*
 * Helper function to parse a boolean (1,0,true,false) from string input.
 */
bool parseBoolFromString(std::string string)
{
	bool bBoolean;
	if (string.compare("1") == 0 || string.compare("0") == 0)
	{
		//If the bool is represented as a int
		std::istringstream(string) >> std::noboolalpha >> bBoolean;
	}
	else
	{
		//If the bool is represented as 'true' or 'false'
		std::istringstream(string) >> std::boolalpha >> bBoolean;
	}
	return bBoolean;
}
