/*
 * params.h
 *
 *  Created on: Jan 6, 2015
 */

#ifndef PARAMS_H_
#define PARAMS_H_

#include <string>
#include <vector>
#include <map>
#include <stdint.h>

#define MZ_START 10
#define MZ_END 2010
#define MZ_RANGE (MZ_END - MZ_START + 1)

/*
 * When adding a parameter it must be added in the following places:
 * 	Parameters::PStruct
 * 	Parameters::PEnum
 * 	initializeEnumMapping()
 * 	showHelpParameters()
 * 	setDefaults()
 * 	printParameters()
 * 	printPepXMLSearchParams()
 * 	setParameter()
 */

namespace Parameters {

/*
 * Output settings with abstraction of the output type
 */
namespace Output {
enum Types
{
	unknown, xml, xls
};
static std::map<Types, std::string> OutMap =
	{
		{ unknown, "unknown" },
		{ xml, "pep.xml" },
		{ xls, "xls" } };
} //PARAMETER::OUTPUT NAMESPACE

/*
 * Parameter struct that holds all the options.
 * They are organized according to SpectraST
 */
struct PStruct
{
		std::vector<std::string> queryFiles;

		/*
		 * CONVERTING OPTIONS
		 * Converting splib to fastpasslib
		 */
		bool convertSpLib;

		/*
		 * TESTING OPTIONS
		 */
		//Confidence interval for counting hits (default: 0.85)
		float testFvalConfidence;
		//These variables define which variable will be varied over a specific range with step size
		std::string testVariable;
		float testMinValueVariable;
		float testMaxValueVariable;
		float testValueVariableStepSize;

		/*
		 * GENERAL OPTIONS
		 */
		//Flag if help needs to be shown (for command line options)
		bool flagShowHelp;
		//Flag if help needs to be shown (for command line options)
		bool showParameters;
		//Flag if normal output must be written to a file
		bool normalOutput;

		//Mandatory unless specified in parameter file. <file> must have .fplib extension.
		std::string libraryFile;
		//Parameter file, if it is omitted, “fastpass.params” is assumed
		std::string parameterFile;
		//Parameter file, if it is omitted, “fastpass.params” is assumed
		uint32_t maxGPUMemoryUsage;

		/*
		 * CANDIDATE SELECTION AND SCORING OPTIONS
		 */
		//Specify precursor m/z tolerance in Th. Monoisotopic mass is assumed (default: 3.0)
		float precursorMzTolerance;
		//Detect homologous lower hits up to <rank>. Looks for lower hits homologous to the top hit and adjust delta accordingly. (Default: is 4)
		unsigned int detectHomologs;
		//Specify the fraction of the normalized delta score (delta/dot) in the F-value formula. (default: 0.4)
		float fvalFractionDelta;
		//Use dot bias to penalize high-scoring matches with massive noise and/or dominant peak. (default: true)
		bool fvalUseDotBias;

		/*
		 * OUTPUT AND DISPLAY OPTIONS
		 */
		//Minimum F value threshold for the top hit. Only top hits having F value greater than <thres> will be printed.	(Default = 0.03)
		float hitListTopHitFvalThreshold;
		//Minimum F value threshold for the lower hits. Only lower hits having F value greater than <thres> will be printed. (Default = 0.45)
		float hitListLowerHitsFvalThreshold;
		//Maximum rank for hits shown for each query (number of hits to show (default: 1)
		unsigned hitListShowMaxRank;
		Output::Types outputExtension;

		/*
		 * SPECTRUM FILTERING OPTIONS
		 */
		//Discard query spectra with fewer than X peaks above threshold (default:10)
		unsigned filterMinPeakCount;
		//Filter peaks in query spectra below this M/Z (default:520)
		unsigned filterAllPeaksBelowMz;
		//Discard query spectra with no peaks with intensity above <inten>.	(Default is 0)
		float filterMaxIntensityBelow;
		//Discard query spectra with m/z range narrower than <range>.	(Default is 350)
		unsigned filterMinMzRange;
		/*
		 * Minimum peak intensity for peaks to be counted. Only peaks with intensity above <thres>
		 * will be counted to meet the requirement for minimum number of peaks (default: 2.01)
		 */
		float filterCountPeakIntensityThreshold;

		/*
		 * SPECTRUM PROCESSING OPTIONS
		 */
		//Noise peak threshold. All peaks with intensities below <thres> will be zeroed. (default: 2.01)
		float filterRemovePeakIntensityThreshold;
		//Remove all but the top X peaks in query spectra (default:150)
		unsigned filterMaxPeaksUsed;
		//Remove all peaks smaller than 1/<num> of the base (highest) peak in query spectra. (default: 1000)
		unsigned filterMaxDynamicRange;
		//Intensity scaling power with respect to the m/z value and the raw intensity. The scaled intensity will be (m/z)^<mzpow> * (raw intensity)^<intpow> (Default is <mzpow> = 0.0, <intpow> = 0.5)
		float peakScalingMzPower;
		float peakScalingIntensityPower;
		//Scaling factor for unassigned peaks in library spectra. Unassigned peaks in the library spectra will be scaled by <factor>.	(Default is 1.0)
		float peakScalingUnassignedPeaks;
		//Fraction of the scaled intensity assigned to neighboring bins (default: 0.5)
		float peakBinningFractionToNeighbor;
		//Remove all but the top <num> peaks in the LIBRARY spectra.	(Default is 50)
		unsigned filterLibMaxPeaksUsed;
		//Remove all light ions with m/z lower than <thres> Th for both (library, not yet) and query spectra.
		unsigned filterLightIonsMzThreshold;
};

/*
 * Parameter enum for mapping arguments to general int for Switch statement
 * They are organized according to SpectraST
 *
 * Put in namespace to avoid collisions
 */

enum PEnum
{
	unknown = 0, //has value of 0, to accommodate for missing arguments

	queryFiles,

	/*
	 * Converting splib to fastpasslib
	 */
	convertSpLib,

	/*
	 * TESTING OPTIONS
	 */
	testFvalConfidence,
	testVariable, testMinValueVariable, testMaxValueVariable,	testValueVariableStepSize,

	/*
	 * GENERAL OPTIONS
	 */
	flagShowHelp,
	showParameters,
	normalOutput,
	libraryFile,
	parameterFile,
	maxGPUMemoryUsage,

	/*
	 * CANDIDATE SELECTION AND SCORING OPTIONS
	 */
	precursorMzTolerance,
	detectHomologs,
	fvalFractionDelta,
	fvalUseDotBias,

	/*
	 * OUTPUT AND DISPLAY OPTIONS
	 */
	hitListTopHitFvalThreshold,
	hitListLowerHitsFvalThreshold,
	hitListShowMaxRank,
	outputExtension,

	/*
	 * SPECTRUM FILTERING OPTIONS
	 */
	filterMinPeakCount,
	filterAllPeaksBelowMz,
	filterMaxIntensityBelow,
	filterMinMzRange,
	filterCountPeakIntensityThreshold,

	/*
	 * SPECTRUM PROCESSING OPTIONS
	 */
	filterRemovePeakIntensityThreshold,
	filterMaxPeaksUsed,
	filterMaxDynamicRange,
	peakScalingMzPower,
	peakScalingIntensityPower,
	peakScalingUnassignedPeaks,
	peakBinningFractionToNeighbor,
	filterLibMaxPeaksUsed,
	filterLightIonsMzThreshold,
};

// Map to associate the strings with the enum values
static std::map<std::string, PEnum> PMap;

static std::string libraryExtension (".fpbin");

} //PARAMETER NAMESPACE

/*
 * Sets a parameter according to a key value
 * Global function to handle the parameters, to keep it all together
 */
void setParameter(Parameters::PStruct *parameters, std::string key, std::string value);

/*
 * Sets the default parameters, reading if set a parameter file and
 * then sets the command line arguments.
 * So the Command line arguments have the highest priority
 */
Parameters::PStruct* getParameters(int argc, char *argv[]);

/*
 * Prints the parameter struct to stdout with printf
 */
void printParameters(Parameters::PStruct *parameters);

/*
 * Prints the parameter struct in pepxml format to a stream
 */
void printPepXMLSearchParams(Parameters::PStruct *parameters, std::ofstream *outputStream);

/*
 * Free the parameter struct
 */
void destroyParameters(Parameters::PStruct *parameters);

/*
 * Shows command line help when options for the possible options
 */
void showHelpParameters();

/*
 * Helper function to parse a boolean (1,0,true,false) from string input.
 */
bool parseBoolFromString(std::string string);

#endif /* PARAMS_H_ */
