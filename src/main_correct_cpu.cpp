#include <version.hpp>
#include <config.hpp>

#include <cxxopts/cxxopts.hpp>
#include <options.hpp>
#include <dispatch_care_correct_cpu.hpp>

#include <threadpool.hpp>

#include <readlibraryio.hpp>
#include <concurrencyhelpers.hpp>

#include <fstream>
#include <iostream>
#include <ios>
#include <string>
#include <omp.h>

#include <experimental/filesystem>

#include "main_correct_cpu.h"

namespace filesys = std::experimental::filesystem;

using namespace care;

void printCommandlineArguments(std::ostream& out, const cxxopts::ParseResult& parseresults){

	const auto args = parseresults.arguments();
	for(const auto& opt : args){
		out << opt.key() << '=' << opt.value() << '\n';
	}
}

bool checkMandatoryArguments(const cxxopts::ParseResult& parseresults){

	const std::vector<std::string> mandatory = {
		"inputfiles", "outdir", "outputfilenames", "coverage",
		"pairmode"
	};

	bool success = true;
	for(const auto& opt : mandatory){
		if(parseresults.count(opt) == 0){
			success = false;
			std::cerr << "Mandatory argument " << opt << " is missing.\n";
		}
	}

	return success;
}

template<class T>
std::string tostring(const T& t){
	return std::to_string(t);
}

template<>
std::string tostring(const bool& b){
	return b ? "true" : "false";
}

int main_correction(int argc, char** argv, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *Q, std::atomic_int *producerDone) {
    for(int i = 0; i < argc; i++) 
        printf("%s ", argv[i]);
    printf("\n");

	bool help = false;
	bool showVersion = false;

	cxxopts::Options commandLineOptions(argv[0], "CARE: Context-Aware Read Error Correction for Illumina reads");

	addMandatoryOptions(commandLineOptions);
	addMandatoryOptionsCorrect(commandLineOptions);
	addMandatoryOptionsCorrectCpu(commandLineOptions);

	addAdditionalOptions(commandLineOptions);
	addAdditionalOptionsCorrect(commandLineOptions);
	addAdditionalOptionsCorrectCpu(commandLineOptions);

	commandLineOptions.add_options("Additional")			
		("help", "Show this help message", cxxopts::value<bool>(help))
		("version", "Print version", cxxopts::value<bool>(showVersion));

	auto parseresults = commandLineOptions.parse(argc, argv);

	if(showVersion){
		std::cout << "CARE version " << CARE_VERSION_STRING << std::endl;
		std::exit(0);
	}

	if(help) {
		std::cout << commandLineOptions.help({"", "Mandatory", "Additional"}) << std::endl;
		std::exit(0);
	}

	const bool mandatoryPresent = checkMandatoryArguments(parseresults);
	if(!mandatoryPresent){
		std::cout << commandLineOptions.help({"Mandatory"}) << std::endl;
		std::exit(0);
	}

	ProgramOptions programOptions(parseresults);

	programOptions.batchsize = 16;

	if(!programOptions.isValid()) throw std::runtime_error("Invalid program options!");

	if(programOptions.useQualityScores){
		const bool hasQ = std::all_of(
			programOptions.inputfiles.begin(),
			programOptions.inputfiles.end(),
			[](const auto& s){
				return hasQualityScores(s);
			}
		);

		if(!hasQ){
			std::cerr << "Quality scores have been disabled because there exist reads in an input file without quality scores.\n";
			
			programOptions.useQualityScores = false;
		}
	}

	if(programOptions.correctionType == CorrectionType::Forest){
		if(programOptions.mlForestfileAnchor == ""){
			std::cerr << "CorrectionType is not set to Classic, but no valid classifier file is provided. Abort!\n";
			std::exit(0);
		}

		if(programOptions.mlForestfileCands == ""){
			programOptions.mlForestfileCands = programOptions.mlForestfileAnchor;
		}
	}

	std::cout << std::boolalpha;
	std::cout << "CARE CORRECT CPU will be started with the following parameters:\n";
	std::cout << "----------------------------------------\n";

	programOptions.printMandatoryOptions(std::cout);
	programOptions.printMandatoryOptionsCorrect(std::cout);
	programOptions.printMandatoryOptionsCorrectCpu(std::cout);

	programOptions.printAdditionalOptions(std::cout);
	programOptions.printAdditionalOptionsCorrect(std::cout);
	programOptions.printAdditionalOptionsCorrectCpu(std::cout);

	std::cout << "----------------------------------------\n";
	std::cout << std::noboolalpha;

    const int numThreads = programOptions.threads;

	omp_set_num_threads(numThreads);

    care::performCorrectionOutToQueue(programOptions, Q, producerDone);

	return 0;
}
