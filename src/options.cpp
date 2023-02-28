#include <options.hpp>
#include <hpc_helpers.cuh>
#include <util.hpp>
#include <config.hpp>
#include <readlibraryio.hpp>
#include <memorymanagement.hpp>
#include <filehelpers.hpp>

#include <iostream>
#include <thread>
#include <string>
#include <stdexcept>

#include <experimental/filesystem>

namespace filesys = std::experimental::filesystem;

namespace care{


    std::string to_string(SequencePairType s){
        switch (s)
        {
        case SequencePairType::Invalid:
            return "Invalid";
        case SequencePairType::SingleEnd:
            return "SingleEnd";
        case SequencePairType::PairedEnd:
            return "PairedEnd";
        default:
            return "Error";
        }
    }

    std::string to_string(CorrectionType t)
    {
        switch (t)
        {
        case CorrectionType::Classic:
            return "Classic";
            break;
        case CorrectionType::Forest:
            return "Forest";
            break;
        case CorrectionType::Print:
            return "Print";
            break;
        default:
            return "Forgot to name correction type";
            break;
        }
    }

    ProgramOptions::ProgramOptions(const cxxopts::ParseResult& pr){
        ProgramOptions& result = *this;

        if(pr.count("correctionQualityLabels")){
            result.outputCorrectionQualityLabels = pr["correctionQualityLabels"].as<bool>();
        }

        if(pr.count("minalignmentoverlap")){
            result.min_overlap = pr["minalignmentoverlap"].as<int>();
        }
        if(pr.count("maxmismatchratio")){
            result.maxErrorRate = pr["maxmismatchratio"].as<float>();
        }
        if(pr.count("minalignmentoverlapratio")){
            result.min_overlap_ratio = pr["minalignmentoverlapratio"].as<float>();
        }

        if(pr.count("excludeAmbiguous")){
            result.excludeAmbiguousReads = pr["excludeAmbiguous"].as<bool>();
        }

        if(pr.count("candidateCorrection")){
            result.correctCandidates = pr["candidateCorrection"].as<bool>();
        }

        if(pr.count("correctionTypeCands")){
            const int val = pr["correctionTypeCands"].as<int>();

            switch(val){
                case 1: result.correctionTypeCands = CorrectionType::Forest; break;
                case 2: result.correctionTypeCands = CorrectionType::Print; break;
                default: result.correctionTypeCands = CorrectionType::Classic; break;
            }
        }

        if(pr.count("useQualityScores")){
            result.useQualityScores = pr["useQualityScores"].as<bool>();
        }

        if(pr.count("gzoutput")){
            result.gzoutput = pr["gzoutput"].as<bool>();
        }

        if(pr.count("enforceHashmapCount")){
            result.mustUseAllHashfunctions = pr["enforceHashmapCount"].as<bool>();
        }

        if(pr.count("singlehash")){
            result.singlehash = pr["singlehash"].as<bool>();
        }

        if(pr.count("coverage")){
            result.estimatedCoverage = pr["coverage"].as<float>();
        }

        if(pr.count("errorfactortuning")){
            result.estimatedErrorrate = pr["errorfactortuning"].as<float>();
        }

        if(pr.count("coveragefactortuning")){
            result.m_coverage = pr["coveragefactortuning"].as<float>();
        }

        if(pr.count("kmerlength")){
            result.kmerlength = pr["kmerlength"].as<int>();
            if(result.kmerlength == 0){
                result.autodetectKmerlength = true;
            }else{
                result.autodetectKmerlength = false;
            }
        }else{
            result.autodetectKmerlength = true;
        }

        if(pr.count("hashmaps")){
            result.numHashFunctions = pr["hashmaps"].as<int>();
        }        

        if(pr.count("batchsize")){
            result.batchsize = pr["batchsize"].as<int>();
        }

        if(pr.count("maxForestTreesAnchor")){
            int n = pr["maxForestTreesAnchor"].as<int>();
            if (n>0) result.maxForestTreesAnchor = n;
        }

        if(pr.count("maxForestTreesCands")){
            int n = pr["maxForestTreesCands"].as<int>();
            if (n>0) result.maxForestTreesCands = n;
        }

        if(pr.count("candidateCorrectionNewColumns")){
            result.new_columns_to_correct = pr["candidateCorrectionNewColumns"].as<int>();
        }

        if(pr.count("correctionType")){
            const int val = pr["correctionType"].as<int>();

            switch(val){
                case 1: result.correctionType = CorrectionType::Forest; break;
                case 2: result.correctionType = CorrectionType::Print; break;
                default: result.correctionType = CorrectionType::Classic; break;
            }
        }

        if(pr.count("thresholdAnchor")){
            float t = pr["thresholdAnchor"].as<float>();
            result.thresholdAnchor = t>=1.0?t/100:t;
        }

        if(pr.count("thresholdCands")){
            float t = pr["thresholdCands"].as<float>();
            result.thresholdCands = t>=1.0?t/100:t;
        }

        if(pr.count("samplingRateAnchor")){
            float t = pr["samplingRateAnchor"].as<float>();
            result.sampleRateAnchor = t>1.0?t/100:t;
        }

        if(pr.count("samplingRateCands")){
            float t = pr["samplingRateCands"].as<float>();
            result.sampleRateCands = t>1.0?t/100:t;
        }

        if(pr.count("pairedFilterThreshold")){
            result.pairedFilterThreshold = pr["pairedFilterThreshold"].as<float>();
        }

        if(pr.count("insertsize")){
            result.insertSize = pr["insertsize"].as<int>();
        }

        if(pr.count("insertsizedev")){
            result.insertSizeStddev = pr["insertsizedev"].as<int>();
        }

        if(pr.count("fixedStddev")){
            result.fixedStddev = pr["fixedStddev"].as<int>();
        }

        if(pr.count("fixedStepsize")){
            result.fixedStepsize = pr["fixedStepsize"].as<int>();
        }

        if(pr.count("allowOutwardExtension")){
            result.allowOutwardExtension = pr["allowOutwardExtension"].as<bool>();
        }

        if(pr.count("sortedOutput")){
            result.sortedOutput = pr["sortedOutput"].as<bool>();
        }

        if(pr.count("outputRemaining")){
            result.outputRemainingReads = pr["outputRemaining"].as<bool>();
        }

        if(pr.count("threads")){
            result.threads = pr["threads"].as<int>();
        }
        result.threads = std::min(result.threads, (int)std::thread::hardware_concurrency());
      
        if(pr.count("showProgress")){
            result.showProgress = pr["showProgress"].as<bool>();
        }

        if(pr.count("gpu")){
            result.deviceIds = pr["gpu"].as<std::vector<int>>();
        }

        result.canUseGpu = result.deviceIds.size() > 0;

        if(pr.count("warpcore")){
            result.warpcore = pr["warpcore"].as<int>();
        }

        if(pr.count("replicateGpuData")){
            result.replicateGpuData = pr["replicateGpuData"].as<bool>();
        }

        if(pr.count("fixedNumberOfReads")){
            result.fixedNumberOfReads = pr["fixedNumberOfReads"].as<std::size_t>();
        }

        auto parseMemoryString = [](const auto& string) -> std::size_t{
            if(string.length() > 0){
                std::size_t factor = 1;
                bool foundSuffix = false;
                switch(string.back()){
                    case 'K':{
                        factor = std::size_t(1) << 10; 
                        foundSuffix = true;
                    }break;
                    case 'M':{
                        factor = std::size_t(1) << 20;
                        foundSuffix = true;
                    }break;
                    case 'G':{
                        factor = std::size_t(1) << 30;
                        foundSuffix = true;
                    }break;
                }
                if(foundSuffix){
                    const auto numberString = string.substr(0, string.size()-1);
                    return factor * std::stoull(numberString);
                }else{
                    return std::stoull(string);
                }
            }else{
                return 0;
            }
        };

        if(pr.count("memTotal")){
            const auto memoryTotalLimitString = pr["memTotal"].as<std::string>();
            const std::size_t parsedMemory = parseMemoryString(memoryTotalLimitString);
            const std::size_t availableMemory = getAvailableMemoryInKB() * 1024;

            // user-provided memory limit could be greater than currently available memory.
            result.memoryTotalLimit = std::min(parsedMemory, availableMemory);
        }else{
            std::size_t availableMemoryInBytes = getAvailableMemoryInKB() * 1024;
            if(availableMemoryInBytes > 2*(std::size_t(1) << 30)){
                availableMemoryInBytes = availableMemoryInBytes - 2*(std::size_t(1) << 30);
            }

            result.memoryTotalLimit = availableMemoryInBytes;
        }

        if(pr.count("memHashtables")){
            const auto memoryForHashtablesString = pr["memHashtables"].as<std::string>();
            result.memoryForHashtables = parseMemoryString(memoryForHashtablesString);
        }else{
            std::size_t availableMemoryInBytes = result.memoryTotalLimit;
            if(availableMemoryInBytes > 1*(std::size_t(1) << 30)){
                availableMemoryInBytes = availableMemoryInBytes - 1*(std::size_t(1) << 30);
            }

            result.memoryForHashtables = availableMemoryInBytes;
        }

        result.memoryForHashtables = std::min(result.memoryForHashtables, result.memoryTotalLimit);

        if(pr.count("hashloadfactor")){
            result.hashtableLoadfactor = pr["hashloadfactor"].as<float>();
        }

        if(pr.count("qualityScoreBits")){
            result.qualityScoreBits = pr["qualityScoreBits"].as<int>();
        }

        if(pr.count("outdir")){
		    result.outputdirectory = pr["outdir"].as<std::string>();
        }

        if(pr.count("pairmode")){
            const std::string arg = pr["pairmode"].as<std::string>();

            if(arg == "se" || arg == "SE"){
                result.pairType = SequencePairType::SingleEnd;
            }else if(arg == "pe" || arg == "PE"){
                result.pairType = SequencePairType::PairedEnd;
            }else{
                result.pairType = SequencePairType::Invalid;
            }
        }  

        if(pr.count("eo")){
            result.extendedReadsOutputfilename = pr["eo"].as<std::string>();
        }

        if(pr.count("nReads")){
		    result.nReads = pr["nReads"].as<std::uint64_t>();
        }

        if(pr.count("min_length")){
            result.minimum_sequence_length = pr["min_length"].as<int>();
        }

        if(pr.count("max_length")){
            result.maximum_sequence_length = pr["max_length"].as<int>();
        }

        if(pr.count("save-preprocessedreads-to")){
            result.save_binary_reads_to = pr["save-preprocessedreads-to"].as<std::string>();
        }

        if(pr.count("load-preprocessedreads-from")){
            result.load_binary_reads_from = pr["load-preprocessedreads-from"].as<std::string>();
        }

        if(pr.count("save-hashtables-to")){
            result.save_hashtables_to = pr["save-hashtables-to"].as<std::string>();
        }

        if(pr.count("load-hashtables-from")){
            result.load_hashtables_from = pr["load-hashtables-from"].as<std::string>();
        }

        if(pr.count("tempdir")){
            result.tempdirectory = pr["tempdir"].as<std::string>();
        }else{
            result.tempdirectory = result.outputdirectory;
        }

        if(pr.count("ml-forestfile")){
            result.mlForestfileAnchor = pr["ml-forestfile"].as<std::string>();
        }

        if(pr.count("ml-cands-forestfile")){
            result.mlForestfileCands = pr["ml-cands-forestfile"].as<std::string>();
        }

        if(pr.count("ml-print-forestfile")){
            result.mlForestfilePrintAnchor = pr["ml-print-forestfile"].as<std::string>();
        }

        if(pr.count("ml-cands-print-forestfile")){
            result.mlForestfilePrintCands = pr["ml-cands-print-forestfile"].as<std::string>();
        }

        if(pr.count("inputfiles")){
            result.inputfiles = pr["inputfiles"].as<std::vector<std::string>>();
        }

        if(pr.count("outputfilenames")){
            result.outputfilenames = pr["outputfilenames"].as<std::vector<std::string>>();
        }

    }

    bool ProgramOptions::isValid() const noexcept{
        const ProgramOptions& opt = *this;
        bool valid = true;

        if(opt.maxErrorRate < 0.0f || opt.maxErrorRate > 1.0f){
            valid = false;
            std::cout << "Error: maxmismatchratio must be in range [0.0, 1.0], is " + std::to_string(opt.maxErrorRate) << std::endl;
        }

        if(opt.min_overlap < 1){
            valid = false;
            std::cout << "Error: min_overlap must be > 0, is " + std::to_string(opt.min_overlap) << std::endl;
        }

        if(opt.min_overlap_ratio < 0.0f || opt.min_overlap_ratio > 1.0f){
            valid = false;
            std::cout << "Error: min_overlap_ratio must be in range [0.0, 1.0], is "
                        + std::to_string(opt.min_overlap_ratio) << std::endl;
        }

        if(opt.estimatedCoverage <= 0.0f){
            valid = false;
            std::cout << "Error: estimatedCoverage must be > 0.0, is " + std::to_string(opt.estimatedCoverage) << std::endl;
        }

        if(opt.estimatedErrorrate <= 0.0f){
            valid = false;
            std::cout << "Error: estimatedErrorrate must be > 0.0, is " + std::to_string(opt.estimatedErrorrate) << std::endl;
        }

        if(opt.batchsize < 1 /*|| corOpts.batchsize > 16*/){
            valid = false;
            std::cout << "Error: batchsize must be in range [1, ], is " + std::to_string(opt.batchsize) << std::endl;
        }

        if(opt.numHashFunctions < 1){
            valid = false;
            std::cout << "Error: Number of hashmaps must be >= 1, is " + std::to_string(opt.numHashFunctions) << std::endl;
        }

        if(opt.kmerlength < 0 || opt.kmerlength > max_k<kmer_type>::value){
            valid = false;
            std::cout << "Error: kmer length must be in range [0, " << max_k<kmer_type>::value 
                << "], is " + std::to_string(opt.kmerlength) << std::endl;
        }

        if(opt.insertSize < 0){
            valid = false;
            std::cout << "Error: insert size must be >= 0, is " 
                << opt.insertSize << std::endl;
        }

        if(opt.insertSizeStddev < 0){
            valid = false;
            std::cout << "Error: insert size deviation must be >= 0, is " 
                << opt.insertSizeStddev << std::endl;
        }

        if(opt.fixedStddev < 0){
            valid = false;
            std::cout << "Error: fixedStddev must be >= 0, is " 
                << opt.fixedStddev << std::endl;
        }

        if(opt.fixedStepsize < 0){
            valid = false;
            std::cout << "Error: fixedStepsize must be >= 0, is " 
                << opt.fixedStepsize << std::endl;
        }

        if(opt.threads < 1){
            valid = false;
            std::cout << "Error: threads must be > 0, is " + std::to_string(opt.threads) << std::endl;
        }

        if(opt.qualityScoreBits != 1 && opt.qualityScoreBits != 2 && opt.qualityScoreBits != 8){
            valid = false;
            std::cout << "Error: qualityScoreBits must be 1,2,or 8, is " + std::to_string(opt.qualityScoreBits) << std::endl;
        }

        if(!filesys::exists(opt.tempdirectory)){
            bool created = filesys::create_directories(opt.tempdirectory);
            if(!created){
                valid = false;
                std::cout << "Error: Could not create temp directory" << opt.tempdirectory << std::endl;
            }
        }

        if(!filesys::exists(opt.outputdirectory)){
            bool created = filesys::create_directories(opt.outputdirectory);
            if(!created){
                valid = false;
                std::cout << "Error: Could not create output directory" << opt.outputdirectory << std::endl;
            }
        }

        {
            for(const auto& inputfile : opt.inputfiles){
                std::ifstream is(inputfile);
                if(!(bool)is){
                    valid = false;
                    std::cout << "Error: cannot find input file " << inputfile << std::endl;
                }
            }            
        }

        {
            std::vector<FileFormat> formats;
            for(const auto& inputfile : opt.inputfiles){
                FileFormat f = getFileFormat(inputfile);
                if(f == FileFormat::FASTQGZ)
                    f = FileFormat::FASTQ;
                if(f == FileFormat::FASTAGZ)
                    f = FileFormat::FASTA;
                formats.emplace_back(f);
            }
            bool sameFormats = std::all_of(
                formats.begin()+1, 
                formats.end(), [&](const auto f){
                    return f == formats[0];
                }
            );
            if(!sameFormats){
                valid = false;
                std::cout << "Error: Must not specify both fasta and fastq files!" << std::endl;
            }
        }

        {
            if(opt.outputfilenames.size() > 1 && opt.inputfiles.size() != opt.outputfilenames.size()){
                valid = false;
                std::cout << "Error: An output file name must be specified for each input file. Number of input files : " << opt.inputfiles.size() << ", number of output file names: " << opt.outputfilenames.size() << "\n";
            }
        }

        {
            //Disallow invalid type
            if(opt.pairType == SequencePairType::Invalid){
                valid = false;
                std::cout << "Error: pairmode is invalid." << std::endl;
            }

            //In paired end mode, there must be a single input file with interleaved reads, or exactly two input files, one per direction.
            if(opt.pairType == SequencePairType::PairedEnd){
                const int countOk = opt.inputfiles.size() == 1 || opt.inputfiles.size() == 2;
                if(!countOk){
                    valid = false;
                    std::cout << "Error: Invalid number of input files for selected pairmode 'PairedEnd'." << std::endl;
                }
            }

            //In single end mode, a single file allowed
            if(opt.pairType == SequencePairType::SingleEnd){
                const int countOk = opt.inputfiles.size() == 1;
                if(!countOk){
                    valid = false;
                    std::cout << "Error: Invalid number of input files for selected pairmode 'SingleEnd'." << std::endl;
                }
            }
        }

        return valid;
    }


    void ProgramOptions::printMandatoryOptions(std::ostream& stream) const{
        stream << "Output directory: " << outputdirectory << "\n";
        stream << "Estimated dataset coverage: " << estimatedCoverage << "\n";
        stream << "Input files: ";
        for(auto& s : inputfiles){
            stream << s << ' ';
        }
        stream << "\n";
        stream << "Output file names: ";
        for(auto& s : outputfilenames){
            stream << s << ' ';
        }
        stream << "\n";
        stream << "Paired mode: " << to_string(pairType) << "\n";
    }

    void ProgramOptions::printMandatoryOptionsCorrect(std::ostream&) const{
        //nothing
    }

    void ProgramOptions::printMandatoryOptionsCorrectCpu(std::ostream&) const{
        //nothing
    }

    void ProgramOptions::printMandatoryOptionsCorrectGpu(std::ostream& stream) const{
        stream << "Can use GPU(s): " << canUseGpu << "\n";
        if(canUseGpu){
            stream << "GPU device ids: [";
            for(int id : deviceIds){
                stream << " " << id;
            }
            stream << " ]\n";
        }
    }

    void ProgramOptions::printMandatoryOptionsExtend(std::ostream& stream) const{
        stream << "Insert size: " << insertSize << "\n";
	    stream << "Insert size deviation: " << insertSizeStddev << "\n";
        stream << "Extended reads output file: " << extendedReadsOutputfilename << "\n";
    }

    void ProgramOptions::printMandatoryOptionsExtendCpu(std::ostream&) const{
        //nothing
    }

    void ProgramOptions::printMandatoryOptionsExtendGpu(std::ostream& stream) const{
        stream << "Can use GPU(s): " << canUseGpu << "\n";
        if(canUseGpu){
            stream << "GPU device ids: [";
            for(int id : deviceIds){
                stream << " " << id;
            }
            stream << " ]\n";
        }
    }

    void ProgramOptions::printAdditionalOptions(std::ostream& stream) const{
        stream << "Number of hash tables / hash functions: " << numHashFunctions << "\n";
        if(autodetectKmerlength){
            stream << "K-mer size for hashing: auto\n";
        }else{
            stream << "K-mer size for hashing: " << kmerlength << "\n";
        }
        stream << "Enforce number of hash tables: " << mustUseAllHashfunctions << "\n";
        stream << "Threads: " << threads << "\n";
        stream << "Use quality scores: " << useQualityScores << "\n";
        stream << "Bits per quality score: " << qualityScoreBits << "\n";
        stream << "Exclude ambigious reads: " << excludeAmbiguousReads << "\n";
        stream << "Alignment absolute required overlap: " << min_overlap << "\n";
        stream << "Alignment relative required overlap: " << min_overlap_ratio << "\n";
        stream << "Alignment max relative number of mismatches in overlap: " << maxErrorRate << "\n";
        stream << "errorfactortuning: " << estimatedErrorrate << "\n";
        stream << "coveragefactortuning: " << m_coverage << "\n";
        stream << "Show progress bar: " << showProgress << "\n";
        stream << "Output directory: " << outputdirectory << "\n";
        stream << "Temporary directory: " << tempdirectory << "\n";
        stream << "Save preprocessed reads to file: " << save_binary_reads_to << "\n";
        stream << "Load preprocessed reads from file: " << load_binary_reads_from << "\n";
        stream << "Save hash tables to file: " << save_hashtables_to << "\n";
        stream << "Load hash tables from file: " << load_hashtables_from << "\n";
        stream << "Maximum memory for hash tables: " << memoryForHashtables << "\n";
        stream << "Maximum memory total: " << memoryTotalLimit << "\n";
        stream << "Hashtable load factor: " << hashtableLoadfactor << "\n";
        stream << "Fixed number of reads: " << fixedNumberOfReads << "\n";
        stream << "GZ compressed output: " << gzoutput << "\n";
        //stream << "singlehash: " << singlehash << "\n";
    
    }

    void ProgramOptions::printAdditionalOptionsCorrect(std::ostream& stream) const{
        stream << "Correct candidate reads: " << correctCandidates << "\n";
        stream << "Output correction quality labels: " << outputCorrectionQualityLabels << "\n";
	    stream << "Max shift for candidate correction: " << new_columns_to_correct << "\n";
        stream << "Correction type (anchor): " << int(correctionType) 
		    << " (" << to_string(correctionType) << ")\n";
	    stream << "Correction type (cands): " << int(correctionTypeCands) 
		    << " (" << to_string(correctionTypeCands) << ")\n";
        stream << "ml-forestfile: " << mlForestfileAnchor << "\n";
        stream << "ml-cands-forestfile: " << mlForestfileCands << "\n";
        stream << "classification thresholds: " << thresholdAnchor << " | " << thresholdCands << "\n";
        stream << "anchor sampling rate: " << sampleRateAnchor << "\n";
        stream << "cands sampling rate: " << sampleRateCands << "\n";
        stream << "pairedFilterThreshold: " << pairedFilterThreshold << "\n";
        stream << "maxForestTreesAnchor: " << maxForestTreesAnchor << "\n";
        stream << "maxForestTreesCands: " << maxForestTreesCands << "\n";
    }

    void ProgramOptions::printAdditionalOptionsExtend(std::ostream& stream) const{
        stream << "Allow extension outside of gap: " << allowOutwardExtension << "\n";
        stream << "Sort extended reads: " << sortedOutput << "\n";
        stream << "Output remaining reads: " << outputRemainingReads << "\n";
        stream << "fixedStddev: " << fixedStddev << "\n";
	    stream << "fixedStepsize: " << fixedStepsize << "\n";
    }

    void ProgramOptions::printAdditionalOptionsCorrectCpu(std::ostream& stream) const{
        stream << "ml-print-forestfile: " << mlForestfilePrintAnchor << "\n";
        stream << "ml-cands-print-forestfile: " << mlForestfilePrintCands << "\n";
    }

    void ProgramOptions::printAdditionalOptionsCorrectGpu(std::ostream& stream) const{
        stream << "Batch size: " << batchsize << "\n";
        stream << "Warpcore: " << warpcore << "\n";
	    stream << "Replicate GPU data: " << replicateGpuData << "\n";
    }

    void ProgramOptions::printAdditionalOptionsExtendCpu(std::ostream&) const{
        //nothing
    }

    void ProgramOptions::printAdditionalOptionsExtendGpu(std::ostream& stream) const{
        stream << "Batch size: " << batchsize << "\n";
        stream << "Warpcore: " << warpcore << "\n";
	    stream << "Replicate GPU data: " << replicateGpuData << "\n";
    }








    template<class T>
    std::string tostring(const T& t){
        return std::to_string(t);
    }

    template<>
    std::string tostring(const bool& b){
        return b ? "true" : "false";
    }


    void addMandatoryOptions(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Mandatory")
            ("d,outdir", "The output directory. Will be created if it does not exist yet.", 
            cxxopts::value<std::string>())
            ("c,coverage", "Estimated coverage of input file. (i.e. number_of_reads * read_length / genome_size)", 
            cxxopts::value<float>())
            ("i,inputfiles", 
                "The file(s) to correct. "
                "Fasta or Fastq format. May be gzip'ed. "
                "Repeat this option for each input file (e.g. -i file1.fastq -i file2.fastq). "
                "Must not mix fasta and fastq files. "
                "The collection of input files is treated as a single read library",
                cxxopts::value<std::vector<std::string>>())
            ("o,outputfilenames", 
                "The names of outputfiles. "
                "Repeat this option for each output file (e.g. -o file1_corrected.fastq -o file2_corrected.fastq). "
                "If a single output file is specified, it will contain the concatenated results of all input files. "
                "If multiple output files are specified, the number of output files must be equal to the number of input files. "
                "In this case, output file i will contain the results of input file i. "
                "Output files are uncompressed.", 
                cxxopts::value<std::vector<std::string>>())
            ("pairmode", 
                "Type of input reads."
                "SE / se : Single-end reads"
                "PE / pe : Paired-end reads",
                cxxopts::value<std::string>());
    }

    void addMandatoryOptionsCorrect(cxxopts::Options&){
        //nothing
    }

    void addMandatoryOptionsExtend(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Mandatory")
            ("insertsize", 
                "Insert size for paired reads. -- explanation how insert size is interpreted ---", 
                cxxopts::value<int>())
            ("insertsizedev", 
                "Insert size deviation for paired reads.", 
                cxxopts::value<int>())
            ("eo", 
                "The name of the output file containing extended reads",
                cxxopts::value<std::string>());
    }

    void addMandatoryOptionsCorrectCpu(cxxopts::Options&){
        //nothing
    }

    void addMandatoryOptionsCorrectGpu(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Mandatory")
		    ("g,gpu", "Comma-separated list of GPU device ids to be used. (Example: --gpu 0,1 to use GPU 0 and GPU 1)", cxxopts::value<std::vector<int>>());
    }

    void addMandatoryOptionsExtendCpu(cxxopts::Options&){
        //nothing
    }

    void addMandatoryOptionsExtendGpu(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Mandatory")        
		    ("g,gpu", "Comma-separated list of GPU device ids to be used. (Example: --gpu 0,1 to use GPU 0 and GPU 1)", cxxopts::value<std::vector<int>>());
    }

    void addAdditionalOptions(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Additional")
            ("h,hashmaps", "The requested number of hash maps. Must be greater than 0. "
                "The actual number of used hash maps may be lower to respect the set memory limit. "
                "Default: " + tostring(ProgramOptions{}.numHashFunctions), 
                cxxopts::value<int>())
            ("k,kmerlength", "The kmer length for minhashing. If 0 or missing, it is automatically determined.", cxxopts::value<int>())
            ("enforceHashmapCount",
                "If the requested number of hash maps cannot be fullfilled, the program terminates without error correction. "
                "Default: " + tostring(ProgramOptions{}.mustUseAllHashfunctions),
                cxxopts::value<bool>()->implicit_value("true")
            )
            ("t,threads", "Maximum number of thread to use. Must be greater than 0", cxxopts::value<int>())            
            ("q,useQualityScores", "If set, quality scores (if any) are considered during read correction. "
                "Default: " + tostring(ProgramOptions{}.useQualityScores),
            cxxopts::value<bool>()->implicit_value("true"))
            ("qualityScoreBits", "How many bits should be used to store a single quality score. Allowed values: 1,2,8. If not 8, a lossy compression via binning is used."
                "Default: " + tostring(ProgramOptions{}.qualityScoreBits), cxxopts::value<int>())
            ("excludeAmbiguous", 
                "If set, reads which contain at least one ambiguous nucleotide will not be corrected. "
                "Default: " + tostring(ProgramOptions{}.excludeAmbiguousReads),
            cxxopts::value<bool>()->implicit_value("true"))
            ("maxmismatchratio", "Overlap between anchor and candidate must contain at "
                "most (maxmismatchratio * overlapsize) mismatches. "
                "Default: " + tostring(ProgramOptions{}.maxErrorRate),
            cxxopts::value<float>())
            ("minalignmentoverlap", "Overlap between anchor and candidate must be at least this long. "
                "Default: " + tostring(ProgramOptions{}.min_overlap),
            cxxopts::value<int>())
            ("minalignmentoverlapratio", "Overlap between anchor and candidate must be at least as "
                "long as (minalignmentoverlapratio * candidatelength). "
                "Default: " + tostring(ProgramOptions{}.min_overlap_ratio),
            cxxopts::value<float>())
            ("errorfactortuning", "errorfactortuning. "
                "Default: " + tostring(ProgramOptions{}.estimatedErrorrate),
            cxxopts::value<float>())
            ("coveragefactortuning", "coveragefactortuning. "
                "Default: " + tostring(ProgramOptions{}.m_coverage),
            cxxopts::value<float>())
            ("p,showProgress", "If set, progress bar is shown during correction",
            cxxopts::value<bool>()->implicit_value("true"))
            ("tempdir", "Directory to store temporary files. Default: output directory", cxxopts::value<std::string>())
            ("save-preprocessedreads-to", "Save binary dump of data structure which stores input reads to disk",
            cxxopts::value<std::string>())
            ("load-preprocessedreads-from", "Load binary dump of read data structure from disk",
            cxxopts::value<std::string>())
            ("save-hashtables-to", "Save binary dump of hash tables to disk. Ignored for GPU hashtables.",
            cxxopts::value<std::string>())
            ("load-hashtables-from", "Load binary dump of hash tables from disk. Ignored for GPU hashtables.",
            cxxopts::value<std::string>())
            ("memHashtables", "Memory limit in bytes for hash tables and hash table construction. Can use suffix K,M,G , e.g. 20G means 20 gigabyte. This option is not a hard limit. Default: A bit less than memTotal.",
            cxxopts::value<std::string>())
            ("m,memTotal", "Total memory limit in bytes. Can use suffix K,M,G , e.g. 20G means 20 gigabyte. This option is not a hard limit. Default: All free memory.",
            cxxopts::value<std::string>())
            ("hashloadfactor", "Load factor of hashtables. 0.0 < hashloadfactor < 1.0. Smaller values can improve the runtime at the expense of greater memory usage."
                "Default: " + std::to_string(ProgramOptions{}.hashtableLoadfactor), cxxopts::value<float>())
            ("fixedNumberOfReads", "Process only the first n reads. Default: " + tostring(ProgramOptions{}.fixedNumberOfReads), cxxopts::value<std::size_t>())
            ("singlehash", "Use 1 hashtables with h smallest unique hashes. Default: " + tostring(ProgramOptions{}.singlehash), cxxopts::value<bool>())
            ("gzoutput", "gz compressed output (very slow). Default: " + tostring(ProgramOptions{}.gzoutput), cxxopts::value<bool>());
            
    }

    void addAdditionalOptionsCorrect(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Additional")
            ("correctionQualityLabels", "If set, correction quality label will be appended to output read headers. "
                "Default: " + tostring(ProgramOptions{}.outputCorrectionQualityLabels),
            cxxopts::value<bool>()->implicit_value("true"))
            ("candidateCorrection", "If set, candidate reads will be corrected,too. "
                "Default: " + tostring(ProgramOptions{}.correctCandidates),
            cxxopts::value<bool>()->implicit_value("true"))
            ("candidateCorrectionNewColumns", "If candidateCorrection is set, a candidates with an absolute shift of candidateCorrectionNewColumns compared to anchor are corrected. "
                "Default: " + tostring(ProgramOptions{}.new_columns_to_correct),
            cxxopts::value<int>())
            ("correctionType", "0: Classic, 1: Forest, 2: Print . Print is only supported in the cpu version",
                cxxopts::value<int>()->default_value("0"))
            ("correctionTypeCands", "0: Classic, 1: Forest, 2: Print. Print is only supported in the cpu version",
                cxxopts::value<int>()->default_value("0"))
            ("ml-forestfile", "Path of the Random Forest classifier (Anchor correction)",
                cxxopts::value<std::string>())
            ("ml-cands-forestfile", "Path of the Random Forest classifier (Candidate correction)",
                cxxopts::value<std::string>())
            ("thresholdAnchor", "Classification threshold for anchor classifier (\"Forest\") mode",
                cxxopts::value<float>())
            ("thresholdCands", "Classification threshold for candidates classifier (\"Forest\") mode",
                cxxopts::value<float>())
            ("samplingRateAnchor", "sampling rate for anchor features (print mode)",
                cxxopts::value<float>())
            ("samplingRateCands", "sampling rate for candidates features (print mode)",
                cxxopts::value<float>())
            ("pairedFilterThreshold", "Controls alignment quality of unpaired candidates which can pass the candidate alignment filter. Candidate alignments with (num_mismatches / overlap_size) > threshold are removed.", cxxopts::value<float>())
            ("maxForestTreesAnchor", "Max. no. of forests to load from anchor forest file. (-1 = all)", cxxopts::value<int>())
            ("maxForestTreesCands", "Max. no. of forests to load from candidate forest file. (-1 = all)", cxxopts::value<int>());
    }

    void addAdditionalOptionsExtend(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Additional")
            ("allowOutwardExtension", "Will try to fill the gap and extend to the outside"
                "Default: " + tostring(ProgramOptions{}.allowOutwardExtension), cxxopts::value<bool>()->implicit_value("true"))	
            ("sortedOutput", "Extended reads in output file will be sorted by read id."
                "Default: " + tostring(ProgramOptions{}.sortedOutput), cxxopts::value<bool>()->implicit_value("true"))	
            ("outputRemaining", "Output remaining reads which could not be extended. Will be sorted by read id."
                "Default: " + tostring(ProgramOptions{}.outputRemainingReads), cxxopts::value<bool>()->implicit_value("true"))
            ("fixedStepsize", "fixedStepsize "
                "Default: " + tostring(ProgramOptions{}.fixedStepsize),
            cxxopts::value<int>())
            ("fixedStddev", "fixedStddev "
                "Default: " + tostring(ProgramOptions{}.fixedStddev),
            cxxopts::value<int>());
    }

    void addAdditionalOptionsCorrectCpu(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Additional")
            ("ml-print-forestfile", "The output file for extracted anchor features when correctionType = Print",
                cxxopts::value<std::string>())
            ("ml-cands-print-forestfile", "The output file for extracted candidate features when correctionTypeCands = Print",
                cxxopts::value<std::string>());
    }

    void addAdditionalOptionsCorrectGpu(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Additional")
            ("batchsize", "Number of reads to correct in a single batch. Must be greater than 0. "
			    "Default: " + tostring(ProgramOptions{}.batchsize),
                cxxopts::value<int>())		
            ("warpcore", "Enable warpcore hash tables. 0: Disabled, 1: Enabled. "
                "Default: " + tostring(ProgramOptions{}.warpcore),
                cxxopts::value<int>())		
            ("replicateGpuData", "If a GPU data structure fits into the memory of a single GPU, allow its replication to other GPUs. This can improve the runtime when multiple GPUs are used."
                "Default: " + std::to_string(ProgramOptions{}.replicateGpuData), cxxopts::value<bool>());		
    }

    void addAdditionalOptionsExtendCpu(cxxopts::Options&){
        //nothing	
    }

    void addAdditionalOptionsExtendGpu(cxxopts::Options& commandLineOptions){
        commandLineOptions.add_options("Additional")
            ("batchsize", "Number of reads in a single batch. Must be greater than 0. "
                "Default: " + tostring(ProgramOptions{}.batchsize),
                cxxopts::value<int>())
            ("warpcore", "Enable warpcore hash tables. 0: Disabled, 1: Enabled. "
                "Default: " + tostring(ProgramOptions{}.warpcore),
                cxxopts::value<int>())
            ("replicateGpuData", "If a GPU data structure fits into the memory of a single GPU, allow its replication to other GPUs. This can improve the runtime when multiple GPUs are used."
                "Default: " + std::to_string(ProgramOptions{}.replicateGpuData), cxxopts::value<bool>());	
    }

} //namespace care