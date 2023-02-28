#include <dispatch_care_correct_cpu.hpp>
#include <hpc_helpers.cuh>
#include <config.hpp>
#include <options.hpp>
#include <singlehashminhasher.hpp>
#include <correct_cpu.hpp>

#include <minhasherlimit.hpp>

#include <correctedsequence.hpp>
#include <correctionresultoutput.hpp>
#include <sequencehelpers.hpp>
#include <cpuminhasherconstruction.hpp>
#include <ordinaryminhasher.hpp>
#include <serializedobjectstorage.hpp>
#include <sortserializedresults.hpp>
#include <chunkedreadstorageconstruction.hpp>
#include <chunkedreadstorage.hpp>

#include <contiguousreadstorage.hpp>
#include <vector>
#include <iostream>
#include <mutex>
#include <thread>
#include <memory>

#include <experimental/filesystem>


namespace filesys = std::experimental::filesystem;



namespace care{

    template<class T>
    void printDataStructureMemoryUsage(const T& datastructure, const std::string& name){
    	auto toGB = [](std::size_t bytes){
    			    double gb = bytes / 1024. / 1024. / 1024.0;
    			    return gb;
    		    };

        auto memInfo = datastructure.getMemoryInfo();
        
        std::cout << name << " memory usage: " << toGB(memInfo.host) << " GB on host\n";
        for(const auto& pair : memInfo.device){
            std::cout << name << " memory usage: " << toGB(pair.second) << " GB on device " << pair.first << '\n';
        }
    }

    void performCorrection(ProgramOptions programOptions){

        std::cout << "Running CARE CPU" << std::endl;

        helpers::CpuTimer step1Timer("STEP1");

        std::cout << "STEP 1: Database construction" << std::endl;

        helpers::CpuTimer buildReadStorageTimer("build_readstorage");

        std::unique_ptr<ChunkedReadStorage> cpuReadStorage = constructChunkedReadStorageFromFiles(programOptions);

        buildReadStorageTimer.print();

        std::cout << "Determined the following read properties:\n";
        std::cout << "----------------------------------------\n";
        std::cout << "Total number of reads: " << cpuReadStorage->getNumberOfReads() << "\n";
        std::cout << "Minimum sequence length: " << cpuReadStorage->getSequenceLengthLowerBound() << "\n";
        std::cout << "Maximum sequence length: " << cpuReadStorage->getSequenceLengthUpperBound() << "\n";
        std::cout << "----------------------------------------\n";

        if(programOptions.save_binary_reads_to != ""){
            std::cout << "Saving reads to file " << programOptions.save_binary_reads_to << std::endl;
            helpers::CpuTimer timer("save_to_file");
            cpuReadStorage->saveToFile(programOptions.save_binary_reads_to);
            timer.print();
            std::cout << "Saved reads" << std::endl;
        }
        
        if(programOptions.autodetectKmerlength){
            const int maxlength = cpuReadStorage->getSequenceLengthUpperBound();

            auto getKmerSizeForHashing = [](int maximumReadLength){
                if(maximumReadLength < 160){
                    return 20;
                }else{
                    return 32;
                }
            };

            programOptions.kmerlength = getKmerSizeForHashing(maxlength);

            std::cout << "Will use k-mer length = " << programOptions.kmerlength << " for hashing.\n";
        }

        std::cout << "Reads with ambiguous bases: " << cpuReadStorage->getNumberOfReadsWithN() << std::endl;        

        printDataStructureMemoryUsage(*cpuReadStorage, "reads");

        //compareMaxRssToLimit(programOptions.memoryTotalLimit, "Error memorylimit after cpureadstorage");


        helpers::CpuTimer buildMinhasherTimer("build_minhasher");

        auto minhasherAndType = constructCpuMinhasherFromCpuReadStorage(
            programOptions,
            *cpuReadStorage,
            CpuMinhasherType::Ordinary
        );

        //compareMaxRssToLimit(programOptions.memoryTotalLimit, "Error memorylimit after cpuminhasher");


        CpuMinhasher* cpuMinhasher = minhasherAndType.first.get();

        buildMinhasherTimer.print();

        std::cout << "Using minhasher type: " << to_string(minhasherAndType.second) << "\n";
        std::cout << "CpuMinhasher can use " << cpuMinhasher->getNumberOfMaps() << " maps\n";

        if(cpuMinhasher->getNumberOfMaps() <= 0){
            std::cout << "Cannot construct a single cpu hashtable. Abort!" << std::endl;
            return;
        }

        if(programOptions.mustUseAllHashfunctions 
            && programOptions.numHashFunctions != cpuMinhasher->getNumberOfMaps()){
            std::cout << "Cannot use specified number of hash functions (" 
                << programOptions.numHashFunctions <<")\n";
            std::cout << "Abort!\n";
            return;
        }

        if(minhasherAndType.second == CpuMinhasherType::Ordinary){

            OrdinaryCpuMinhasher* ordinaryCpuMinhasher = dynamic_cast<OrdinaryCpuMinhasher*>(cpuMinhasher);
            assert(ordinaryCpuMinhasher != nullptr);

            if(programOptions.save_hashtables_to != "") {
                std::cout << "Saving minhasher to file " << programOptions.save_hashtables_to << std::endl;
                std::ofstream os(programOptions.save_hashtables_to);
                assert((bool)os);
                helpers::CpuTimer timer("save_to_file");
                ordinaryCpuMinhasher->writeToStream(os);
                timer.print();

                std::cout << "Saved minhasher" << std::endl;
            }

        }

        printDataStructureMemoryUsage(*cpuMinhasher, "hash tables");

        step1Timer.print();

        std::cout << "STEP 2: Error correction" << std::endl;

        helpers::CpuTimer step2Timer("STEP2");

        auto partialResults = cpu::correct_cpu(
            programOptions, 
            *cpuMinhasher, 
            *cpuReadStorage
        );

        step2Timer.print();

        std::cout << "Correction throughput : ~" << (cpuReadStorage->getNumberOfReads() / step2Timer.elapsed()) << " reads/second.\n";

        std::cerr << "Constructed " << partialResults.size() << " corrections. ";
        std::cerr << "They occupy a total of " << (partialResults.dataBytes() + partialResults.offsetBytes()) << " bytes\n";

        //compareMaxRssToLimit(programOptions.memoryTotalLimit, "Error memorylimit after correction");

        minhasherAndType.first.reset();
        cpuMinhasher = nullptr;        
        cpuReadStorage.reset();

        if(programOptions.correctionType != CorrectionType::Print && programOptions.correctionTypeCands != CorrectionType::Print){

            //Merge corrected reads with input file to generate output file

            const std::size_t availableMemoryInBytes = getAvailableMemoryInKB() * 1024;
            const auto partialResultMemUsage = partialResults.getMemoryInfo();

            // std::cerr << "availableMemoryInBytes = " << availableMemoryInBytes << "\n";
            // std::cerr << "memoryLimitOption = " << programOptions.memoryTotalLimit << "\n";
            // std::cerr << "partialResultMemUsage = " << partialResultMemUsage.host << "\n";

            std::size_t memoryForSorting = std::min(
                availableMemoryInBytes,
                programOptions.memoryTotalLimit - partialResultMemUsage.host
            );

            if(memoryForSorting > 1*(std::size_t(1) << 30)){
                memoryForSorting = memoryForSorting - 1*(std::size_t(1) << 30);
            }
            //std::cerr << "memoryForSorting = " << memoryForSorting << "\n"; 

            std::cout << "STEP 3: Constructing output file(s)" << std::endl;

            helpers::CpuTimer step3Timer("STEP3");

            helpers::CpuTimer sorttimer("sort_results_by_read_id");

            sortSerializedResultsByReadIdAscending<EncodedTempCorrectedSequence>(
                partialResults,
                memoryForSorting
            );

            sorttimer.print();
            
            std::vector<FileFormat> formats;
            for(const auto& inputfile : programOptions.inputfiles){
                formats.emplace_back(getFileFormat(inputfile));
            }
            std::vector<std::string> outputfiles;
            for(const auto& outputfilename : programOptions.outputfilenames){
                outputfiles.emplace_back(programOptions.outputdirectory + "/" + outputfilename);
            }
            FileFormat outputFormat = formats[0];
            if(programOptions.gzoutput){
                if(outputFormat == FileFormat::FASTQ){
                    outputFormat = FileFormat::FASTQGZ;
                }else if(outputFormat == FileFormat::FASTA){
                    outputFormat = FileFormat::FASTAGZ;
                }
            }else{
                if(outputFormat == FileFormat::FASTQGZ){
                    outputFormat = FileFormat::FASTQ;
                }else if(outputFormat == FileFormat::FASTAGZ){
                    outputFormat = FileFormat::FASTA;
                }
            }
            constructOutputFileFromCorrectionResults(
                programOptions.inputfiles, 
                partialResults, 
                outputFormat,
                outputfiles,
                programOptions.showProgress,
                programOptions
            );

            step3Timer.print();

            //compareMaxRssToLimit(programOptions.memoryTotalLimit, "Error memorylimit after output construction");

            std::cout << "Construction of output file(s) finished." << std::endl;
        }

    }

}
