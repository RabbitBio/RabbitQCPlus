#include <readlibraryio.hpp>

#include <config.hpp>
#include <sequencehelpers.hpp>
#include <util.hpp>

#include <iterator>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <cassert>
#include <future>

#include <zlib.h>
#include <fcntl.h> // open

namespace care{


//###### BEGIN WRITER IMPLEMENTATION

void SequenceFileWriter::writeRead(const std::string& name, const std::string& comment, const std::string& sequence, const std::string& quality){
    writeReadImpl(name, comment, sequence, quality);
}

void SequenceFileWriter::writeRead(const std::string& header, const std::string& sequence, const std::string& quality){
    writeReadImpl(header, sequence, quality);
}

void SequenceFileWriter::writeRead(const Read& read){
    //std::cerr << "Write " << header << "\n" << sequence << " " << "\n" << quality << "\n";
    writeRead(read.header, read.sequence, read. quality);
}

UncompressedWriter::UncompressedWriter(const std::string& filename, FileFormat format)
        : SequenceFileWriter(filename, format){

    assert(format == FileFormat::FASTA || format == FileFormat::FASTQ);

    ofs = std::ofstream(filename);
    if(!ofs){
        throw std::runtime_error("Cannot open file " + filename + " for writing.");
    }

    isFastq = format == FileFormat::FASTQ || format == FileFormat::FASTQGZ;
    delimHeader = '>';
    if(isFastq){
        delimHeader = '@';
    }

}

void UncompressedWriter::writeReadImpl(const std::string& name, const std::string& comment, const std::string& sequence, const std::string& quality){
    ofs << delimHeader << name;
    if(comment.length() > 0){
        ofs << ' ' << comment;
    }
    ofs << '\n' << sequence << '\n';
    if(format == FileFormat::FASTQ){
        ofs << '+' << '\n'
            << quality << '\n';
    }
}

void UncompressedWriter::writeReadImpl(const std::string& header, const std::string& sequence, const std::string& quality){
    ofs << delimHeader << header;
    ofs << '\n' << sequence << '\n';
    if(format == FileFormat::FASTQ){
        ofs << '+' << '\n'
            << quality << '\n';
    }
}

void UncompressedWriter::writeImpl(const std::string& data){
    ofs << data;
}

GZipWriter::GZipWriter(const std::string& filename, FileFormat format)
        : SequenceFileWriter(filename, format){

    assert(format == FileFormat::FASTAGZ || format == FileFormat::FASTQGZ);

    fp = gzopen(filename.c_str(), "w");
    if(fp == NULL){
        throw std::runtime_error("Cannot open file " + filename + " for writing.");
    }

    isFastq = format == FileFormat::FASTQ || format == FileFormat::FASTQGZ;
    delimHeader = '>';
    if(isFastq){
        delimHeader = '@';
    }

    charbuffer.reserve(2*1024*1024);
}

GZipWriter::~GZipWriter(){
    // if(numBufferedReads > 0){
    //     writeBufferedReads();
    // }
    // numBufferedReads = 0;
    flush();
    gzclose(fp);
}

void GZipWriter::writeReadImpl(const std::string& header, const std::string& sequence, const std::string& quality){
    writeReadImpl(header, "", sequence, quality);
}

void GZipWriter::writeReadImpl(const std::string& name, const std::string& comment, const std::string& sequence, const std::string& quality){
    write({&delimHeader, 1});
    write(name);
    if(comment.size() > 0){
        write(" ");
        write(comment);
    }
    write("\n");
    write(sequence);
    write("\n");
    if(isFastq){
        write("+");
        write("\n");
        write(quality);
        write("\n");
    }
}

void GZipWriter::writeImpl(const std::string& data){
    write(data);
}


//###### END WRITER IMPLEMENTATION



    bool hasQualityScores(const std::string& filename){
        kseqpp::KseqPP reader(filename);

        const int n = 5;
        int i = 0;
        int count = 0;
        while(reader.next() >= 0 && i < n){
            if(reader.getCurrentQuality().size() > 0){
                count++;
            }
            i++;
        }
        if(count > 0 && count == i){
            return true;
        }else if(count == 0){
            return false;
        }else{
            throw std::runtime_error("Error. Some reads do not have quality scores");
        }
    }

    FileFormat getFileFormat(const std::string& filename){
        const bool gzip = kseqpp::hasGzipHeader(filename);
        const bool qscore = hasQualityScores(filename);

        if(gzip){
            if(qscore){
                return FileFormat::FASTQGZ;
            }else{
                return FileFormat::FASTAGZ;
            }
        }else{
            if(qscore){
                return FileFormat::FASTQ;
            }else{
                return FileFormat::FASTA;
            }
        }
    }

    std::unique_ptr<SequenceFileWriter> makeSequenceWriter(const std::string& filename, FileFormat fileFormat){
        switch (fileFormat) {
        case FileFormat::FASTA:
        case FileFormat::FASTQ:
            return std::make_unique<UncompressedWriter>(filename, fileFormat);
        case FileFormat::FASTAGZ:
        case FileFormat::FASTQGZ:
            return std::make_unique<GZipWriter>(filename, fileFormat);
        case FileFormat::NONE:
    	default:
    		throw std::runtime_error("makeSequenceWriter: invalid format.");
    	}
    }

    SequenceFileProperties getSequenceFileProperties(const std::string& filename){
        return getSequenceFileProperties(filename, false);
    }

    SequenceFileProperties getSequenceFileProperties(const std::string& filename, bool printProgress){
        SequenceFileProperties prop;

        prop.maxSequenceLength = 0;
        prop.minSequenceLength = std::numeric_limits<int>::max();

        auto showProgress = [&](auto totalCount, auto seconds){
            if(printProgress){
                std::cout << "Found " << totalCount << " reads. Elapsed time: " << seconds << " seconds.\n";
            }
        };

        auto updateShowProgressInterval = [](auto duration){
            return duration * 2;
        };

        ProgressThread<read_number> progressThread(
            std::numeric_limits<read_number>::max(), 
            showProgress, 
            updateShowProgressInterval
        );

		std::uint64_t totalCount = 0;

        forEachReadInFile(
            filename, 
            [&](auto /*readNumber*/, const auto& read){
                int len = read.sequence.length();
                if(len > prop.maxSequenceLength)
                    prop.maxSequenceLength = len;
                if(len < prop.minSequenceLength)
                    prop.minSequenceLength = len;

                ++totalCount;

                progressThread.addProgress(1);
            }
        );

        progressThread.finished();

        prop.nReads = totalCount;


        return prop;
    }

	std::uint64_t getNumberOfReads(const std::string& filename){

        std::uint64_t count = 0;
        forEachReadInFile(
            filename, 
            [&](auto /*readNumber*/, auto /*read*/){
                count++;
            }
        );

        return count;
	}



}
