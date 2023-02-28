#ifndef CARE_SEQUENCEFILEIO_HPP
#define CARE_SEQUENCEFILEIO_HPP

#include <config.hpp>

#include <kseqpp/kseqpp.hpp>

#include <hpc_helpers.cuh>

#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <exception>
#include <memory>
#include <sstream>



namespace care{

enum class FileFormat {FASTA, FASTQ, FASTAGZ, FASTQGZ, NONE};

struct SequenceFileProperties{
        std::uint64_t nReads{};
        int minSequenceLength{};
        int maxSequenceLength{};
};


struct Read {
	std::string header = "";
	std::string sequence = "";
	std::string quality = "";

	bool operator==(const Read& other) const
	{
		return (header == other.header 
                 && sequence == other.sequence 
                && quality == other.quality);
	}
	bool operator!=(const Read& other) const
	{
		return !(*this == other);
	}

	void reset()
	{
		header.clear();
		sequence.clear();
		quality.clear();
	}
};

struct ReadWithId{
    int fileId;
    std::uint64_t readIdInFile;
    std::uint64_t globalReadId;
    Read read;
};

struct MultiInputReader{
    int inputFileId{};
    std::int64_t readIdInFile{};
    std::int64_t globalReadId{};
    ReadWithId current{};
    std::vector<kseqpp::KseqPP> readerVector{};
    std::vector<std::string> filenames{};

    MultiInputReader() = default;

    MultiInputReader(std::vector<std::string> inputfilenames)
        : filenames(std::move(inputfilenames))
    {
        for(const auto& inputfile : filenames){
            readerVector.emplace_back(std::move(kseqpp::KseqPP{inputfile}));
        }
    }

    int next(){
        const int numFiles = readerVector.size();
        if(inputFileId >= numFiles){ //all files have been processed
            return -1;
        }

        //repeat until a read was retrieved or all files are processed
        while(true){
            const int status = readerVector[inputFileId].next();

            if(status >= 0){
                std::swap(current.read.header, readerVector[inputFileId].getCurrentHeader());
                std::swap(current.read.sequence, readerVector[inputFileId].getCurrentSequence());
                std::swap(current.read.quality, readerVector[inputFileId].getCurrentQuality());
                current.fileId = inputFileId;
                current.readIdInFile = readIdInFile;
                current.globalReadId = globalReadId;

                readIdInFile++;
                globalReadId++;

                return status;
            }else{
                inputFileId++;
                readIdInFile = 0;

                if(inputFileId >= numFiles){ //all files have been processed
                    return -1;
                }
                
            }
        }
    }

    ReadWithId& getCurrent(){
        return current;
    }
};


struct PairedInputReader{
    std::int64_t readIdInFile{};
    std::int64_t globalReadId{};
    ReadWithId current1{};
    ReadWithId current2{};
    std::vector<kseqpp::KseqPP> readerVector{};
    std::vector<std::string> filenames{};

    PairedInputReader() = default;

    PairedInputReader(std::vector<std::string> inputfilenames)
        : filenames(std::move(inputfilenames))
    {
        assert(filenames.size() > 0);
        assert(filenames.size() <= 2);

        for(const auto& inputfile : filenames){
            readerVector.emplace_back(std::move(kseqpp::KseqPP{inputfile}));
        }
    }

    int next(){

        const int status1 = readerVector[0].next();
        if(status1 >= 0){
            std::swap(current1.read.header, readerVector[0].getCurrentHeader());
            std::swap(current1.read.sequence, readerVector[0].getCurrentSequence());
            std::swap(current1.read.quality, readerVector[0].getCurrentQuality());
            current1.fileId = 0;
            current1.readIdInFile = readIdInFile;
            current1.globalReadId = globalReadId;

            globalReadId++;
        }else{
            return -1;
        }

        const int fileIdForMate = readerVector.size() > 1 ? 1 : 0;

        if(fileIdForMate == 0){
            readIdInFile++;
        }

        const int status2 = readerVector[fileIdForMate].next();
        if(status2 >= 0){
            std::swap(current2.read.header, readerVector[fileIdForMate].getCurrentHeader());
            std::swap(current2.read.sequence, readerVector[fileIdForMate].getCurrentSequence());
            std::swap(current2.read.quality, readerVector[fileIdForMate].getCurrentQuality());
            current2.fileId = fileIdForMate;
            current2.readIdInFile = readIdInFile;
            current2.globalReadId = globalReadId;

            globalReadId++;
        }else{
            return -1;
        }

        readIdInFile++;

        return 0;
    }

    ReadWithId& getCurrent1(){
        return current1;
    }

    ReadWithId& getCurrent2(){
        return current2;
    }
};

struct SequenceFileWriter{

    SequenceFileWriter(const std::string& filename_, FileFormat format_) : filename(filename_), format(format_)
	{

	};
	virtual ~SequenceFileWriter()
	{

	}

    void writeRead(const Read& read);

    void writeRead(const std::string& name, const std::string& comment, const std::string& sequence, const std::string& quality);

    void writeRead(const std::string& header, const std::string& sequence, const std::string& quality);

    virtual void writeReadImpl(const std::string& name, const std::string& comment, const std::string& sequence, const std::string& quality) = 0;

    virtual void writeReadImpl(const std::string& header, const std::string& sequence, const std::string& quality) = 0;

    virtual void writeImpl(const std::string& data) = 0;

protected:

    std::string filename;
    FileFormat format;

};

struct UncompressedWriter : public SequenceFileWriter{
    UncompressedWriter(const std::string& filename, FileFormat format);

    void writeReadImpl(const std::string& name, const std::string& comment, const std::string& sequence, const std::string& quality) override;
    void writeReadImpl(const std::string& header, const std::string& sequence, const std::string& quality) override;
    void writeImpl(const std::string& data) override;

    bool isFastq;
    char delimHeader;

    std::ofstream ofs;
};

struct GZipWriter : public SequenceFileWriter{
    static constexpr int maxBufferedBytes = 4*1024*1024;

    GZipWriter(const std::string& filename, FileFormat format);
    ~GZipWriter();

    void writeReadImpl(const std::string& name, const std::string& comment, const std::string& sequence, const std::string& quality) override;
    void writeReadImpl(const std::string& header, const std::string& sequence, const std::string& quality) override;
    void writeImpl(const std::string& data) override;

private:

    void flush(){
        gzwrite(fp, charbuffer.data(), charbuffer.size());
        charbuffer.clear();
    }

    void write(std::string_view data){
        std::size_t remaining = data.size();
        auto src = data.begin();
        while(remaining > 0){
            const std::size_t free = charbuffer.size() < maxBufferedBytes ? maxBufferedBytes - charbuffer.size() : 0;
            const std::size_t toCopy = std::min(free, remaining);
            charbuffer.insert(charbuffer.end(), src, src + toCopy);
            src += toCopy;
            remaining -= toCopy;
            if(charbuffer.size() >= maxBufferedBytes){
                flush();
            }
        }
    }



    bool isFastq;
    char delimHeader;

    std::vector<char> charbuffer;

    gzFile fp;
};



std::unique_ptr<SequenceFileWriter> makeSequenceWriter(const std::string& filename, FileFormat fileFormat);

bool hasQualityScores(const std::string& filename);
FileFormat getFileFormat(const std::string& filename);


SequenceFileProperties getSequenceFileProperties(const std::string& filename);
SequenceFileProperties getSequenceFileProperties(const std::string& filename, bool printProgress);

std::uint64_t getNumberOfReads(const std::string& filename);

template<class Func>
void forEachReadInFile(const std::string& filename, Func f){

    kseqpp::KseqPP reader(filename);

    Read read;

    std::int64_t readNumber = 0;

    auto getNextRead = [&](){
        const int status = reader.next();
        //std::cerr << "parser status = 0 in file " << filenames[i] << '\n';
        if(status >= 0){
            #if 0
                read.name = reader.getCurrentName();
                read.comment = reader.getCurrentComment();
                read.sequence = reader.getCurrentSequence();
                read.quality = reader.getCurrentQuality();
            #else
                std::swap(read.header, reader.getCurrentHeader());
                std::swap(read.sequence, reader.getCurrentSequence());
                std::swap(read.quality, reader.getCurrentQuality());
            #endif
        }else if(status < -1){
            std::cerr << "parser error status " << status << " in file " << filename << '\n';
        }

        bool success = (status >= 0);

        return success;
    };

    bool success = getNextRead();

    while(success){

        f(readNumber, read);
        
        readNumber++;

        success = getNextRead();
    }
}



template<class Func>
void forEachReadInPairedFiles(const std::string& file1, const std::string& file2, Func f){

    kseqpp::KseqPP reader1(file1);
    kseqpp::KseqPP reader2(file2);

    int which = 0;

    Read read;

    std::int64_t readNumber = 0;

    auto getNextRead = [&](){
        kseqpp::KseqPP* ptr = &reader1;
        if(which == 1){
            ptr = &reader2;
        }

        auto& reader = *ptr;
        const int status = reader.next();
        //std::cerr << "parser status = 0 in file " << filenames[i] << '\n';
        if(status >= 0){
            #if 0
                read.name = reader.getCurrentName();
                read.comment = reader.getCurrentComment();
                read.sequence = reader.getCurrentSequence();
                read.quality = reader.getCurrentQuality();
            #else
                std::swap(read.header, reader.getCurrentHeader());
                std::swap(read.sequence, reader.getCurrentSequence());
                std::swap(read.quality, reader.getCurrentQuality());
            #endif
        }else if(status < -1){
            std::cerr << "parser error status " << status << " in file " << (which == 0 ? file1 : file2) << '\n';
        }

        bool success = (status >= 0);

        if(which == 0){
            which = 1;
        }else{
            which = 0;
        }

        return success;
    };

    bool success = getNextRead();

    while(success){

        f(readNumber, read);
        
        readNumber++;

        success = getNextRead();
    }
}







} //end namespace

#endif
