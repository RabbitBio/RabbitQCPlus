#ifndef CARE_CORRECTEDSEQUENCE_HPP
#define CARE_CORRECTEDSEQUENCE_HPP

#include <config.hpp>
#include <hpc_helpers.cuh>

#include <sequencehelpers.hpp>

#include <sstream>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace care{

    enum class TempCorrectedSequenceType : int {Anchor, Candidate};

    struct TempCorrectedSequence; //forward declaration

    struct EncodedCorrectionEdit{
        // char b;
        // int p;
        std::uint16_t data;

        EncodedCorrectionEdit() = default;
        
        HOSTDEVICEQUALIFIER
        explicit EncodedCorrectionEdit(int p, char b){
            data = p;
            data = data << 2;
            std::uint16_t enc = SequenceHelpers::encodeBase(b);
            data |= enc;

            // this->p = p;
            // this->b = b;
        }

        HOSTDEVICEQUALIFIER
        explicit EncodedCorrectionEdit(int p, std::uint8_t b){
            data = p;
            data = data << 2;
            data |= b;

            // this->p = p;
            // this->b = convertIntToDNACharNoIf(b);
        }

        HOSTDEVICEQUALIFIER
        bool operator==(const EncodedCorrectionEdit& rhs) const{
            return data == rhs.data;
        }

        HOSTDEVICEQUALIFIER
        bool operator!=(const EncodedCorrectionEdit& rhs) const{
            return !(operator==(rhs));
        }

        HOSTDEVICEQUALIFIER
        int pos() const{
            return data >> 2;
            //return p;
        }

        HOSTDEVICEQUALIFIER
        char base() const{
            std::uint8_t enc = data & 0x03;
            return SequenceHelpers::decodeBase(enc);
            //return b;
        }

        HOSTDEVICEQUALIFIER
        void pos(int i){
            std::uint16_t a = i;
            data = (a << 2) | (data & 0x03);
            //p = i;
        }

        HOSTDEVICEQUALIFIER
        void base(char b){
            std::uint16_t enc = SequenceHelpers::encodeBase(b);
            data = (data & ~0x03) | enc;
            //this->b = b;
        }

        HOSTDEVICEQUALIFIER
        void base(std::uint8_t b){
            std::uint16_t enc = b;
            data = (data & ~0x03) | enc;

            //this->b = convertIntToDNACharNoIf(b);
        }
    };

    struct CorrectionEdit{
        char base_;
        int pos_;

        CorrectionEdit() = default;
        HOSTDEVICEQUALIFIER
        CorrectionEdit(int p, char b) : base_(b), pos_(p){}

        HOSTDEVICEQUALIFIER
        CorrectionEdit(const EncodedCorrectionEdit& rhs) : base_(rhs.base()), pos_(rhs.pos()){}

        HOSTDEVICEQUALIFIER
        bool operator==(const CorrectionEdit& rhs) const{
            return base() == rhs.base() && pos() == rhs.pos();
        }

        HOSTDEVICEQUALIFIER
        bool operator!=(const CorrectionEdit& rhs) const{
            return !(operator==(rhs));
        }

        CorrectionEdit& operator=(const CorrectionEdit& rhs) = default;

        HOSTDEVICEQUALIFIER
        CorrectionEdit& operator=(const EncodedCorrectionEdit& rhs){
            base(rhs.base());
            pos(rhs.pos());
            return *this;
        }

        HOSTDEVICEQUALIFIER
        char base() const noexcept{
            return base_;
        }

        HOSTDEVICEQUALIFIER
        int pos() const noexcept{
            return pos_;
        }

        HOSTDEVICEQUALIFIER
        void base(char b) noexcept{
            base_ = b;
        }

        HOSTDEVICEQUALIFIER
        void pos(int i) noexcept{
            pos_ = i;
        }
    };


    struct EncodedTempCorrectedSequence{
        std::uint32_t encodedflags{}; //contains size of data in bytes, and boolean flags
        read_number readId{};
        std::unique_ptr<std::uint8_t[]> data{};

        EncodedTempCorrectedSequence() = default;
        EncodedTempCorrectedSequence(EncodedTempCorrectedSequence&& rhs){
            *this = std::move(rhs);
        }

        EncodedTempCorrectedSequence& operator=(EncodedTempCorrectedSequence&& rhs){
            encodedflags = std::exchange(rhs.encodedflags, 0);
            readId = std::exchange(rhs.readId, 0);
            data = std::move(rhs.data);

            return *this;
        }

        EncodedTempCorrectedSequence(const EncodedTempCorrectedSequence& rhs){
            *this = rhs;
        }

        EncodedTempCorrectedSequence& operator=(const EncodedTempCorrectedSequence& rhs){
            encodedflags = rhs.encodedflags;
            readId = rhs.readId;

            const int numBytes = rhs.getNumBytes();
            data = std::make_unique<std::uint8_t[]>(numBytes);
            std::memcpy(data.get(), rhs.data.get(), numBytes);

            return *this;
        }

        // EncodedTempCorrectedSequence& operator=(const TempCorrectedSequence& rhs){
        //     rhs.encodeInto(*this);

        //     return *this;
        // }

        bool writeToBinaryStream(std::ostream& os) const{
            //assert(bool(os)); 
            os.write(reinterpret_cast<const char*>(&readId), sizeof(read_number));
            //assert(bool(os));
            os.write(reinterpret_cast<const char*>(&encodedflags), sizeof(std::uint32_t));
            //assert(bool(os));
            const int numBytes = getNumBytes();
            os.write(reinterpret_cast<const char*>(data.get()), sizeof(std::uint8_t) * numBytes);
            //assert(bool(os));
            return bool(os);
        }

        bool readFromBinaryStream(std::istream& is){
            is.read(reinterpret_cast<char*>(&readId), sizeof(read_number));
            is.read(reinterpret_cast<char*>(&encodedflags), sizeof(std::uint32_t));
            const int numBytes = getNumBytes();

            data = std::make_unique<std::uint8_t[]>(numBytes);

            is.read(reinterpret_cast<char*>(data.get()), sizeof(std::uint8_t) * numBytes);

            return bool(is);
        }

        std::uint8_t* copyToContiguousMemory(std::uint8_t* ptr, std::uint8_t* endPtr) const{
            const int dataBytes = getNumBytes();

            const std::size_t availableBytes = std::distance(ptr, endPtr);
            const std::size_t requiredBytes = sizeof(read_number) + sizeof(std::uint32_t) + dataBytes;
            if(requiredBytes <= availableBytes){
                std::memcpy(ptr, &readId, sizeof(read_number));
                ptr += sizeof(read_number);
                std::memcpy(ptr, &encodedflags, sizeof(std::uint32_t));
                ptr += sizeof(std::uint32_t);
                std::memcpy(ptr, data.get(), dataBytes);
                ptr += dataBytes;
                return ptr;
            }else{
                return nullptr;
            }        
        }

        void copyFromContiguousMemory(const std::uint8_t* ptr){
            std::memcpy(&readId, ptr, sizeof(read_number));
            ptr += sizeof(read_number);
            std::memcpy(&encodedflags, ptr, sizeof(std::uint32_t));
            ptr += sizeof(read_number);

            const int numBytes = getNumBytes();
            data = std::make_unique<std::uint8_t[]>(numBytes);

            std::memcpy(data.get(), ptr, numBytes);
            //ptr += numBytes;
        }

        bool operator==(const EncodedTempCorrectedSequence& rhs) const{
            const std::uint32_t numBytes = getNumBytes();
            return encodedflags == rhs.encodedflags && readId == rhs.readId 
                    && std::memcmp(data.get(), rhs.data.get(), numBytes);
        }

        bool operator!=(const EncodedTempCorrectedSequence& rhs) const{
            return !(operator==(rhs));
        }

        int getNumBytes() const{
            constexpr std::uint32_t mask = (std::uint32_t(1) << 29)-1;
            return (encodedflags & mask);
        }

        int getSerializedNumBytes() const noexcept{
            return sizeof(read_number) + sizeof(std::uint32_t) + getNumBytes();
        }

        //from serialized object beginning at ptr, return the read id of this object
        static read_number parseReadId(const std::uint8_t* ptr){
            read_number id;
            std::memcpy(&id, ptr, sizeof(read_number));
            return id;
        }

        read_number getReadId() const noexcept{
            return readId;
        }

        bool isHQ() const noexcept{
            return (encodedflags >> 31) & std::uint32_t(1);
        }

        bool useEdits() const noexcept{
            return (encodedflags >> 30) & std::uint32_t(1);
        }

        int getNumEdits() const noexcept{
            if(useEdits()){
                //num edits is stored in the first int of encoded data
                int num;
                std::memcpy(&num, data.get(), sizeof(int));

                return num;
            }else{
                return 0;
            }
        }

        TempCorrectedSequenceType getType() const noexcept{
            return TempCorrectedSequenceType((encodedflags >> 29) & std::uint32_t(1));
        }

        static void encodeDataIntoEncodedCorrectedSequence(
            EncodedTempCorrectedSequence& target,
            read_number readId,
            bool hq,
            bool useEdits,
            TempCorrectedSequenceType type,
            int shift,
            int numEdits,
            const CorrectionEdit* edits,
            int sequenceLength,
            const char* sequence
        ){
            const std::uint32_t oldNumBytes = target.getNumBytes(); 

            target.readId = readId;

            target.encodedflags = (std::uint32_t(hq) << 31);
            target.encodedflags |= (std::uint32_t(useEdits) << 30);
            target.encodedflags |= (std::uint32_t(int(type)) << 29);


            std::uint32_t numBytes = 0;
            if(useEdits){
                numBytes += sizeof(int);
                numBytes += numEdits * (sizeof(int) + sizeof(char));
            }else{
                numBytes += sizeof(int);
                numBytes += sizeof(char) * sequenceLength;
            }

            if(type == TempCorrectedSequenceType::Anchor){
                ; //nothing
            }else{
                //candidate shift
                numBytes += sizeof(int);
            }

            #ifndef NDEBUG
            //flags use 3 bits, remainings bit can be used for encoding
            constexpr std::uint32_t maxNumBytes = (std::uint32_t(1) << 29)-1;
            assert(numBytes <= maxNumBytes);
            #endif

            target.encodedflags |= numBytes;

            if(numBytes > oldNumBytes){
                target.data = std::make_unique<std::uint8_t[]>(numBytes);
            }else{
                ; //reuse buffer
            }

            //fill buffer

            std::uint8_t* ptr = target.data.get();

            if(useEdits){
                std::memcpy(ptr, &numEdits, sizeof(int));
                ptr += sizeof(int);
                for(int i = 0; i < numEdits; i++){
                    const auto& edit = edits[i];
                    const int p = edit.pos();
                    std::memcpy(ptr, &p, sizeof(int));
                    ptr += sizeof(int);
                }
                for(int i = 0; i < numEdits; i++){
                    const auto& edit = edits[i];
                    const char c = edit.base();
                    std::memcpy(ptr, &c, sizeof(char));
                    ptr += sizeof(char);
                }
            }else{
                std::memcpy(ptr, &sequenceLength, sizeof(int));
                ptr += sizeof(int);
                std::memcpy(ptr, sequence, sizeof(char) * sequenceLength);
                ptr += sizeof(char) * sequenceLength;
            }

            if(type == TempCorrectedSequenceType::Anchor){
                ; //nothing
            }else{
                std::memcpy(ptr, &shift, sizeof(int));
                ptr += sizeof(int);
            }
        }

    };

    // represents a sequence produced by the correction of a read.
    // Will be saved to file during correction.
    // Will be loaded from file during mergeResultFiles
    struct TempCorrectedSequence{
        
        TempCorrectedSequence() = default;
        TempCorrectedSequence(const TempCorrectedSequence&) = default;
        TempCorrectedSequence(TempCorrectedSequence&&) = default;
        TempCorrectedSequence& operator=(const TempCorrectedSequence&) = default;
        TempCorrectedSequence& operator=(TempCorrectedSequence&&) = default;

        TempCorrectedSequence(const EncodedTempCorrectedSequence& encoded){
            decode(encoded);
        }

        TempCorrectedSequence& operator=(const EncodedTempCorrectedSequence& encoded){
            decode(encoded);
            return *this;
        }

        bool operator==(const TempCorrectedSequence& rhs) const{
            return hq == rhs.hq && useEdits == rhs.useEdits && type == rhs.type && shift == rhs.shift && readId == rhs.readId
                && sequence == rhs.sequence && edits == rhs.edits;
        }

        bool operator!=(const TempCorrectedSequence& rhs) const{
            return !(operator==(rhs));
        }

        void encodeInto(EncodedTempCorrectedSequence& target) const{
            EncodedTempCorrectedSequence::encodeDataIntoEncodedCorrectedSequence(
                target,
                readId,
                hq,
                useEdits,
                type,
                shift,
                edits.size(),
                edits.data(),
                sequence.size(),
                sequence.data()
            );
        }

        EncodedTempCorrectedSequence encode() const{
            EncodedTempCorrectedSequence encoded;
            encodeInto(encoded);

            return encoded;
        }

        void decode(const EncodedTempCorrectedSequence& encoded){

            readId = encoded.getReadId();

            hq = encoded.isHQ();
            useEdits = encoded.useEdits();
            type = encoded.getType();

            const std::uint8_t* ptr = encoded.data.get();
        

            if(useEdits){
                int size;
                std::memcpy(&size, ptr, sizeof(int));
                ptr += sizeof(int);

                edits.resize(size);

                for(auto& edit : edits){
                    int p;
                    std::memcpy(&p, ptr, sizeof(int));
                    edit.pos(p);
                    ptr += sizeof(int);
                }
                for(auto& edit : edits){
                    char c;
                    std::memcpy(&c, ptr, sizeof(char));
                    edit.base(c);
                    ptr += sizeof(char);
                }
            }else{
                int length;
                std::memcpy(&length, ptr, sizeof(int));
                ptr += sizeof(int);

                sequence.resize(length);
                sequence.replace(0, length, (const char*)ptr, length);

                ptr += sizeof(char) * length;
            }

            if(type == TempCorrectedSequenceType::Anchor){
                ; //nothing
            }else{
                std::memcpy(&shift, ptr, sizeof(int));
                ptr += sizeof(int);
            }
        }

        bool writeToBinaryStream(std::ostream& os) const{
            os.write(reinterpret_cast<const char*>(&readId), sizeof(read_number));
            
            std::uint8_t data = bool(hq);
            data = (data << 1) | bool(useEdits);
            data = (data << 6) | std::uint8_t(int(type));

            os.write(reinterpret_cast<const char*>(&data), sizeof(std::uint8_t));

            if(useEdits){
                os << edits.size() << ' ';
                for(const auto& edit : edits){
                    os << edit.pos() << ' ';
                }
                for(const auto& edit : edits){
                    os << edit.base();
                }
                if(edits.size() > 0){
                    os << ' ';
                }
            }else{
                os << sequence << ' ';
            }

            if(type == TempCorrectedSequenceType::Anchor){
                ; // nothing
            }else{
                os << shift;
            }

            return bool(os);
        }

        bool readFromBinaryStream(std::istream& is){
            std::uint8_t data = 0;

            is.read(reinterpret_cast<char*>(&readId), sizeof(read_number));
            is.read(reinterpret_cast<char*>(&data), sizeof(std::uint8_t));

            std::string line;
            if(std::getline(is, line)){
                std::stringstream sstream(line);
                auto& stream = sstream;

                hq = (data >> 7) & 1;
                useEdits = (data >> 6) & 1;
                type = TempCorrectedSequenceType(int(data & 0x3F));

                if(useEdits){
                    size_t size;
                    stream >> size;
                    int numEdits = size;
                    edits.resize(size);
                    for(int i = 0; i < numEdits; i++){
                        int p;
                        stream >> p;
                        edits[i].pos(p);
                    }
                    for(int i = 0; i < numEdits; i++){
                        char c;
                        stream >> c;
                        edits[i].base(c);
                    }
                }else{
                    stream >> sequence;
                }

                if(type == TempCorrectedSequenceType::Anchor){
                    ; //nothing
                }else{
                    stream >> shift;
                    shift = std::abs(shift);
                }
            }

            return bool(is);
        }


        bool hq = false; //if anchor
        bool useEdits = false;
        TempCorrectedSequenceType type = TempCorrectedSequenceType::Anchor;
        int shift = 0; //if candidate
        read_number readId = 0;

        std::string sequence = "";
        std::vector<CorrectionEdit> edits;


        
    };






}

#endif