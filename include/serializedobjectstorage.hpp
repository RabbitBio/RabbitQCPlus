#ifndef CARE_SERIALIZED_OBJECT_STORAGE_HPP
#define CARE_SERIALIZED_OBJECT_STORAGE_HPP

#include <mmapbuffer.hpp>
#include <filehelpers.hpp>
#include <memorymanagement.hpp>

#include <iostream>
#include <algorithm>
#include <memory>
namespace care{


class SerializedObjectStorage{
private:
    std::unique_ptr<FileBackedUVector<std::uint8_t>> databuffer;
    std::unique_ptr<FileBackedUVector<std::size_t>> offsetbuffer;
public:
    SerializedObjectStorage(std::size_t memoryLimitData_, std::size_t memoryLimitOffsets_, std::string filedirectory){
        std::string nametemplate = filedirectory + "serializedobjectstorage-XXXXXX";
        std::string databufferfilename = filehelpers::makeRandomFile(nametemplate);
        std::string offsetbufferfilename = filehelpers::makeRandomFile(nametemplate);

        databuffer = std::make_unique<FileBackedUVector<std::uint8_t>>(
            0, memoryLimitData_, databufferfilename);

        offsetbuffer = std::make_unique<FileBackedUVector<std::size_t>>(
            0, memoryLimitOffsets_, offsetbufferfilename);
    }

    void insert(const std::uint8_t* begin, const std::uint8_t* end){
        auto first = databuffer->insert(databuffer->end(), begin, end);
        offsetbuffer->push_back(std::distance(databuffer->begin(), first));
    }

    MemoryUsage getMemoryInfo() const{
        MemoryUsage result;
        result.host += databuffer->getCapacityInMemoryInBytes();
        result.host += offsetbuffer->getCapacityInMemoryInBytes();
        return result;
    }

    std::uint8_t* getPointer(std::size_t i) noexcept{
        const std::size_t offset = getOffset(i);
        return databuffer->data() + offset;
    }

    const std::uint8_t* getPointer(std::size_t i) const noexcept{
        const std::size_t offset = getOffset(i);
        return databuffer->data() + offset;
    }

    std::size_t getOffset(std::size_t i) const noexcept{
        return (*offsetbuffer)[i];
    }

    std::uint8_t* getDataBuffer() noexcept{
        return databuffer->data();
    }

    const std::uint8_t* getDataBuffer() const noexcept{
        return databuffer->data();
    }

    std::size_t* getOffsetBuffer() noexcept{
        return offsetbuffer->data();
    }

    const std::size_t* getOffsetBuffer() const noexcept{
        return offsetbuffer->data();
    }

    std::size_t size() const noexcept{
        return offsetbuffer->size();
    }

    std::size_t getNumElements() const noexcept{
        return offsetbuffer->size();
    }    

    std::size_t dataBytes() const noexcept{
        return sizeof(std::uint8_t) * databuffer->size();
    }

    std::size_t offsetBytes() const noexcept{
        return sizeof(std::size_t) * offsetbuffer->size();
    }


    void saveToStream(std::ostream& out) const{
        std::size_t dbytes = dataBytes();
        std::size_t obytes = offsetBytes();
        std::size_t delemns = databuffer->size();
        std::size_t oelemns = offsetbuffer->size();

        out.write(reinterpret_cast<const char*>(&dbytes), sizeof(std::size_t));
        out.write(reinterpret_cast<const char*>(&obytes), sizeof(std::size_t));
        out.write(reinterpret_cast<const char*>(&delemns), sizeof(std::size_t));
        out.write(reinterpret_cast<const char*>(&oelemns), sizeof(std::size_t));

        out.write(reinterpret_cast<const char*>(getDataBuffer()), dbytes);
        out.write(reinterpret_cast<const char*>(getOffsetBuffer()), obytes);
    }

    void loadFromStream(std::istream& in){
        std::size_t dbytes = 0;
        std::size_t obytes = 0;
        std::size_t delemns = 0;
        std::size_t oelemns = 0;

        in.read(reinterpret_cast<char*>(&dbytes), sizeof(std::size_t));
        in.read(reinterpret_cast<char*>(&obytes), sizeof(std::size_t));
        in.read(reinterpret_cast<char*>(&delemns), sizeof(std::size_t));
        in.read(reinterpret_cast<char*>(&oelemns), sizeof(std::size_t));

        databuffer->resize(delemns);
        offsetbuffer->resize(oelemns);

        in.read(reinterpret_cast<char*>(getDataBuffer()), dbytes);
        in.read(reinterpret_cast<char*>(getOffsetBuffer()), obytes);
    }
};


}





#endif