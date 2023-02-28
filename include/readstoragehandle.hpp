#ifndef CARE_READSTORAGEHANDLE_HPP
#define CARE_READSTORAGEHANDLE_HPP

#include <limits>

namespace care{

    class CpuReadStorage; //forward declaration
    namespace gpu{
        class GpuReadStorage; //forward declaration
    }

    class ReadStorageHandle{

        friend class CpuReadStorage;
        friend class gpu::GpuReadStorage;

    public:

        ReadStorageHandle() = default;
        ReadStorageHandle(int i) : id(i){}

        int getId() const noexcept{
            return id;
        }

    private:

        int id = std::numeric_limits<int>::max();
    };


}


#endif