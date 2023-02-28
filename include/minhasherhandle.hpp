#ifndef CARE_MINHASHERHANDLE_HPP
#define CARE_MINHASHERHANDLE_HPP

#include <limits>

namespace care{

    class CpuMinhasher; //forward declaration
    namespace gpu{
        class GpuMinhasher; //forward declaration
    }

    class MinhasherHandle{

        friend class CpuMinhasher;
        friend class gpu::GpuMinhasher;

    public:

        MinhasherHandle() = default;
        MinhasherHandle(int i) : id(i){}

        int getId() const noexcept{
            return id;
        }

    private:

        int id = std::numeric_limits<int>::max();
    };


}


#endif