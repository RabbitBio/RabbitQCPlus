#ifndef CARE_DISPATCH_CARE_CORRECT_GPU_CUH
#define CARE_DISPATCH_CARE_CORRECT_GPU_CUH

#include <options.hpp>

#include <vector>

namespace care{

    namespace correction{

        std::vector<int> getUsableDeviceIds(std::vector<int> deviceIds);

    }

    void performCorrection(
        ProgramOptions programOptions
    );

}




#endif