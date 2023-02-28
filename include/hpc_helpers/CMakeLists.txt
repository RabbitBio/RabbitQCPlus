cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(hpc_helpers VERSION 0.0.1 LANGUAGES CXX CUDA)

###################################################################################################
# - build type ------------------------------------------------------------------------------------

# Set a default build type if none was specified
set(DEFAULT_BUILD_TYPE "Release")

###################################################################################################
# - find packages we depend on --------------------------------------------------------------------

find_package(CUDAToolkit 9.0 REQUIRED)

###################################################################################################
# - define target   -------------------------------------------------------------------------------
add_library(helpers INTERFACE)
target_include_directories(helpers INTERFACE
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                           $<INSTALL_INTERFACE:include>)

###################################################################################################
# - build options ---------------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")

option(HELPERS_BUILD_EXAMPLES "Configure CMake to build examples" OFF)
if(HELPERS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif(HELPERS_BUILD_EXAMPLES)

option(HELPERS_BUILD_TESTS "Configure CMake to build tests" OFF)
if(HELPERS_BUILD_TESTS)
    add_subdirectory(tests)
endif(HELPERS_BUILD_TESTS)