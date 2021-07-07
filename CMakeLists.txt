cmake_minimum_required(VERSION 3.18)
project(RabbitQCPlus)
#set(CMAKE_CXX_COMPILER "/usr/local/Cellar/gcc@8/8.4.0_1/bin/gcc-8")

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-O3 -g -ffast-math -lz ${CMAKE_CXX_FLAGS}")

aux_source_directory(./src/ DIR_SRCS)
add_executable(RabbitQCPlus ${DIR_SRCS})
#target_link_libraries(RabbitQCPlus z)