cmake_minimum_required(VERSION 3.10)

project(ReRabbitQC)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-DVec256)
set(CMAKE_CXX_FLAGS "-O3 -g -w -ffast-math -funroll-loops -flto -fopenmp  -march=native ${CMAKE_CXX_FLAGS}")
#add_definitions(-DVec512)

#set(CMAKE_CXX_FLAGS " -g  -lz ${CMAKE_CXX_FLAGS}")

aux_source_directory(./src/ DIR_SRCS)
add_executable(ReRabbitQC ${DIR_SRCS})
target_link_libraries(ReRabbitQC z)
target_link_libraries(ReRabbitQC pthread)
