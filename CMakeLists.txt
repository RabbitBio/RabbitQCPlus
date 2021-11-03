cmake_minimum_required(VERSION 3.10)

#set(CMAKE_C_COMPILER "")


project(ReRabbitQC)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-DVec512)
set(CMAKE_CXX_FLAGS "-I../ -I../common -O3 -g -w -ffast-math -funroll-loops -flto=jobserver   -fopenmp  -march=native -mtune=native -D_POSIX_C_SOURCE=200809L -D_FILE_OFFSET_BITS=64 ${CMAKE_CXX_FLAGS}")

#set(CMAKE_CXX_FLAGS " -g  -lz ${CMAKE_CXX_FLAGS}")

aux_source_directory(./src/ DIR_SRCS)
add_executable(ReRabbitQC ${DIR_SRCS})
target_link_libraries(ReRabbitQC z)
target_link_libraries(ReRabbitQC pthread)
target_link_libraries(ReRabbitQC rt)