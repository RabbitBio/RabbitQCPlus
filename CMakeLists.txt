cmake_minimum_required(VERSION 3.10)

project(RabbitQCPlus)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-O3 -g -w -ffast-math ${CMAKE_CXX_FLAGS}")
#add_definitions(-DVec512)

#set(CMAKE_CXX_FLAGS " -g  -lz ${CMAKE_CXX_FLAGS}")

aux_source_directory(./src/ DIR_SRCS)
add_executable(RabbitQCPlus ${DIR_SRCS})
target_link_libraries(RabbitQCPlus z)
target_link_libraries(RabbitQCPlus pthread)
