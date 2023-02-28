cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

set(EXAMPLES_BINARY_DIR "${CMAKE_BINARY_DIR}/examples")

###################################################################################################
# - example targets -------------------------------------------------------------------------------

add_executable(helpers_example "${CMAKE_CURRENT_SOURCE_DIR}/helpers_example.cu")
set_target_properties(helpers_example PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXAMPLES_BINARY_DIR})
target_link_libraries(helpers_example PRIVATE helpers)

add_executable(packed_types_example "${CMAKE_CURRENT_SOURCE_DIR}/packed_types_example.cu")
set_target_properties(packed_types_example PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXAMPLES_BINARY_DIR})
target_link_libraries(packed_types_example PRIVATE helpers)