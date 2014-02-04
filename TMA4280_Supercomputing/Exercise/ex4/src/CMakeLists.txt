project(ex4)
cmake_minimum_required(VERSION 2.6)

add_subdirectory(../common common)

add_definitions(${COMMON_DEFINITIONS})
list(APPEND CMAKE_C_FLAGS ${COMMON_C_FLAGS})
include_directories(${COMMON_INCLUDE_PATHS})

add_executable(ex4 ex4.c)
target_link_libraries(ex4 ${COMMON_LIBRARIES})

