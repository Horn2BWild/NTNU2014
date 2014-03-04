project(ex6)
cmake_minimum_required(VERSION 2.6)

add_subdirectory(../common common)

add_definitions(${COMMON_DEFINITIONS})
list(APPEND CMAKE_C_FLAGS ${COMMON_C_FLAGS})
include_directories(${COMMON_INCLUDE_PATHS})

add_executable(ex6parallel poissonparallel.c)
target_link_libraries(ex6parallel ${COMMON_LIBRARIES})

add_executable(ex6seriell poisson.c)
target_link_libraries(ex6seriell ${COMMON_LIBRARIES})

