project(ex3)
cmake_minimum_required(VERSION 2.6)

enable_language(C)

option(ENABLE_OPENMP "Enable OpenMP support?" ON)
option(ENABLE_MPI    "Enable MPI support?" ON)

list(APPEND COMMON_LIBRARIES common)

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    list(APPEND COMMON_DEFINITIONS -DHAVE_OPENMP)
    list(APPEND COMMON_C_FLAGS ${OpenMP_C_FLAGS})
  endif()
endif()

if(ENABLE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    list(APPEND COMMON_DEFINITIONS -DHAVE_MPI)
    list(APPEND COMMON_C_FLAGS ${MPI_C_FLAGS})
    list(APPEND COMMON_LIBRARIES ${MPI_C_LIBRARIES})
    list(APPEND COMMON_INCLUDE_PATHS ${MPI_C_INCLUDE_PATH})
  endif()
endif()

list(APPEND CMAKE_C_FLAGS ${COMMON_C_FLAGS})
add_definitions(${COMMON_DEFINITIONS})
include_directories(${COMMON_INCLUDE_PATHS})

list(APPEND COMMON_INCLUDE_PATHS ${PROJECT_SOURCE_DIR})


add_executable(ex3 ex3.c)
