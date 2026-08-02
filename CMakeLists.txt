cmake_minimum_required(VERSION 3.18.1)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project (pimdb)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Search the directory tree for files that match the pattern
file(GLOB_RECURSE PIMDB_SRC src/*.cpp)
file(GLOB_RECURSE PIMDB_HEADERS src/*.h)
file(GLOB_RECURSE PIMDB_LIBS libs/*.cpp libs/*.h libs/*.c)

# Include directories
include_directories(include)
include_directories(libs)

# Create the executable
add_executable(pimdb ${PIMDB_SRC} ${PIMDB_HEADERS} ${PIMDB_LIBS})
# Configure inih library for longer configuration lines
target_compile_definitions(pimdb PRIVATE INI_USE_STACK=0 INI_ALLOW_REALLOC=1 INI_MAX_LINE=2048)

# MPI
find_package(MPI REQUIRED)
#include_directories(${FFTW3_INCLUDE_DIR})

include_directories(${MPI_INCLUDE_PATH})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS "${CMAKE_FLAGS} ${MPI_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "GCC detected, adding compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "Clang detected, adding compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message(STATUS "Apple Clang detected, adding compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(STATUS "MSVC detected, adding compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2")
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(STATUS "Intel compiler detected, adding compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    
else()
    message(WARNING "Unknown compiler: ${CMAKE_CXX_COMPILER_ID}")
    # Fallback: try generic flags that might work
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_EXE_LINKER_FLAGS}")

target_link_libraries(pimdb ${MPI_C_LIBRARIES})
#target_link_libraries(pimdb ${MPI_C_LIBRARIES} ${FFTW3_LIBRARIES})

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET pimdb PROPERTY CXX_STANDARD 20)
endif()

# Determine the dimension of space (default = 1)
if (NOT NDIM)
    set(NDIM 3)
endif()
add_definitions(-D NDIM=${NDIM})

# Determine if we want to compile with the naive (factorially-scaling) bosonic algorithm
set(FACTORIAL_BOSONIC_ALGORITHM false CACHE BOOL "Set to true to compile with the naive bosonic algorithm")

if (FACTORIAL_BOSONIC_ALGORITHM)
    add_compile_definitions(FACTORIAL_BOSONIC_ALGORITHM)
endif()
