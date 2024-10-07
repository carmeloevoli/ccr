cmake_minimum_required(VERSION 3.5)
project(CCRH CXX)
set(VERSION 2.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CCRH_EXTRA_SOURCES)
set(CCRH_EXTRA_INCLUDES)
set(CCRH_EXTRA_LIBRARIES)

set(CMAKE_CXX_STANDARD 14)
cmake_policy(SET CMP0048 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-std=c++14 -Wall -pedantic -ftree-vectorize -O3")

# Set default build-type to release to enable performance improvements
if (NOT CMAKE_BUILD_TYPE)
    #set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# C++ Threads required
#find_package(Threads REQUIRED)
#list(APPEND CCRH_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

# GSL (required)
find_package(GSL REQUIRED)
list(APPEND CCRH_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
list(APPEND CCRH_EXTRA_LIBRARIES ${GSL_LIBRARIES})

# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
include_directories(include ${CCRH_EXTRA_INCLUDES})

add_library(CCRH SHARED
    #src/cosmo_progs.cpp
    #src/dump.cpp
    #src/misc.cpp
    #src/ps.cpp
    #src/reionization.cpp
    #src/SFR.cpp
    #src/tridiag.cpp
    src/utilities.cpp
    )
    target_link_libraries(CCRH ${CCRH_EXTRA_LIBRARIES})

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# copy data files in build
#file(GLOB CCRH_DATA "data/*.txt")
#file(COPY ${CCRH_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

# make apps
add_executable(run apps/main.cpp)
target_link_libraries (run CCRH ${CCRH_EXTRA_LIBRARIES})

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "CCRH_INCLUDES ...... = ${CCRH_EXTRA_INCLUDES}")
message (STATUS "CCRH_LIBRARIES ..... = ${CCRH_EXTRA_LIBRARIES}")

