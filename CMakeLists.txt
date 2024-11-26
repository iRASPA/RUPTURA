cmake_minimum_required(VERSION 3.14)

project(ruptura LANGUAGES CXX)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(src)

# Compiler flags
set(CXX_COMPILE_FLAGS
    -g
    -O3
    -march=native
    -ffast-math
    -Wall
    -Wextra
    -Wshadow
    -Wnon-virtual-dtor
    -Wold-style-cast
    -Wcast-align
    -Wunused
    -Woverloaded-virtual
    -Wpedantic
    -Wconversion
    -Wsign-conversion
    -Wnull-dereference
    -Wdouble-promotion
    -Wformat=2
    -Werror
    -fomit-frame-pointer
    -ftree-vectorize
    -fno-stack-check
    -funroll-loops
)

# Collect source files
set(SOURCES
    src/breakthrough.cpp
    src/component.cpp
    src/fitting.cpp
    src/inputreader.cpp
    src/isotherm.cpp
    src/mixture_prediction.cpp
    src/multi_site_isotherm.cpp
    src/random_numbers.cpp
    src/special_functions.cpp
)

# -------------------------------
# Build the ruptura executable
# -------------------------------
add_executable(ruptura src/main.cpp ${SOURCES})

target_compile_options(ruptura PRIVATE ${CXX_COMPILE_FLAGS})

# -------------------------------
# Doxygen Documentation
# -------------------------------
add_subdirectory(docs)
