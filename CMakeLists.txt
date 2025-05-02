cmake_minimum_required(VERSION 3.24)
project(master_algorithm)

set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/capdMake/cmake/")

include_directories("/home/agnieszka/CAPD/capdDynSys/include")
include_directories("/home/agnieszka/CAPD/capdAlg/include")
include_directories("/home/agnieszka/CAPD/capdAux/include")
include_directories("/home/agnieszka/CAPD/capdExt/include")
include_directories("/home/agnieszka/CAPD/capdDynSys/lib")
include_directories("/home/agnieszka/CAPD/capdAlg/lib")
include_directories("/home/agnieszka/CAPD/capdAux/lib")
include_directories("/home/agnieszka/CAPD/capdExt/lib")

find_package(CAPD REQUIRED PATHS /home/agnieszka/CAPD/build)

add_executable(master_algorithm main.cpp
        source/ChebyshevSeries/ChebyshevSeries.hpp
        source/ChebyshevSeries/ChebyshevSeries.tpp
        source/ChebyshevSeries/norm.hpp
        source/ChebyshevSeries/norm.tpp
        source/ChebyshevSeries/ChebyshevOperatorFinite.hpp
        source/ChebyshevSeries/ChebyshevOperatorFinite.tpp source/ChebyshevSeries/ChebyshevOperatorInfinite.hpp source/ChebyshevSeries/ChebyshevOperatorInfinite.tpp)

target_link_libraries(master_algorithm
        /home/agnieszka/CAPD/build/libcapd.a
        /home/agnieszka/CAPD/build/capdExt/filibsrc/libfilib.a
        )
