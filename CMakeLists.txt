cmake_minimum_required(VERSION 3.24)
project(master_algorithm)

set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/capdMake/cmake/")

include_directories(
        "/home/agnieszka/CAPD/capdDynSys/include"
        "/home/agnieszka/CAPD/capdAlg/include"
        "/home/agnieszka/CAPD/capdAux/include"
        "/home/agnieszka/CAPD/capdExt/include"
        "/home/agnieszka/CAPD/capdMake/libcapd/include"
        "/home/agnieszka/CAPD/capdMake/bin/include"
)

link_directories(
        "/home/agnieszka/CAPD/capdDynSys/lib"
        "/home/agnieszka/CAPD/capdAlg/lib"
        "/home/agnieszka/CAPD/capdAux/lib"
        "/home/agnieszka/CAPD/capdExt/lib"
        "/home/agnieszka/CAPD/capdMake/libcapd/lib"
        "/home/agnieszka/CAPD/capdMake/bin/lib"
)


find_package(CAPD REQUIRED PATHS /home/agnieszka/CAPD/build)

add_executable(master_algorithm main.cpp
        source/ChebyshevSeries/ChebyshevSeries.hpp
        source/ChebyshevSeries/ChebyshevSeries.tpp
        source/ChebyshevSeries/Norm.hpp
        source/ChebyshevSeries/norm.tpp
        source/ChebyshevSeries/ChebyshevOperatorFinite.hpp
        source/ChebyshevSeries/ChebyshevOperatorFinite.tpp source/ChebyshevSeries/RadiiPolynomials.hpp source/ChebyshevSeries/RadiiPolynomials.tpp helpers.hpp testing.hpp)

target_link_libraries(master_algorithm
        /home/agnieszka/CAPD/build/libcapd.a
        /home/agnieszka/CAPD/build/capdExt/filibsrc/libfilib.a
        )
