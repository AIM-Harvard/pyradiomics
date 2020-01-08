cmake_minimum_required(VERSION 3.7)

project(radiomics)

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(PythonExtensions REQUIRED)

add_subdirectory(radiomics/src)
