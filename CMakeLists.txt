cmake_minimum_required(VERSION 3.15...3.26)

project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C)

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

add_subdirectory(radiomics/src)
