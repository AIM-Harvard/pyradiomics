find_package(NumPy REQUIRED)

add_library(_cmatrices MODULE _cmatrices.c cmatrices.c)
python_extension_module(_cmatrices)
target_include_directories(_cmatrices PRIVATE ${NumPy_INCLUDE_DIR})

add_library(_cshape MODULE _cshape.c cshape.c)
python_extension_module(_cshape)
target_include_directories(_cshape PRIVATE ${NumPy_INCLUDE_DIR})

install(TARGETS _cmatrices _cshape LIBRARY DESTINATION radiomics)
