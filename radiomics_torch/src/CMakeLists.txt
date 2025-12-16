Python_add_library(_cmatrices MODULE _cmatrices.c cmatrices.c)
target_link_libraries(_cmatrices PRIVATE Python::NumPy)

Python_add_library(_cshape MODULE _cshape.c cshape.c)
target_link_libraries(_cshape PRIVATE Python::NumPy)

install(TARGETS _cmatrices _cshape DESTINATION radiomics)
