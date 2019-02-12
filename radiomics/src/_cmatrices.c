#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <stdlib.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "cmatrices.h"

static char module_docstring[] = ("This module links to C-compiled code for efficient calculation of various matrices "
                                 "in the pyRadiomics package. It provides fast calculation for GLCM, GLDM, NGTDM, "
                                 "GLRLM and GLSZM. All functions are given names as follows: ""calculate_<Matrix>"", "
                                 "where <Matrix> is the name of the matix, in lowercase. Arguments for these functions "
                                 "are positional and start with 2-3 numpy arrays (image, mask and [distances]) and 3 integers "
                                 "(Ng number of gray levels, force2D and force2Ddimension). Optionally extra arguments "
                                 "may be required, see function docstrings for detailed information. "
                                 "Functions return a tuple with the calculated matrix and angles.");
static char glcm_docstring[] = "Arguments: Image, Mask, Ng, force2D, for2Ddimension.";
static char glszm_docstring[] = "Arguments: Image, Mask, Ng, Ns, force2D, for2Ddimension, matrix is cropped to maximum size encountered.";
static char glrlm_docstring[] = "Arguments: Image, Mask, Ng, Nr, force2D, for2Ddimension.";
static char ngtdm_docstring[] = "Arguments: Image, Mask, Ng, force2D, for2Ddimension.";
static char gldm_docstring[] = "Arguments: Image, Mask, Ng, Alpha, force2D, for2Ddimension.";
static char generate_angles_docstring[] = "Arguments: Boundingbox Size, distances, bidirectional, force2Ddimension.";

static PyObject *cmatrices_calculate_glcm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_glszm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_glrlm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_ngtdm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_gldm(PyObject *self, PyObject *args);
static PyObject *cmatrices_generate_angles(PyObject *self, PyObject *args);

// Function to check if array input is valid. Additionally extracts size and stride values
int build_angles_arr(PyObject *distances_obj, PyArrayObject **angles_arr, int *size, int Nd, int force2Ddimension, int bidirectional, int *Na);
int try_parse_arrays(PyObject *image_obj, PyObject *mask_obj, PyArrayObject **image_arr, PyArrayObject **mask_arr, int *Nd, int **size, int **strides, int mask_flags);

static PyMethodDef module_methods[] = {
  //{"calculate_", cmatrices_, METH_VARARGS, _docstring},
  { "calculate_glcm", cmatrices_calculate_glcm, METH_VARARGS, glcm_docstring },
  { "calculate_glszm", cmatrices_calculate_glszm, METH_VARARGS, glszm_docstring },
  { "calculate_glrlm", cmatrices_calculate_glrlm, METH_VARARGS, glrlm_docstring },
  { "calculate_ngtdm", cmatrices_calculate_ngtdm, METH_VARARGS, ngtdm_docstring },
  { "calculate_gldm", cmatrices_calculate_gldm, METH_VARARGS, gldm_docstring },
  { "generate_angles", cmatrices_generate_angles, METH_VARARGS, generate_angles_docstring},
  { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_cmatrices",        /* m_name */
  module_docstring,    /* m_doc */
  -1,                  /* m_size */
  module_methods,      /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

#endif

static PyObject *
moduleinit(void)
{
    PyObject *m;

#if PY_MAJOR_VERSION >= 3
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule3("_cmatrices",
                       module_methods, module_docstring);
#endif

  if (!m)
      return NULL;

  return m;
}

#if PY_MAJOR_VERSION < 3
  PyMODINIT_FUNC
  init_cmatrices(void)
  {
    // Initialize numpy functionality
    import_array();

    moduleinit();

  }
#else
  PyMODINIT_FUNC
  PyInit__cmatrices(void)
  {
    // Initialize numpy functionality
    import_array();

    return moduleinit();
  }
#endif

static PyObject *cmatrices_calculate_glcm(PyObject *self, PyObject *args)
{
  int Ng, force2D, force2Ddimension;
  PyObject *image_obj, *mask_obj, *distances_obj;
  PyArrayObject *image_arr, *mask_arr;
  int Nd, Na;
  int *size, *bb, *strides;
  npy_intp dims[3];
  PyArrayObject *glcm_arr, *angles_arr;
  int *image;
  char *mask;
  int *angles;
  double *glcm;
  int d;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOOiii", &image_obj, &mask_obj, &distances_obj, &Ng, &force2D, &force2Ddimension))
    return NULL;

  // Check if array input is valid and extract sizes and strides of image and mask
  // Returns 0 if successful, -1 if failed.
  if(try_parse_arrays(image_obj, mask_obj, &image_arr, &mask_arr, &Nd, &size, &strides, 0) < 0)
    return NULL;

  // If extraction is not forced 2D, ensure the dimension is set to a non-existent one (ensuring 3D angles when possible)
  if(!force2D) force2Ddimension = -1;

  if(build_angles_arr(distances_obj, &angles_arr, size, Nd, force2Ddimension, 0, &Na) < 0)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);

    free(size);
    free(strides);
    return NULL;
  }

  angles = (int *)PyArray_DATA(angles_arr);

  // Initialize output array (elements not set)
  dims[0] = Ng;
  dims[1] = Ng;
  dims[2] = Na;

  // Check that the maximum size of the array won't overflow the index variable (int32)
  if (dims[0] * dims[1] * dims[2] > INT_MAX)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    PyErr_SetString(PyExc_RuntimeError, "Number of elements in GLCM would overflow index variable! Increase bin width or decrease number of angles to prevent this error.");
    return NULL;
  }

  glcm_arr = (PyArrayObject *)PyArray_SimpleNew(3, dims, NPY_DOUBLE);
  if (!glcm_arr)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    return PyErr_NoMemory();
  }

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (char *)PyArray_DATA(mask_arr);
  glcm = (double *)PyArray_DATA(glcm_arr);

  // Set all elements to 0
  memset(glcm, 0, sizeof *glcm * Ng * Ng * Na);

  // initialize bb
  bb = (int *)malloc(sizeof *bb * Nd * 2);
  if (!bb)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(glcm_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    return PyErr_NoMemory();
  }
  memset(bb, 0, sizeof *bb * Nd);
  for (d = 0; d < Nd; d++)
    bb[Nd + d] = size[d] - 1;

  //Calculate GLCM
  if (!calculate_glcm(image, mask, size, bb, strides, angles, Na, Nd, glcm, Ng))
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(glcm_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);
    free(bb);

    PyErr_SetString(PyExc_IndexError, "Calculation of GLCM Failed.");
    return NULL;
  }

  // Clean up
  Py_XDECREF(image_arr);
  Py_XDECREF(mask_arr);

  free(size);
  free(strides);
  free(bb);

  return Py_BuildValue("NN", PyArray_Return(glcm_arr), PyArray_Return(angles_arr));
}

static PyObject *cmatrices_calculate_glszm(PyObject *self, PyObject *args)
{
  int Ng, Ns, force2D, force2Ddimension;
  PyObject *image_obj, *mask_obj;
  PyArrayObject *image_arr, *mask_arr;
  int Na, Nd;
  int *size, *bb, *strides;
  npy_intp dims[2];
  PyArrayObject *glszm_arr, *angles_arr;
  int *image;
  char *mask;
  int *angles;
  int *tempData;
  int maxRegion;
  double *glszm;
  int d;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOiiii", &image_obj, &mask_obj, &Ng, &Ns, &force2D, &force2Ddimension))
    return NULL;

  // Check if array input is valid and extract sizes and strides of image and mask
  // Returns 0 if successful, 1-3 if failed.
  if(try_parse_arrays(image_obj, mask_obj, &image_arr, &mask_arr, &Nd, &size, &strides, NPY_ARRAY_ENSURECOPY) < 0)
    return NULL;

  // If extraction is not forced 2D, ensure the dimension is set to a non-existent one (ensuring 3D angles when possible)
  if(!force2D) force2Ddimension = -1;

  if(build_angles_arr(NULL, &angles_arr, size, Nd, force2Ddimension, 1, &Na) < 0)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);

    free(size);
    free(strides);
    return NULL;
  }

  angles = (int *)PyArray_DATA(angles_arr);

  // Initialize temporary output array (elements not set)
  // add +1 to the size so in the case every voxel represents a separate region,
  // tempData still contains a -1 element at the end
  tempData = (int *)malloc(sizeof *tempData * (2 * Ns + 1));
  if (!tempData)  // No memory allocated
  {
	  Py_XDECREF(image_arr);
	  Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

	  free(size);
    free(strides);

    return PyErr_NoMemory();
  }

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (char *)PyArray_DATA(mask_arr);

  // initialize bb
  bb = (int *)malloc(sizeof *bb * Nd * 2);
  if (!bb)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(glszm_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    return PyErr_NoMemory();
  }
  memset(bb, 0, sizeof *bb * Nd);
  for (d = 0; d < Nd; d++)
    bb[Nd + d] = size[d] - 1;

  //Calculate GLSZM
  maxRegion = calculate_glszm(image, mask, size, bb, strides, angles, Na, Nd, tempData, Ng, Ns);
  if (maxRegion < 0) // Error occured
  {
	  free(tempData);
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);
    free(bb);

    PyErr_SetString(PyExc_IndexError, "Calculation of GLSZM Failed.");
    return NULL;
  }

  // Clean up image, mask and angles arrays (not needed anymore)
  Py_XDECREF(image_arr);
  Py_XDECREF(mask_arr);
  Py_XDECREF(angles_arr);

  free(size);
  free(strides);
  free(bb);

  // Initialize output array (elements not set)
  if (maxRegion == 0) maxRegion = 1;
  dims[0] = Ng;
  dims[1] = maxRegion;

  // Check that the maximum size of the array won't overflow the index variable (int32)
  if (dims[0] * dims[1] > INT_MAX)
  {
    free(tempData);
    PyErr_SetString(PyExc_RuntimeError, "Number of elements in GLSZM would overflow index variable! Increase bin width to prevent this error.");
    return NULL;
  }

  glszm_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  if (!glszm_arr)
  {
    free(tempData);

    return PyErr_NoMemory();
  }

  glszm = (double *)PyArray_DATA(glszm_arr);

  // Set all elements to 0
  memset(glszm, 0, sizeof *glszm * maxRegion * Ng);

  if (!fill_glszm(tempData, glszm, Ng, maxRegion))
  {
    free(tempData);
    Py_XDECREF(glszm_arr);
    PyErr_SetString(PyExc_IndexError, "Error filling GLSZM.");
    return NULL;
  }

  // Clean up
  free(tempData);

  return PyArray_Return(glszm_arr);
}

static PyObject *cmatrices_calculate_glrlm(PyObject *self, PyObject *args)
{
  int Ng, Nr, force2D, force2Ddimension;
  PyObject *image_obj, *mask_obj;
  PyArrayObject *image_arr, *mask_arr;
  int Nd, Na;
  int *size, *bb, *strides;
  npy_intp dims[3];
  PyArrayObject *glrlm_arr, *angles_arr;
  int *image;
  char *mask;
  int *angles;
  double *glrlm;
  int d;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOiiii", &image_obj, &mask_obj, &Ng, &Nr, &force2D, &force2Ddimension))
    return NULL;

  // Check if array input is valid and extract sizes and strides of image and mask
  // Returns 0 if successful, -1 if failed.
  if(try_parse_arrays(image_obj, mask_obj, &image_arr, &mask_arr, &Nd, &size, &strides, 0) < 0)
    return NULL;

  // If extraction is not forced 2D, ensure the dimension is set to a non-existent one (ensuring 3D angles when possible)
  if(!force2D) force2Ddimension = -1;

  if(build_angles_arr(NULL, &angles_arr, size, Nd, force2Ddimension, 0, &Na) < 0)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);

    free(size);
    free(strides);
    return NULL;
  }

  angles = (int *)PyArray_DATA(angles_arr);

  // Initialize output array (elements not set)
  dims[0] = Ng;
  dims[1] = Nr;
  dims[2] = Na;

  // Check that the maximum size of the array won't overflow the index variable (int32)
  if (dims[0] * dims[1] * dims[2] > INT_MAX)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    PyErr_SetString(PyExc_RuntimeError, "Number of elements in GLRLM would overflow index variable! Increase bin width or decrease number of angles to prevent this error.");
    return NULL;
  }

  glrlm_arr = (PyArrayObject *)PyArray_SimpleNew(3, dims, NPY_DOUBLE);
  if (!glrlm_arr)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    return PyErr_NoMemory();
  }

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (char *)PyArray_DATA(mask_arr);
  glrlm = (double *)PyArray_DATA(glrlm_arr);

  // Set all elements to 0
  memset(glrlm, 0, sizeof *glrlm * Ng * Nr * Na);

  // initialize bb
  bb = (int *)malloc(sizeof *bb * Nd * 2);
  if (!bb)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(glrlm_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    return PyErr_NoMemory();
  }
  memset(bb, 0, sizeof *bb * Nd);
  for (d = 0; d < Nd; d++)
    bb[Nd + d] = size[d] - 1;

  //Calculate GLRLM
  if (!calculate_glrlm(image, mask, size, bb, strides, angles, Na, Nd, glrlm, Ng, Nr))
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(glrlm_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);
    free(bb);

    PyErr_SetString(PyExc_IndexError, "Calculation of GLRLM Failed.");
    return NULL;
  }

  // Clean up
  Py_XDECREF(image_arr);
  Py_XDECREF(mask_arr);

  free(size);
  free(strides);
  free(bb);

  return Py_BuildValue("NN", PyArray_Return(glrlm_arr), PyArray_Return(angles_arr));
}

static PyObject *cmatrices_calculate_ngtdm(PyObject *self, PyObject *args)
{
  int Ng, force2D, force2Ddimension;
  PyObject *image_obj, *mask_obj, *distances_obj;
  PyArrayObject *image_arr, *mask_arr;
  int Nd, Na;
  int *size, *bb, *strides;
  npy_intp dims[2];
  PyArrayObject *ngtdm_arr, *angles_arr;
  int *image;
  char *mask;
  int *angles;
  double *ngtdm;
  int d;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOOiii", &image_obj, &mask_obj, &distances_obj, &Ng, &force2D, &force2Ddimension))
    return NULL;

  // Check if array input is valid and extract sizes and strides of image and mask
  // Returns 0 if successful, -1 if failed.
  if(try_parse_arrays(image_obj, mask_obj, &image_arr, &mask_arr, &Nd, &size, &strides, 0) < 0)
    return NULL;

  // If extraction is not forced 2D, ensure the dimension is set to a non-existent one (ensuring 3D angles when possible)
  if(!force2D) force2Ddimension = -1;

  if(build_angles_arr(distances_obj, &angles_arr, size, Nd, force2Ddimension, 1, &Na) < 0)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);

    free(size);
    free(strides);
    return NULL;
  }

  angles = (int *)PyArray_DATA(angles_arr);

  // Initialize output array (elements not set)
  dims[0] = Ng;
  dims[1] = 3;

  // Check that the maximum size of the array won't overflow the index variable (int32)
  if (dims[0] * dims[1] > INT_MAX)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    PyErr_SetString(PyExc_RuntimeError, "Number of elements in NGTDM would overflow index variable! Increase bin width to prevent this error.");
    return NULL;
  }

  ngtdm_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  if (!ngtdm_arr)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    return PyErr_NoMemory();
  }

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (char *)PyArray_DATA(mask_arr);
  ngtdm = (double *)PyArray_DATA(ngtdm_arr);

  // Set all elements to 0
  memset(ngtdm, 0, sizeof *ngtdm * Ng * 3);

  // initialize bb
  bb = (int *)malloc(sizeof *bb * Nd * 2);
  if (!bb)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(ngtdm_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    return PyErr_NoMemory();
  }
  memset(bb, 0, sizeof *bb * Nd);
  for (d = 0; d < Nd; d++)
    bb[Nd + d] = size[d] - 1;

  //Calculate NGTDM
  if (!calculate_ngtdm(image, mask, size, bb, strides, angles, Na, Nd, ngtdm, Ng))
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(ngtdm_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);
    free(bb);

    PyErr_SetString(PyExc_IndexError, "Calculation of NGTDM Failed.");
    return NULL;
  }

  // Clean up
  Py_XDECREF(image_arr);
  Py_XDECREF(mask_arr);
  Py_XDECREF(angles_arr);

  free(size);
  free(strides);
  free(bb);

  return PyArray_Return(ngtdm_arr);
}

static PyObject *cmatrices_calculate_gldm(PyObject *self, PyObject *args)
{
  int Ng, alpha, force2D, force2Ddimension;
  PyObject *image_obj, *mask_obj, *distances_obj;
  PyArrayObject *image_arr, *mask_arr;
  int Nd, Na;
  int *size, *bb, *strides;
  npy_intp dims[2];
  PyArrayObject *gldm_arr, *angles_arr;
  int *image;
  char *mask;
  int *angles;
  double *gldm;
  int d;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOOiiii", &image_obj, &mask_obj, &distances_obj, &Ng, &alpha, &force2D, &force2Ddimension))
    return NULL;

  // Check if array input is valid and extract sizes and strides of image and mask
  // Returns 0 if successful, -1 if failed.
  if(try_parse_arrays(image_obj, mask_obj, &image_arr, &mask_arr, &Nd, &size, &strides, 0) < 0)
    return NULL;

  // If extraction is not forced 2D, ensure the dimension is set to a non-existent one (ensuring 3D angles when possible)
  if(!force2D) force2Ddimension = -1;

  if(build_angles_arr(distances_obj, &angles_arr, size, Nd, force2Ddimension, 1, &Na) < 0)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);

    free(size);
    free(strides);
    return NULL;
  }

  angles = (int *)PyArray_DATA(angles_arr);

  // Initialize output array (elements not set)
  dims[0] = Ng;
  dims[1] = Na * 2 + 1;  // No of possible dependency values = Na *2 + 1 (Na angels, 2 directions and +1 for no dependency)

  // Check that the maximum size of the array won't overflow the index variable (int32)
  if (dims[0] * dims[1] > INT_MAX)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    PyErr_SetString(PyExc_RuntimeError, "Number of elements in GLDM would overflow index variable! Increase bin width or decrease number of angles to prevent this error.");
    return NULL;
  }

  gldm_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  if (!gldm_arr)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    return PyErr_NoMemory();
  }

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (char *)PyArray_DATA(mask_arr);
  gldm = (double *)PyArray_DATA(gldm_arr);

  // Set all elements to 0
  memset(gldm, 0, sizeof *gldm * Ng * (Na * 2 + 1));

  // initialize bb
  bb = (int *)malloc(sizeof *bb * Nd * 2);
  if (!bb)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(gldm_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);

    return PyErr_NoMemory();
  }
  memset(bb, 0, sizeof *bb * Nd);
  for (d = 0; d < Nd; d++)
    bb[Nd + d] = size[d] - 1;

  //Calculate GLDM
  if (!calculate_gldm(image, mask, size, bb, strides, angles, Na, Nd, gldm, Ng, alpha))
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(gldm_arr);
    Py_XDECREF(angles_arr);

    free(size);
    free(strides);
    free(bb);

    PyErr_SetString(PyExc_IndexError, "Calculation of GLDM Failed.");
    return NULL;
  }

  // Clean up
  Py_XDECREF(image_arr);
  Py_XDECREF(mask_arr);
  Py_XDECREF(angles_arr);

  free(size);
  free(strides);
  free(bb);

  return PyArray_Return(gldm_arr);
}

static PyObject *cmatrices_generate_angles(PyObject *self, PyObject *args)
{
  PyObject *size_obj, *distances_obj;
  PyArrayObject *size_arr;
  char bidirectional;
  int force2D, force2Ddimension;
  int Nd, Na;
  int *size;
  npy_intp dims[2];
  PyArrayObject *angles_arr;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOiii", &size_obj, &distances_obj, &bidirectional, &force2D, &force2Ddimension))
    return NULL;

  // Interpret the size input as numpy array
  size_arr = (PyArrayObject *)PyArray_FROM_OTF(size_obj, NPY_INT, NPY_ARRAY_FORCECAST | NPY_ARRAY_IN_ARRAY);

  if (!size_arr)
    return NULL;

  // Check if Size has 1 dimension
  if (PyArray_NDIM(size_arr) != 1)
  {
    Py_XDECREF(size_arr);
    PyErr_SetString(PyExc_ValueError, "Expected a 1D array for size");
    return NULL;
  }

  Nd = (int)PyArray_DIM(size_arr, 0);  // Number of dimensions
  size = (int *)PyArray_DATA(size_arr);

  if (!force2D) force2Ddimension = -1;

  if(build_angles_arr(distances_obj, &angles_arr, size, Nd, force2Ddimension, bidirectional, &Na) < 0)
  {
    Py_XDECREF(size_arr);
    return NULL;
  }

  Py_XDECREF(size_arr);

  return PyArray_Return(angles_arr);
}

int build_angles_arr(PyObject *distances_obj, PyArrayObject **angles_arr, int *size, int Nd, int force2Ddimension, int bidirectional, int *Na)
{
  PyArrayObject *distances_arr;
  int *distances, *angles;
  int Ndist;
  npy_intp dims[2];

  if (distances_obj)
  {
    // Interpret the distance object as numpy array
    distances_arr = (PyArrayObject *)PyArray_FROM_OTF(distances_obj, NPY_INT, NPY_ARRAY_FORCECAST | NPY_ARRAY_IN_ARRAY);

    if (!distances_arr)
    {
      PyErr_SetString(PyExc_RuntimeError, "Error parsing distances array.");
      return -1;
    }

    if (PyArray_NDIM(distances_arr) != 1)
    {
      Py_XDECREF(distances_arr);

      PyErr_SetString(PyExc_ValueError, "Expecting distances array to be 1-dimensional.");
      return -1;
    }

    // Get the number of distances and the distances array data
    Ndist = (int)PyArray_DIM(distances_arr, 0);
    distances = (int *)PyArray_DATA(distances_arr);
  }
  else
  {
    Ndist = 1;
    distances = (int *)malloc(sizeof *distances);
    if (!distances)
    {
      PyErr_NoMemory();
      return -1;
    }
    distances[0] = 1;
  }

  // Generate the angles needed for texture matrix calculation
  *Na = get_angle_count(size, distances, Nd, Ndist, bidirectional, force2Ddimension);  // 3D, mono-directional
  if ((*Na) == 0)
  {
    if (distances_obj)
      Py_XDECREF(distances_arr);
    else
      free(distances);

    PyErr_SetString(PyExc_RuntimeError, "Error getting angle count.");
    return -1;
  }

  // Wrap angles in a numpy array
  dims[0] = *Na;
  dims[1] = Nd;

  *angles_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_INT);
  if (!(*angles_arr))
  {
    if (distances_obj)
      Py_XDECREF(distances_arr);
    else
      free(distances);

    PyErr_NoMemory();
    return -1;
  }
  angles = (int *)PyArray_DATA(*angles_arr);

  if(build_angles(size, distances, Nd, Ndist, force2Ddimension, *Na, angles) > 0)
  {
    if (distances_obj)
      Py_XDECREF(distances_arr);
    else
      free(distances);

    Py_XDECREF(angles_arr);

    PyErr_SetString(PyExc_RuntimeError, "Error building angles.");
    return -1;
  }

  // Clean up distances array
  if (distances_obj)
    Py_XDECREF(distances_arr);
  else
    free(distances);

  return 0;
}

int try_parse_arrays(PyObject *image_obj, PyObject *mask_obj, PyArrayObject **image_arr, PyArrayObject **mask_arr, int *Nd, int **size, int **strides, int mask_flags)
{
  int d;

  // Interpret the image and mask objects as numpy arrays
  *image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_ARRAY_FORCECAST | NPY_ARRAY_IN_ARRAY);
  if (!(*image_arr))
    return -1;

  *mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BOOL, NPY_ARRAY_FORCECAST | NPY_ARRAY_IN_ARRAY | mask_flags);
  if (!(*mask_arr))
  {
    Py_XDECREF(*image_arr);
    return -1;
  }

  *Nd = (int)PyArray_NDIM(*image_arr);

  // Check if Image and Mask have 3 dimensions, and if Angles has 2 dimensions
  if (*Nd != PyArray_NDIM(*mask_arr))
  {
    Py_XDECREF(*image_arr);
    Py_XDECREF(*mask_arr);
    PyErr_SetString(PyExc_ValueError, "Expected a image and mask to have equal number of dimensions.");
    return -1;
  }

  *size = (int *)malloc(sizeof **size * *Nd);
  if (!(*size))
  {
    Py_XDECREF(*image_arr);
    Py_XDECREF(*mask_arr);
    PyErr_NoMemory();
    return -1;
  }

  *strides = (int *)malloc(sizeof **strides * *Nd);
  if (!(*strides))
  {
    free(*size);

    Py_XDECREF(*image_arr);
    Py_XDECREF(*mask_arr);
    PyErr_NoMemory();
    return -1;
  }

  for (d = 0; d < *Nd; d++)
  {
    (*size)[d] = (int)PyArray_DIM(*image_arr, d);
    if ((*size)[d] != (int)PyArray_DIM(*mask_arr, d))
    {
      free(*size);
      free(*strides);

      Py_XDECREF(*image_arr);
      Py_XDECREF(*mask_arr);
      PyErr_SetString(PyExc_ValueError, "Dimensions of image and mask do not match.");
      return -1;
    }
    (*strides)[d] = (int)(PyArray_STRIDE(*image_arr, d) / PyArray_ITEMSIZE(*image_arr));
  }
  return 0;
}
