#include <Python.h>
#include <numpy/arrayobject.h>
#include "cmatrices.h"

static char module_docstring[] = ("This module links to C-compiled code for efficient calculation of various matrices "
                                 "in the pyRadiomics package. It provides fast calculation for GLCM, GLDM, NGTDM, "
                                 "GLRLM and GLSZM. All functions are given names as follows: ""calculate_<Matrix>"", "
                                 "where <Matrix> is the name of the matix, in lowercase. Arguments for these functions "
                                 "are positional and start with 3 numpy arrays (image, mask and angles) and 1 integer "
                                 "(Ng, number of gray levels). Optionally extra arguments may be required, see function "
                                 "docstrings for detailed information.");
static char glcm_docstring[] = "Arguments: Image, Mask, Angles, Ng.";
static char gldm_docstring[] = "Arguments: Image, Mask, Angles, Ng, Alpha.";
static char ngtdm_docstring[] = "Arguments: Image, Mask, Angles, Ng.";
static char glszm_docstring[] = "Arguments: Image, Mask, Angles, Ng, Ns, matrix is cropped to maximum size encountered.";
static char gldzm_docstring[] = "Arguments: Image, Mask, DistanceMap, Angles, Ng, Nd.";
static char glrlm_docstring[] = "Arguments: Image, Mask, Angles, Ng, Nr.";

static PyObject *cmatrices_calculate_glcm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_gldm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_ngtdm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_glszm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_gldzm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_glrlm(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
  //{"calculate_", cmatrices_, METH_VARARGS, _docstring},
  { "calculate_glcm", cmatrices_calculate_glcm, METH_VARARGS, glcm_docstring },
  { "calculate_gldm", cmatrices_calculate_gldm, METH_VARARGS, gldm_docstring },
  { "calculate_ngtdm", cmatrices_calculate_ngtdm, METH_VARARGS, ngtdm_docstring },
  { "calculate_glszm", cmatrices_calculate_glszm, METH_VARARGS, glszm_docstring },
  { "calculate_gldzm", cmatrices_calculate_gldzm, METH_VARARGS, gldzm_docstring },
  { "calculate_glrlm", cmatrices_calculate_glrlm, METH_VARARGS, glrlm_docstring },
  { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC
init_cmatrices(void) {
  PyObject *m = Py_InitModule3("_cmatrices", module_methods, module_docstring);
  if (m == NULL)
    return;

  // Initialize numpy functionality
  import_array();
}

static PyObject *cmatrices_calculate_glcm(PyObject *self, PyObject *args)
{
  int Ng;
  PyObject *image_obj, *mask_obj, *angles_obj;
  PyArrayObject *image_arr, *mask_arr, *angles_arr;
  int Sx, Sy, Sz, Na;
  int dims[3];
  PyArrayObject *glcm_arr;
  int *image;
  char *mask;
  int *angles;
  double *glcm;
  int glcm_size, k;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOOi", &image_obj, &mask_obj, &angles_obj, &Ng))
    return NULL;

  // Interpret the input as numpy arrays
  image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BOOL, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

  if (image_arr == NULL || mask_arr == NULL || angles_arr == NULL)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);
    return NULL;
  }

  // Check if Image and Mask have 3 dimensions
  if (PyArray_NDIM(image_arr) != 3 || PyArray_NDIM(mask_arr) != 3)
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Expected a 3D array for image and mask.");
    return NULL;
  }

  // Get sizes of the arrays
  Sx = (int)PyArray_DIM(image_arr, 2);
  Sy = (int)PyArray_DIM(image_arr, 1);
  Sz = (int)PyArray_DIM(image_arr, 0);

  Na = (int)PyArray_DIM(angles_arr, 0);

  // Check if image and mask are the same size
  if (Sx != (int)PyArray_DIM(mask_arr, 2) || Sy != (int)PyArray_DIM(mask_arr, 1) || Sz != (int)PyArray_DIM(mask_arr, 0))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Dimensions of image and mask do not match.");
    return NULL;
  }

  // Initialize output array (elements not set)
  dims[0] = Ng;
  dims[1] = Ng;
  dims[2] = Na;
  glcm_arr = (PyArrayObject *)PyArray_FromDims(3, (int*)dims, NPY_DOUBLE);

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (char *)PyArray_DATA(mask_arr);
  angles = (int *)PyArray_DATA(angles_arr);
  glcm = (double *)PyArray_DATA(glcm_arr);

  // Set all elements to 0
  glcm_size = Ng * Ng * Na;
  for (k = 0; k < glcm_size; k++) glcm[k] = 0;

  //Calculate GLCM
  if (!calculate_glcm(image, mask, Sx, Sy, Sz, angles, Na, glcm, Ng))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    Py_DECREF(glcm_arr);
    PyErr_SetString(PyExc_RuntimeError, "Calculation of GLCM Failed.");
    return NULL;
  }

  // Clean up
  Py_DECREF(image_arr);
  Py_DECREF(mask_arr);
  Py_DECREF(angles_arr);

  return PyArray_Return(glcm_arr);
}

static PyObject *cmatrices_calculate_gldm(PyObject *self, PyObject *args)
{
  int Ng, alpha;
  PyObject *image_obj, *mask_obj, *angles_obj;
  PyArrayObject *image_arr, *mask_arr, *angles_arr;
  int Sx, Sy, Sz, Na;
  int dims[2];
  PyArrayObject *gldm_arr;
  int *image;
  char *mask;
  int *angles;
  double *gldm;
  int gldm_size, k;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOOii", &image_obj, &mask_obj, &angles_obj, &Ng, &alpha))
    return NULL;

  // Interpret the input as numpy arrays
  image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BOOL, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

  if (image_arr == NULL || mask_arr == NULL || angles_arr == NULL)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);
    return NULL;
  }

  // Check if Image and Mask have 3 dimensions
  if (PyArray_NDIM(image_arr) != 3 || PyArray_NDIM(mask_arr) != 3)
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Expected a 3D array for image and mask.");
    return NULL;
  }

  // Get sizes of the arrays
  Sx = (int)PyArray_DIM(image_arr, 2);
  Sy = (int)PyArray_DIM(image_arr, 1);
  Sz = (int)PyArray_DIM(image_arr, 0);

  Na = (int)PyArray_DIM(angles_arr, 0);

  // Check if image and mask are the same size
  if (Sx != (int)PyArray_DIM(mask_arr, 2) || Sy != (int)PyArray_DIM(mask_arr, 1) || Sz != (int)PyArray_DIM(mask_arr, 0))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Dimensions of image and mask do not match.");
    return NULL;
  }

  // Initialize output array (elements not set)
  dims[0] = Ng;
  dims[1] = Na * 2 + 1;
  gldm_arr = (PyArrayObject *)PyArray_FromDims(2, (int*)dims, NPY_DOUBLE);

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (char *)PyArray_DATA(mask_arr);
  angles = (int *)PyArray_DATA(angles_arr);
  gldm = (double *)PyArray_DATA(gldm_arr);

  // Set all elements to 0
  gldm_size = Ng * (Na * 2 + 1);
  for (k = 0; k < gldm_size; k++) gldm[k] = 0;

  //Calculate GLDM
  if (!calculate_gldm(image, mask, Sx, Sy, Sz, angles, Na, gldm, Ng, alpha))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    Py_DECREF(gldm_arr);
    PyErr_SetString(PyExc_RuntimeError, "Calculation of GLDM Failed.");
    return NULL;
  }

  // Clean up
  Py_DECREF(image_arr);
  Py_DECREF(mask_arr);
  Py_DECREF(angles_arr);

  return PyArray_Return(gldm_arr);
}

static PyObject *cmatrices_calculate_ngtdm(PyObject *self, PyObject *args)
{
  int Ng;
  PyObject *image_obj, *mask_obj, *angles_obj;
  PyArrayObject *image_arr, *mask_arr, *angles_arr;
  int Sx, Sy, Sz, Na;
  int dims[2];
  PyArrayObject *ngtdm_arr;
  int *image;
  char *mask;
  int *angles;
  double *ngtdm;
  int ngtdm_size, k;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOOi", &image_obj, &mask_obj, &angles_obj, &Ng))
    return NULL;

  // Interpret the input as numpy arrays
  image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BOOL, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

  if (image_arr == NULL || mask_arr == NULL || angles_arr == NULL)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);
    return NULL;
  }

  // Check if Image and Mask have 3 dimensions
  if (PyArray_NDIM(image_arr) != 3 || PyArray_NDIM(mask_arr) != 3)
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Expected a 3D array for image and mask.");
    return NULL;
  }

  // Get sizes of the arrays
  Sx = (int)PyArray_DIM(image_arr, 2);
  Sy = (int)PyArray_DIM(image_arr, 1);
  Sz = (int)PyArray_DIM(image_arr, 0);

  Na = (int)PyArray_DIM(angles_arr, 0);

  // Check if image and mask are the same size
  if (Sx != (int)PyArray_DIM(mask_arr, 2) || Sy != (int)PyArray_DIM(mask_arr, 1) || Sz != (int)PyArray_DIM(mask_arr, 0))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Dimensions of image and mask do not match.");
    return NULL;
  }

  // Initialize output array (elements not set)
  dims[0] = Ng;
  dims[1] = 3;
  ngtdm_arr = (PyArrayObject *)PyArray_FromDims(2, (int*)dims, NPY_DOUBLE);

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (char *)PyArray_DATA(mask_arr);
  angles = (int *)PyArray_DATA(angles_arr);
  ngtdm = (double *)PyArray_DATA(ngtdm_arr);

  // Set all elements to 0
  ngtdm_size = Ng * 3;
  for (k = 0; k < ngtdm_size; k++) ngtdm[k] = 0;

  //Calculate NGTDM
  if (!calculate_ngtdm(image, mask, Sx, Sy, Sz, angles, Na, ngtdm, Ng))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    Py_DECREF(ngtdm_arr);
    PyErr_SetString(PyExc_RuntimeError, "Calculation of NGTDM Failed.");
    return NULL;
  }

  // Clean up
  Py_DECREF(image_arr);
  Py_DECREF(mask_arr);
  Py_DECREF(angles_arr);

  return PyArray_Return(ngtdm_arr);
}

static PyObject *cmatrices_calculate_glszm(PyObject *self, PyObject *args)
{
  int Ng, Ns;
  PyObject *image_obj, *mask_obj, *angles_obj;
  PyArrayObject *image_arr, *mask_arr, *angles_arr;
  int Sx, Sy, Sz, Na;
  int tempDims[2];
  PyArrayObject *temp_arr;
  int *image;
  signed char *mask;
  int *angles;
  int * tempData;
  int tempData_size, k;
  int maxRegion;
  int dims[2];
  PyArrayObject *glszm_arr;
  double *glszm;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOOii", &image_obj, &mask_obj, &angles_obj, &Ng, &Ns))
    return NULL;

  // Interpret the input as numpy arrays
  image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY |NPY_IN_ARRAY);
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BYTE, NPY_FORCECAST | NPY_ENSURECOPY | NPY_IN_ARRAY);
  angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

  if (image_arr == NULL || mask_arr == NULL || angles_arr == NULL)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);
    return NULL;
  }

  // Check if Image and Mask have 3 dimensions
  if (PyArray_NDIM(image_arr) != 3 || PyArray_NDIM(mask_arr) != 3)
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Expected a 3D array for image and mask.");
    return NULL;
  }

  // Get sizes of the arrays
  Sx = (int)PyArray_DIM(image_arr, 2);
  Sy = (int)PyArray_DIM(image_arr, 1);
  Sz = (int)PyArray_DIM(image_arr, 0);

  Na = (int)PyArray_DIM(angles_arr, 0);

  // Check if image and mask are the same size
  if (Sx != (int)PyArray_DIM(mask_arr, 2) || Sy != (int)PyArray_DIM(mask_arr, 1) || Sz != (int)PyArray_DIM(mask_arr, 0))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Dimensions of image and mask do not match.");
    return NULL;
  }

  // Initialize temporary output array (elements not set)
  // add +1 to the first dimension so in the case every voxel represents a separate region,
  // tempData still contains a -1 element at the end
  tempDims[0] = Ns + 1;
  tempDims[1] = 2;
  temp_arr = (PyArrayObject *)PyArray_FromDims(2, (int*)tempDims, NPY_INT);

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (signed char *)PyArray_DATA(mask_arr);
  angles = (int *)PyArray_DATA(angles_arr);
  tempData = (int *)PyArray_DATA(temp_arr);

  // Set all elements to 0
  tempData_size = tempDims[0] * tempDims[1];
  for (k = 0; k < tempData_size; k++) tempData[k] = -1;

  //Calculate GLSZM
  maxRegion = 0;
  maxRegion = calculate_glszm(image, mask, Sx, Sy, Sz, angles, Na, tempData, Ng, Ns);
  if (maxRegion == -1) // Error occured
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    Py_DECREF(temp_arr);
    PyErr_SetString(PyExc_RuntimeError, "Calculation of GLSZM Failed.");
    return NULL;
  }

  // Clean up image, mask and angles arrays (not needed anymore)
  Py_DECREF(image_arr);
  Py_DECREF(mask_arr);
  Py_DECREF(angles_arr);

  // Initialize temporary output array (elements not set)
  if (maxRegion == 0) maxRegion = 1;
  dims[0] = Ng;
  dims[1] = maxRegion;
  glszm_arr = (PyArrayObject *)PyArray_FromDims(2, (int*)dims, NPY_DOUBLE);
  glszm = (double *)PyArray_DATA(glszm_arr);

  if (!fill_glszm(tempData, glszm, Ng, maxRegion))
  {
    Py_DECREF(temp_arr);
    Py_DECREF(glszm_arr);
    PyErr_SetString(PyExc_RuntimeError, "Error filling GLSZM.");
  return NULL;
  }

  // Clean up tempData
  Py_DECREF(temp_arr);

  return PyArray_Return(glszm_arr);
}

static PyObject *cmatrices_calculate_gldzm(PyObject *self, PyObject *args)
{
  int Ng, Nd;
  PyObject *image_obj, *mask_obj, *distmap_obj, *angles_obj;
  PyArrayObject *image_arr, *mask_arr, *distmap_arr, *angles_arr;
  int size[3];
  int Na;
  int dims[2];
  PyArrayObject *gldzm_arr;

  int *image;
  signed char *mask;
  int *distmap;
  int *angles;

  int gldzm_size, k;
  double *gldzm;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOOOii", &image_obj, &mask_obj, &distmap_obj, &angles_obj, &Ng, &Nd))
    return NULL;

  // Interpret the input as numpy arrays
  image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BYTE, NPY_FORCECAST | NPY_ENSURECOPY | NPY_IN_ARRAY);
  distmap_arr = (PyArrayObject *)PyArray_FROM_OTF(distmap_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

  if (image_arr == NULL || mask_arr == NULL || distmap_arr == NULL || angles_arr == NULL)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(distmap_arr);
    Py_XDECREF(angles_arr);
    return NULL;
  }

  // Check if Image and Mask have 3 dimensions
  if (PyArray_NDIM(image_arr) != 3 || PyArray_NDIM(mask_arr) != 3 || PyArray_NDIM(distmap_arr) != 3)
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(distmap_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Expected a 3D array for image, mask and distance map.");
    return NULL;
  }

  // Get sizes of the arrays
  size[2] = (int)PyArray_DIM(image_arr, 2);
  size[1] = (int)PyArray_DIM(image_arr, 1);
  size[0] = (int)PyArray_DIM(image_arr, 0);

  Na = (int)PyArray_DIM(angles_arr, 0);

  // Check if image and mask are the same size
  if (size[2] != (int)PyArray_DIM(mask_arr, 2) || size[1] != (int)PyArray_DIM(mask_arr, 1) || size[0] != (int)PyArray_DIM(mask_arr, 0) &&
    size[2] != (int)PyArray_DIM(distmap_arr, 2) || size[1] != (int)PyArray_DIM(distmap_arr, 1) || size[0] != (int)PyArray_DIM(distmap_arr, 0))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(distmap_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Dimensions of image, mask and distance map do not match.");
    return NULL;
  }

  // Initialize temporary output array (elements not set)
  dims[0] = Ng;
  dims[1] = Nd;
  gldzm_arr = (PyArrayObject *)PyArray_FromDims(2, (int*)dims, NPY_DOUBLE);

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (signed char *)PyArray_DATA(mask_arr);
  distmap = (int *)PyArray_DATA(distmap_arr);
  angles = (int *)PyArray_DATA(angles_arr);
  gldzm = (double *)PyArray_DATA(gldzm_arr);

  // Set all elements to 0
  gldzm_size = Nd * Ng;
  for (k = 0; k < gldzm_size; k++) gldzm[k] = 0;

  //Calculate GLDZM
  if (!calculate_gldzm(image, mask, distmap, size, angles, Na, gldzm, Ng, Nd)) // Error occured
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(distmap_arr);
    Py_DECREF(angles_arr);
    Py_DECREF(gldzm_arr);
    PyErr_SetString(PyExc_RuntimeError, "Calculation of GLDZM Failed.");
    return NULL;
  }

  // Clean up
  Py_DECREF(image_arr);
  Py_DECREF(mask_arr);
  Py_DECREF(distmap_arr);
  Py_DECREF(angles_arr);

  return PyArray_Return(gldzm_arr);
}

static PyObject *cmatrices_calculate_glrlm(PyObject *self, PyObject *args)
{
  int Ng, Nr;
  PyObject *image_obj, *mask_obj, *angles_obj;
  PyArrayObject *image_arr, *mask_arr, *angles_arr;
  int Na;
  int size[3];
  int strides[3];
  int dims[3];
  PyArrayObject *glrlm_arr;
  int *image;
  char *mask;
  int *angles;
  double *glrlm;
  int glrlm_size, k;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOOii", &image_obj, &mask_obj, &angles_obj, &Ng, &Nr))
    return NULL;

  // Interpret the input as numpy arrays
  image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BOOL, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
  angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

  if (image_arr == NULL || mask_arr == NULL || angles_arr == NULL)
  {
    Py_XDECREF(image_arr);
    Py_XDECREF(mask_arr);
    Py_XDECREF(angles_arr);
    return NULL;
  }

  // Check if Image and Mask have 3 dimensions
  if (PyArray_NDIM(image_arr) != 3 || PyArray_NDIM(mask_arr) != 3)
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Expected a 3D array for image and mask.");
    return NULL;
  }

  // Get sizes of the arrays
  size[2] = (int)PyArray_DIM(image_arr, 2);
  size[1] = (int)PyArray_DIM(image_arr, 1);
  size[0] = (int)PyArray_DIM(image_arr, 0);

  strides[2] = (int)(image_arr->strides[2] / image_arr->descr->elsize);
  strides[1] = (int)(image_arr->strides[1] / image_arr->descr->elsize);
  strides[0] = (int)(image_arr->strides[0] / image_arr->descr->elsize);

  Na = (int)PyArray_DIM(angles_arr, 0);

  // Check if image and mask are the same size
  if (size[2] != (int)PyArray_DIM(mask_arr, 2) || size[1] != (int)PyArray_DIM(mask_arr, 1) || size[0] != (int)PyArray_DIM(mask_arr, 0))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    PyErr_SetString(PyExc_RuntimeError, "Dimensions of image and mask do not match.");
    return NULL;
  }

  // Initialize output array (elements not set)
  dims[0] = Ng;
  dims[1] = Nr;
  dims[2] = Na;
  glrlm_arr = (PyArrayObject *)PyArray_FromDims(3, (int*)dims, NPY_DOUBLE);

  // Get arrays in Ctype
  image = (int *)PyArray_DATA(image_arr);
  mask = (char *)PyArray_DATA(mask_arr);
  angles = (int *)PyArray_DATA(angles_arr);
  glrlm = (double *)PyArray_DATA(glrlm_arr);

  // Set all elements to 0
  glrlm_size = Ng * Nr * Na;
  for (k = 0; k < glrlm_size; k++) glrlm[k] = 0;

  //Calculate GLRLM
  if (!calculate_glrlm(image, mask, size, strides, angles, Na, glrlm, Ng, Nr))
  {
    Py_DECREF(image_arr);
    Py_DECREF(mask_arr);
    Py_DECREF(angles_arr);
    Py_DECREF(glrlm_arr);
    PyErr_SetString(PyExc_RuntimeError, "Calculation of GLRLM Failed.");
    return NULL;
  }

  // Clean up
  Py_DECREF(image_arr);
  Py_DECREF(mask_arr);
  Py_DECREF(angles_arr);

  return PyArray_Return(glrlm_arr);
}
