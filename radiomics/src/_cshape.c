#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <stdlib.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "cshape.h"

static char module_docstring[] = "This module links to C-compiled code for efficient calculation of the surface area "
                                 "in the pyRadiomics package. It provides fast calculation using a marching cubes "
                                 "algortihm, accessed via ""calculate_surfacearea"". Arguments for this function"
                                 "are positional and consist of two numpy arrays, mask and pixelspacing, which must "
                                 "be supplied in this order. Pixelspacing is a 3 element vector containing the pixel"
                                 "spacing in z, y and x dimension, respectively. All non-zero elements in mask are "
                                 "considered to be a part of the segmentation and are included in the algorithm.";
static char surface_docstring[] = "Arguments: Mask, PixelSpacing, uses a marching cubes algorithm to calculate an "
                                  "approximation to the total surface area. The isovalue is considered to be situated "
                                  "midway between a voxel that is part of the segmentation and a voxel that is not.";
static char diameter_docstring[] = "Arguments: Mask, PixelSpacing, ROI size.";

static PyObject *cshape_calculate_surfacearea(PyObject *self, PyObject *args);
static PyObject *cshape_calculate_diameter(PyObject *self, PyObject *args);

int check_arrays(PyArrayObject *mask_arr, PyArrayObject *spacing_arr, int *size, int *strides);

static PyMethodDef module_methods[] = {
  //{"calculate_", cmatrices_, METH_VARARGS, _docstring},
  { "calculate_surfacearea", cshape_calculate_surfacearea, METH_VARARGS, surface_docstring },
  { "calculate_diameter", cshape_calculate_diameter,METH_VARARGS, diameter_docstring},
  { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_cshape",           /* m_name */
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
    m = Py_InitModule3("_cshape",
                       module_methods, module_docstring);
#endif

  if (m == NULL)
      return NULL;

  return m;
}

#if PY_MAJOR_VERSION < 3
  PyMODINIT_FUNC
  init_cshape(void)
  {
    // Initialize numpy functionality
    import_array();

    moduleinit();
  }
#else
  PyMODINIT_FUNC
  PyInit__cshape(void)
  {
    // Initialize numpy functionality
    import_array();

    return moduleinit();
  }
#endif

static PyObject *cshape_calculate_surfacearea(PyObject *self, PyObject *args)
{
  PyObject *mask_obj, *spacing_obj;
  PyArrayObject *mask_arr, *spacing_arr;
  int size[3];
  int strides[3];
  char *mask;
  double *spacing;
  double SA;
  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OO", &mask_obj, &spacing_obj))
    return NULL;

  // Interpret the input as numpy arrays
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BYTE, NPY_ARRAY_FORCECAST | NPY_ARRAY_UPDATEIFCOPY | NPY_ARRAY_IN_ARRAY);
  spacing_arr = (PyArrayObject *)PyArray_FROM_OTF(spacing_obj, NPY_DOUBLE, NPY_ARRAY_FORCECAST | NPY_ARRAY_UPDATEIFCOPY | NPY_ARRAY_IN_ARRAY);

  if (check_arrays(mask_arr, spacing_arr, size, strides) > 0) return NULL;

  // Get arrays in Ctype
  mask = (char *)PyArray_DATA(mask_arr);
  spacing = (double *)PyArray_DATA(spacing_arr);

  //Calculate Surface Area
  SA = calculate_surfacearea(mask, size, strides, spacing);

  // Clean up
  Py_XDECREF(mask_arr);
  Py_XDECREF(spacing_arr);

  if (SA < 0) // if SA < 0, an error has occurred
  {
    PyErr_SetString(PyExc_RuntimeError, "Calculation of Surface Area Failed.");
    return NULL;
  }

  return Py_BuildValue("f", SA);
}

static PyObject *cshape_calculate_diameter(PyObject *self, PyObject *args)
{
  PyObject *mask_obj, *spacing_obj;
  PyArrayObject *mask_arr, *spacing_arr;
  int Ns;
  int size[3];
  int strides[3];
  char *mask;
  double *spacing;
  double *diameters;
  PyObject *rslt;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOi", &mask_obj, &spacing_obj, &Ns))
    return NULL;

  // Interpret the input as numpy arrays
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BYTE, NPY_ARRAY_FORCECAST | NPY_ARRAY_UPDATEIFCOPY | NPY_ARRAY_IN_ARRAY);
  spacing_arr = (PyArrayObject *)PyArray_FROM_OTF(spacing_obj, NPY_DOUBLE, NPY_ARRAY_FORCECAST | NPY_ARRAY_UPDATEIFCOPY | NPY_ARRAY_IN_ARRAY);

  if (check_arrays(mask_arr, spacing_arr, size, strides) > 0) return NULL;

  // Get arrays in Ctype
  mask = (char *)PyArray_DATA(mask_arr);
  spacing = (double *)PyArray_DATA(spacing_arr);

  // Initialize output array (elements not set)
  diameters = (double *)calloc(4, sizeof(double));

  // Calculating Max 3D Diameter
  if (!calculate_diameter(mask, size, strides, spacing, Ns, diameters))
  {
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    free(diameters);
    PyErr_SetString(PyExc_RuntimeError, "Calculation of maximum 3D diameter failed.");
    return NULL;
  }

  rslt = Py_BuildValue("ffff", diameters[0], diameters[1], diameters[2], diameters[3]);

  // Clean up
  Py_XDECREF(mask_arr);
  Py_XDECREF(spacing_arr);
  free(diameters);

  return rslt;
}

int check_arrays(PyArrayObject *mask_arr, PyArrayObject *spacing_arr, int *size, int *strides)
{
  if (mask_arr == NULL || spacing_arr == NULL)
  {
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_SetString(PyExc_RuntimeError, "Error parsing array arguments.");
    return 1;
  }

  if (PyArray_NDIM(mask_arr) != 3 || PyArray_NDIM(spacing_arr) != 1)
  {
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_SetString(PyExc_RuntimeError, "Expected a 3D array for mask, 1D for spacing.");
    return 2;
  }

  if ( !PyArray_IS_C_CONTIGUOUS(mask_arr) || !PyArray_IS_C_CONTIGUOUS(spacing_arr))
  {
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_SetString(PyExc_RuntimeError, "Expecting input arrays to be C-contiguous.");
    return 3;
  }

  if (PyArray_DIM(spacing_arr, 0) != 3)
  {
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_SetString(PyExc_RuntimeError, "Expecting spacing array to have shape (3,).");
    return 4;
  }

  // Get sizes of the arrays
  size[2] = (int)PyArray_DIM(mask_arr, 2);
  size[1] = (int)PyArray_DIM(mask_arr, 1);
  size[0] = (int)PyArray_DIM(mask_arr, 0);

  strides[2] = (int)(PyArray_STRIDE(mask_arr, 2) / PyArray_ITEMSIZE(mask_arr));
  strides[1] = (int)(PyArray_STRIDE(mask_arr, 1) / PyArray_ITEMSIZE(mask_arr));
  strides[0] = (int)(PyArray_STRIDE(mask_arr, 0) / PyArray_ITEMSIZE(mask_arr));

  return 0;
}
