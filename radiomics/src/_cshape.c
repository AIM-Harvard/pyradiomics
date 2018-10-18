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
static char coefficients_docstring[] = "Arguments: Mask, PixelSpacing. Uses a marching cubes algorithm to calculate an "
                                       "approximation to the total surface area, volume and maximum diameters. "
                                       "The isovalue is considered to be situated midway between a voxel that is part "
                                       "of the segmentation and a voxel that is not.";
static char coefficients2D_docstring[] = "Arguments: Mask, PixelSpacing. Uses an adapted 2D marching cubes algorithm "
                                         "to calculate an approximation to the total perimeter, surface and maximum "
                                         "diameter. The isovalue is considered to be situated midway between a pixel "
                                         "that is part of the segmentation and a pixel that is not.";

static PyObject *cshape_calculate_coefficients(PyObject *self, PyObject *args);
static PyObject *cshape_calculate_coefficients2D(PyObject *self, PyObject *args);

int check_arrays(PyArrayObject *mask_arr, PyArrayObject *spacing_arr, int *size, int *strides, int dimension);

static PyMethodDef module_methods[] = {
  //{"calculate_", cmatrices_, METH_VARARGS, _docstring},
  { "calculate_coefficients", cshape_calculate_coefficients, METH_VARARGS, coefficients_docstring },
  { "calculate_coefficients2D", cshape_calculate_coefficients2D, METH_VARARGS, coefficients2D_docstring },
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

static PyObject *cshape_calculate_coefficients(PyObject *self, PyObject *args)
{
  PyObject *mask_obj, *spacing_obj;
  PyArrayObject *mask_arr, *spacing_arr;
  int size[3];
  int strides[3];
  char *mask;
  double *spacing;
  double SA, Volume;
  double diameters[4];
  PyObject *diameter_obj;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OO", &mask_obj, &spacing_obj))
    return NULL;

  // Interpret the input as numpy arrays
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BYTE, NPY_ARRAY_FORCECAST | NPY_ARRAY_IN_ARRAY);
  spacing_arr = (PyArrayObject *)PyArray_FROM_OTF(spacing_obj, NPY_DOUBLE, NPY_ARRAY_FORCECAST | NPY_ARRAY_IN_ARRAY);

  if (check_arrays(mask_arr, spacing_arr, size, strides, 3) > 0) return NULL;

  // Get arrays in Ctype
  mask = (char *)PyArray_DATA(mask_arr);
  spacing = (double *)PyArray_DATA(spacing_arr);

  //Calculate Surface Area and volume
  if (calculate_coefficients(mask, size, strides, spacing, &SA, &Volume, diameters))
  {
    // An error has occurred
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_SetString(PyExc_RuntimeError, "Calculation of Shape coefficients failed.");
    return NULL;
  }

  // Clean up
  Py_XDECREF(mask_arr);
  Py_XDECREF(spacing_arr);

  diameter_obj = Py_BuildValue("ffff", diameters[0], diameters[1], diameters[2], diameters[3]);
  return Py_BuildValue("ffN", SA, Volume, diameter_obj);
}

static PyObject *cshape_calculate_coefficients2D(PyObject *self, PyObject *args)
{
  PyObject *mask_obj, *spacing_obj;
  PyArrayObject *mask_arr, *spacing_arr;
  int size[2];
  int strides[2];
  char *mask;
  double *spacing;
  double perimeter, surface, diameter;

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OO", &mask_obj, &spacing_obj))
    return NULL;

  // Interpret the input as numpy arrays
  mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_BYTE, NPY_ARRAY_FORCECAST | NPY_ARRAY_IN_ARRAY);
  spacing_arr = (PyArrayObject *)PyArray_FROM_OTF(spacing_obj, NPY_DOUBLE, NPY_ARRAY_FORCECAST | NPY_ARRAY_IN_ARRAY);

  if (check_arrays(mask_arr, spacing_arr, size, strides, 2) > 0) return NULL;

  // Get arrays in Ctype
  mask = (char *)PyArray_DATA(mask_arr);
  spacing = (double *)PyArray_DATA(spacing_arr);

  //Calculate Surface Area and volume
  if (calculate_coefficients2D(mask, size, strides, spacing, &perimeter, &surface, &diameter))
  {
    // An error has occurred
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_SetString(PyExc_RuntimeError, "Calculation of Shape coefficients failed.");
    return NULL;
  }

  // Clean up
  Py_XDECREF(mask_arr);
  Py_XDECREF(spacing_arr);

  return Py_BuildValue("fff", perimeter, surface, diameter);
}

int check_arrays(PyArrayObject *mask_arr, PyArrayObject *spacing_arr, int *size, int *strides, int dimension)
{
  int i;

  if (mask_arr == NULL || spacing_arr == NULL)
  {
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_SetString(PyExc_RuntimeError, "Error parsing array arguments.");
    return 1;
  }

  if (PyArray_NDIM(mask_arr) != dimension || PyArray_NDIM(spacing_arr) != 1)
  {
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_Format(PyExc_ValueError, "Expected a %iD array for mask, 1D for spacing.", dimension);
    return 2;
  }

  if ( !PyArray_IS_C_CONTIGUOUS(mask_arr) || !PyArray_IS_C_CONTIGUOUS(spacing_arr))
  {
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_SetString(PyExc_ValueError, "Expecting input arrays to be C-contiguous.");
    return 3;
  }

  if (PyArray_DIM(spacing_arr, 0) != PyArray_NDIM(mask_arr))
  {
    Py_XDECREF(mask_arr);
    Py_XDECREF(spacing_arr);
    PyErr_SetString(PyExc_ValueError, "Expecting spacing array to have shape (3,).");
    return 4;
  }

  // Get sizes and strides of the arrays
  for (i = 0; i < dimension; i++)
  {
    size[i] = (int)PyArray_DIM(mask_arr, i);
    strides[i] = (int)(PyArray_STRIDE(mask_arr, i) / PyArray_ITEMSIZE(mask_arr));
  }

  return 0;
}
