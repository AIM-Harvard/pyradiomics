#include <Python.h>
#include <numpy/arrayobject.h>
#include "cmatrices.h"

static char module_docstring[] = "";
static char glcm_docstring[] = "";
static char gldm_docstring[] = "";
static char ngtdm_docstring[] = "";
static char glszm_docstring[] = "";

static PyObject *cmatrices_calculate_glcm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_gldm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_ngtdm(PyObject *self, PyObject *args);
static PyObject *cmatrices_calculate_glszm(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
	//{"calculate_", cmatrices_, METH_VARARGS, _docstring},
	{ "calculate_glcm", cmatrices_calculate_glcm, METH_VARARGS, glcm_docstring },
	{ "calculate_gldm", cmatrices_calculate_gldm, METH_VARARGS, gldm_docstring },
	{ "calculate_ngtdm", cmatrices_calculate_ngtdm, METH_VARARGS, ngtdm_docstring },
	{ "calculate_glszm", cmatrices_calculate_glszm, METH_VARARGS, glszm_docstring },
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
	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "OOOi", &image_obj, &mask_obj, &angles_obj, &Ng))
		return NULL;

	// Interpret the input as numpy arrays
	PyArrayObject *image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
	PyArrayObject *mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
	PyArrayObject *angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

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
	int Sx = (int)PyArray_DIM(image_arr, 2);
	int Sy = (int)PyArray_DIM(image_arr, 1);
	int Sz = (int)PyArray_DIM(image_arr, 0);

	int Na = (int)PyArray_DIM(angles_arr, 0);

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
	int dims[3] = { Ng, Ng, Na };
	PyArrayObject *glcm_arr = (PyArrayObject *)PyArray_FromDims(3, (int*)dims, NPY_DOUBLE);

	// Get arrays in Ctype
	int *image = (int *)PyArray_DATA(image_arr);
	int *mask = (int *)PyArray_DATA(mask_arr);
	int *angles = (int *)PyArray_DATA(angles_arr);
	double *glcm = (double *)PyArray_DATA(glcm_arr);

	// Set all elements to 0
	int glcm_size = Ng * Ng * Na;
	for (int k = 0; k < glcm_size; k++) glcm[k] = 0;

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
	int Ng, Nd, alpha;
	PyObject *image_obj, *mask_obj, *angles_obj;
	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "OOOiii", &image_obj, &mask_obj, &angles_obj, &Ng, &Nd, &alpha))
		return NULL;

	// Interpret the input as numpy arrays
	PyArrayObject *image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
	PyArrayObject *mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
	PyArrayObject *angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

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
	int Sx = (int)PyArray_DIM(image_arr, 2);
	int Sy = (int)PyArray_DIM(image_arr, 1);
	int Sz = (int)PyArray_DIM(image_arr, 0);

	int Na = (int)PyArray_DIM(angles_arr, 0);

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
	int dims[2] = { Ng, Nd };
	PyArrayObject *gldm_arr = (PyArrayObject *)PyArray_FromDims(2, (int*)dims, NPY_DOUBLE);

	// Get arrays in Ctype
	int *image = (int *)PyArray_DATA(image_arr);
	int *mask = (int *)PyArray_DATA(mask_arr);
	int *angles = (int *)PyArray_DATA(angles_arr);
	double *gldm = (double *)PyArray_DATA(gldm_arr);

	// Set all elements to 0
	int gldm_size = Ng * Nd;
	for (int k = 0; k < gldm_size; k++) gldm[k] = 0;

	//Calculate GLDM
	if (!calculate_gldm(image, mask, Sx, Sy, Sz, angles, Na, gldm, Ng, Nd, alpha))
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
	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "OOOi", &image_obj, &mask_obj, &angles_obj, &Ng))
		return NULL;

	// Interpret the input as numpy arrays
	PyArrayObject *image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
	PyArrayObject *mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
	PyArrayObject *angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

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
	int Sx = (int)PyArray_DIM(image_arr, 2);
	int Sy = (int)PyArray_DIM(image_arr, 1);
	int Sz = (int)PyArray_DIM(image_arr, 0);

	int Na = (int)PyArray_DIM(angles_arr, 0);

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
	int dims[2] = { Ng, 3 };
	PyArrayObject *ngtdm_arr = (PyArrayObject *)PyArray_FromDims(2, (int*)dims, NPY_DOUBLE);

	// Get arrays in Ctype
	int *image = (int *)PyArray_DATA(image_arr);
	int *mask = (int *)PyArray_DATA(mask_arr);
	int *angles = (int *)PyArray_DATA(angles_arr);
	double *ngtdm = (double *)PyArray_DATA(ngtdm_arr);

	// Set all elements to 0
	int ngtdm_size = Ng * 3;
	for (int k = 0; k < ngtdm_size; k++) ngtdm[k] = 0;

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
	int Ng;
	PyObject *image_obj, *mask_obj, *angles_obj;
	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "OOOi", &image_obj, &mask_obj, &angles_obj, &Ng))
		return NULL;

	// Interpret the input as numpy arrays
	PyArrayObject *image_arr = (PyArrayObject *)PyArray_FROM_OTF(image_obj, NPY_INT, NPY_FORCECAST | NPY_ENSURECOPY |NPY_IN_ARRAY);
	PyArrayObject *mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);
	PyArrayObject *angles_arr = (PyArrayObject *)PyArray_FROM_OTF(angles_obj, NPY_INT, NPY_FORCECAST | NPY_UPDATEIFCOPY | NPY_IN_ARRAY);

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
	int Sx = (int)PyArray_DIM(image_arr, 2);
	int Sy = (int)PyArray_DIM(image_arr, 1);
	int Sz = (int)PyArray_DIM(image_arr, 0);
	int Ns = Sx * Sy * Sz;

	int Na = (int)PyArray_DIM(angles_arr, 0);

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
	int dims[2] = { Ng, Ns };
	PyArrayObject *glszm_arr = (PyArrayObject *)PyArray_FromDims(2, (int*)dims, NPY_DOUBLE);

	// Get arrays in Ctype
	int *image = (int *)PyArray_DATA(image_arr);
	int *mask = (int *)PyArray_DATA(mask_arr);
	int *angles = (int *)PyArray_DATA(angles_arr);
	double *glszm = (double *)PyArray_DATA(glszm_arr);

	// Set all elements to 0
	int ngtdm_size = Ng * Ns;
	for (int k = 0; k < ngtdm_size; k++) glszm[k] = 0;

	//Calculate NGTDM
	if (!calculate_glszm(image, mask, Sx, Sy, Sz, angles, Na, glszm, Ng, Ns))
	{
		Py_DECREF(image_arr);
		Py_DECREF(mask_arr);
		Py_DECREF(angles_arr);
		Py_DECREF(glszm_arr);
		PyErr_SetString(PyExc_RuntimeError, "Calculation of GLSZM Failed.");
		return NULL;
	}

	// Clean up
	Py_DECREF(image_arr);
	Py_DECREF(mask_arr);
	Py_DECREF(angles_arr);

	return PyArray_Return(glszm_arr);
}
