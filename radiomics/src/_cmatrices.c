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
	mask_arr = (PyArrayObject *)PyArray_FROM_OTF(mask_obj, NPY_UBYTE, NPY_FORCECAST | NPY_ENSURECOPY | NPY_IN_ARRAY);
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
	tempDims[0] = Ns;
	tempDims[1] = 2;
	temp_arr = (PyArrayObject *)PyArray_FromDims(2, (int*)tempDims, NPY_INT);

	// Get arrays in Ctype
	image = (int *)PyArray_DATA(image_arr);
	mask = (signed char *)PyArray_DATA(mask_arr);
	angles = (int *)PyArray_DATA(angles_arr);
	tempData = (int *)PyArray_DATA(temp_arr);

	// Set all elements to 0
	tempData_size = Ns * 2;
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
