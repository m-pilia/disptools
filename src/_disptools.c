/*!
 * Python C extension.
 *
 * This file defines a Python C extension module that wraps the C library.
 */

// NumPy API version
#define NPY_NO_DEPRECATED_API NPY_1_14_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

#include "headers/field.h"
#include "headers/jacobian.h"
#include "headers/shape_descriptors.h"
#include "headers/displacement_field_gradient.h"
#include "headers/displacement_field_greedy.h"
#include "headers/VolumeMatching3D.h"

// String literals
#define CUBENESS "cubeness"
#define OCTAHEDRONESS "octahedroness"
#define ALGORITHM_GRADIENT "gradient"
#define ALGORITHM_GREEDY "greedy"
#define ALGORITHM_MATCHING "matching"

// Float size
#if FLOAT_SIZE == 32
    #define PY_FMT_FLOATING "f"
    #define PY_FMT_REGULARISE "Of"
    #define PY_FMT_SHAPE_DESCRIPTOR "OOOs"
    #define PY_FMT_JACOBIAN "OOO"
    #define PY_FMT_DISPLACEMENT "OOOffffffffplOs"
    #define NUMPY_FLOATING_TYPE NPY_FLOAT
#else
    #define PY_FMT_FLOATING "d"
    #define PY_FMT_REGULARISE "Od"
    #define PY_FMT_SHAPE_DESCRIPTOR "OOOs"
    #define PY_FMT_JACOBIAN "OOO"
    #define PY_FMT_DISPLACEMENT "OOOddddddddplOs"
    #define NUMPY_FLOATING_TYPE NPY_DOUBLE
#endif


/*!
 * Read a tuple of 3 floats and return an array of FLOATING.
 */
bool read_triplet(PyObject *tuple, FLOATING *result, const int expected_length, const char *parameter_name) {

    int len = PyTuple_Size(tuple);
    if (expected_length != len) {
        PyErr_Format(PyExc_ValueError, "Wrong number of components (%d) for %s.", len, parameter_name);
        return true;
    }

    for (int i = 0; i < len; ++i) {
        result[i] = (FLOATING) PyFloat_AsDouble(PyTuple_GetItem(tuple, i));
        if (PyErr_Occurred()) {
            PyErr_Format(PyExc_ValueError, "%s must be a tuple of floats.", parameter_name);
            return true;
        }
    }

    return false;
}

/*!
 * Method that returns the size (in bit) of the underlying C floating
 * point type.
 */
static PyObject *method_get_float_type_size(PyObject *module) {
    (void) module;
    return Py_BuildValue("i", 8 * sizeof (FLOATING));
}

/*!
 * Method to apply a lower threshold to the Jacobian.
 */
static PyObject *method_regularise(PyObject *self, PyObject *args)
{
    const FLOATING epsilon;
    PyArrayObject *jacobian = NULL;

    (void) self;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, PY_FMT_REGULARISE, &jacobian, &epsilon)) {
        return NULL;
    }

    /* Error check */
    if (!jacobian) {
        return NULL;
    }

    /* Get the shape */
    const int nz = PyArray_DIM(jacobian, 0);
    const int ny = PyArray_DIM(jacobian, 1);
    const int nx = PyArray_DIM(jacobian, 2);

    /* Get the data pointers */
    FLOATING *jacobian_data = (FLOATING*) PyArray_DATA(jacobian);

    /* Call the library function to compute the Jacobian */
    regularise(nx, ny, nz, (void*)jacobian_data, epsilon);

    return Py_None;
}

/*!
 * Method to compute a shape descriptor of a segment.
 */
static PyObject *method_shape_descriptor(PyObject *self, PyObject *args)
{
    FLOATING result;
    FLOATING centroid[3];
    FLOATING spacing[3];
    PyArrayObject *image = NULL;
    PyObject *centroid_tuple = NULL, *spacing_tuple = NULL;
    const char *descriptor;

    (void) self;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, PY_FMT_SHAPE_DESCRIPTOR, &image, &centroid_tuple, &spacing_tuple, &descriptor)) {
        return NULL;
    }

    /* Sanity check */
    if (strcmp(descriptor, CUBENESS) && strcmp(descriptor, OCTAHEDRONESS)) {
        PyErr_Format(PyExc_ValueError, "Invalid shape descriptor '%s'.", descriptor);
        return NULL;
    }

    /* Error check */
    if (!image || !centroid_tuple || !spacing_tuple) {
        return NULL;
    }

    /* Get the centroid */
    if (read_triplet(centroid_tuple, centroid, 3, "centroid")) {
        return NULL;
    }

    /* Get the spacing */
    if (read_triplet(spacing_tuple, spacing, 3, "spacing")) {
        return NULL;
    }

    /* Get the shape */
    const int nz = PyArray_DIM(image, 0);
    const int ny = PyArray_DIM(image, 1);
    const int nx = PyArray_DIM(image, 2);

    /* Get the data pointers */
    FLOATING *image_data = (FLOATING*) PyArray_DATA(image);

    /* Call the library function to compute the shape descriptor */
    if (!strcmp(descriptor, CUBENESS)) {
        result = cubeness((void*) image_data,
                          nx, ny, nz,
                          centroid[0], centroid[1], centroid[2],
                          spacing[0], spacing[1], spacing[2]);
    }
    else if (!strcmp(descriptor, OCTAHEDRONESS)) {
        result = octahedroness((void*) image_data,
                               nx, ny, nz,
                               centroid[0], centroid[1], centroid[2],
                               spacing[0], spacing[1], spacing[2]);
    }

    return Py_BuildValue(PY_FMT_FLOATING, result);
}

/*!
 * Method to compute the Jacobian of a given displacement field.
 */
static PyObject *method_jacobian(PyObject *self, PyObject *args)
{
    FLOATING spacing[3];
    PyArrayObject *field = NULL, *jacobian = NULL;
    PyObject *spacing_tuple = NULL;

    (void) self;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, PY_FMT_JACOBIAN, &spacing_tuple, &field, &jacobian)) {
        return NULL;
    }

    /* Error check */
    if (!field || !jacobian || !spacing_tuple) {
        return NULL;
    }

    /* Get the spacing */
    if (read_triplet(spacing_tuple, spacing, 3, "spacing")) {
        return NULL;
    }

    /* Get the shape (the ITK index order is zyx)*/
    const int nz = PyArray_DIM(jacobian, 0);
    const int ny = PyArray_DIM(jacobian, 1);
    const int nx = PyArray_DIM(jacobian, 2);

    /* Get the data pointers */
    FLOATING *field_data = (FLOATING*) PyArray_DATA(field);
    FLOATING *jacobian_data = (FLOATING*) PyArray_DATA(jacobian);

    /* Call the library function to compute the Jacobian */
    jacobian_dynamic(nx, ny, nz, spacing[0], spacing[1], spacing[2], field_data, jacobian_data);

    return Py_None;
}

/*!
 * Method to compute a displacement with given Jacobian.
 */
static PyObject *method_displacement(PyObject *self, PyObject *args)
{
    FLOATING spacing[3];
    PyObject *spacing_tuple = NULL;
    PyArrayObject *field = NULL, *jacobian = NULL, *mask = NULL;
    const FLOATING epsilon, tolerance, eta, alpha, beta, gamma, delta, zeta;
    const bool strict;
    const long it_max;
    const char *algorithm;

    (void) self;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args,
                          PY_FMT_DISPLACEMENT,
                          &spacing_tuple,
                          &jacobian,
                          &mask,
                          &epsilon,
                          &tolerance,
                          &eta,
                          &alpha,
                          &beta,
                          &gamma,
                          &delta,
                          &zeta,
                          &strict,
                          &it_max,
                          &field,
                          &algorithm)) {
        return NULL;
    }

    /* Sanity check */
    if (strcmp(algorithm, ALGORITHM_GRADIENT) &&
        strcmp(algorithm, ALGORITHM_GREEDY) &&
        strcmp(algorithm, ALGORITHM_MATCHING))
    {
        PyErr_Format(PyExc_ValueError, "Invalid algorithm '%s'.", algorithm);
        return NULL;
    }

    /* Error check */
    if (!jacobian || !mask || !field || !spacing_tuple) {
        return NULL;
    }

    /* Get the spacing */
    if (read_triplet(spacing_tuple, spacing, 3, "spacing")) {
        return NULL;
    }

    /* Get the shape */
    const int nz = PyArray_DIM(jacobian, 0);
    const int ny = PyArray_DIM(jacobian, 1);
    const int nx = PyArray_DIM(jacobian, 2);

    /* Get the data pointers */
    FLOATING *jacobian_data = (FLOATING*) PyArray_DATA(jacobian);
    bool *mask_data = (bool*) PyArray_DATA(mask);
    FLOATING *field_data = (FLOATING*) PyArray_DATA(field);

    /* Call the library function to compute the displacement */
    if (!strcmp(algorithm, ALGORITHM_GRADIENT)) {
        generate_displacement_gradient(
                nx, ny, nz,
                spacing[0], spacing[1], spacing[2],
                (void*)jacobian_data,
                (void*)mask_data,
                epsilon,
                tolerance,
                eta,
                alpha,
                beta,
                gamma,
                delta,
                zeta,
                strict,
                it_max,
                (void*)field_data);
    }
    else if (!strcmp(algorithm, ALGORITHM_GREEDY)) {
        generate_displacement_greedy(
                nx, ny, nz,
                spacing[0], spacing[1], spacing[2],
                (void*)jacobian_data,
                (void*)mask_data,
                epsilon,
                tolerance,
                eta,
                alpha,
                beta,
                gamma,
                delta,
                zeta,
                strict,
                it_max,
                (void*)field_data);
    }
    else if (!strcmp(algorithm, ALGORITHM_MATCHING)) {
        volume_matching_3d(
                nx, ny, nz,
                spacing[0], spacing[1], spacing[2],
                (void*)jacobian_data,
                (void*)mask_data,
                epsilon,
                tolerance,
                eta,
                alpha,
                beta,
                gamma,
                delta,
                zeta,
                strict,
                it_max,
                (void*)field_data);
    }

    return Py_None;
}

/*!
 * Module state handling.
 */

struct module_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

static PyObject *
error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static int myextension_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int myextension_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

/*!
 * Methods exported by the extension.
 */

static char module_docstring[] =
    "Python C extension module.\n"
    "\n"
    "This module wraps the C library for the generation of displacement fields.";

static char docstring_get_float_type_size[] =
    "Get the size of the floating point type.\n"
    "\n"
    "Returns: int";

static char docstring_regularise[] =
    "Apply a lower threshold to the Jacobian.\n"
    "\n"
    "Note: the Jacobian argument will be modified in place.\n"
    "\n"
    "Parameters\n"
    "    Volume image of the Jacobian (numpy.ndarray of type float and dimension 3, indices zyx)\n"
    "    Lower threshold epsilon (float)\n"
    "\n"
    "Returns\n"
    "    None";

static char docstring_shape_descriptor[] =
    "Compute a shape descriptor.\n"
    "\n"
    "Parameters\n"
    "    Binary image encoded as float (numpy.ndarray of type float and dimension 3, indices zyx)\n"
    "    Centroid (tuple of three floats)\n"
    "    Spacing (tuple of three floats)\n"
    "    Descriptor (str with value 'cubeness' or 'octahedroness')\n"
    "\n"
    "Returns\n"
    "    Value for the shape descriptor (float)";

static char docstring_jacobian[] =
    "Compute the Jacobian of a displacement field.\n"
    "\n"
    "Note: the Jacobian argument will be modified in place.\n"
    "\n"
    "Parameters"
    "    Spacing (tuple of three floats)\n"
    "    Input displacement field (numpy.ndarray of type float and dimension 4, indices dzyx)\n"
    "    Output Jacobian (numpy.ndarray of type float and dimension 3, indices zyx)\n"
    "\n"
    "Returns\n"
    "    None";

static char docstring_displacement[] =
    "Compute a displacement with given Jacobian\n."
    "\n"
    "Note: the field argument will be modified in place.\n"
    "\n"
    "Parameters"
    "    Spacing (tuple of three floats)\n"
    "    Input Jacobian (numpy.ndarray of type float and dimension 3, indices zyx)\n"
    "    Input mask (numpy.ndarray of type bool and dimension 3, indices zyx)\n"
    "    Epsilon parameter (float)\n"
    "    Tolerance parameter (float)\n"
    "    Eta parameter (float)\n"
    "    Alpha parameter (float)\n"
    "    Beta parameter (float)\n"
    "    Gamma parameter (float)\n"
    "    Delta parameter (float)\n"
    "    Zeta parameter (float)\n"
    "    Strict parameter (bool)\n"
    "    it_max parameter (int)\n"
    "    Output displacement field (numpy.ndarray of type float and dimension 4, indices dzyx)\n"
    "    Algorithm parameter (string with value 'gradient', 'greedy', or 'matching')\n"
    "\n"
    "Returns\n"
    "    None";

static PyMethodDef myextension_methods[] = {
    {"error_out",           (PyCFunction) error_out,                  METH_NOARGS,  NULL},
    {"get_float_type_size", (PyCFunction) method_get_float_type_size, METH_NOARGS,  docstring_get_float_type_size},
    {"regularise",          (PyCFunction) method_regularise,          METH_VARARGS, docstring_regularise},
    {"shape_descriptor",    (PyCFunction) method_shape_descriptor,    METH_VARARGS, docstring_shape_descriptor},
    {"jacobian",            (PyCFunction) method_jacobian,            METH_VARARGS, docstring_jacobian},
    {"displacement",        (PyCFunction) method_displacement,        METH_VARARGS, docstring_displacement},
    {NULL, NULL, 0, NULL} // sentinel method
};

/*!
 * Module definition attributes.
 */
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_disptools",
        module_docstring,
        sizeof(struct module_state),
        myextension_methods,
        NULL,
        myextension_traverse,
        myextension_clear,
        NULL
};

/*!
 * Module initialisation function.
 */
PyMODINIT_FUNC
PyInit__disptools(void)
{
    PyObject *module = PyModule_Create(&moduledef);

    if (module == NULL) {
        return NULL;
    }

    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("_disptools.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        return NULL;
    }

    return module;
}

