//
// Created by redwards on 8/4/20.
//
#include <Python.h>
#include "pyprimer-trimming.h"
#include "primer-trimming.h"
#include "pyprinseq.h"

PyObject *
pyprimer_trimming(PyObject *self, PyObject *args) {
    // COMMAND LINE OPTIONS
    char *infile = NULL;
    char **primersL = NULL;
    char **primersR = NULL;

    char *leftPrimer = NULL;
    char *rightPrimer = NULL;

    if(!PyArg_ParseTuple(args, "sss", &infile, &leftPrimer, &rightPrimer)) {
        PyErr_SetString(PyExc_RuntimeError, "Could not parse the arguments to python_input");
        return NULL;
    }

    if (leftPrimer != NULL)
        primersL = load_primers(leftPrimer);
    if (rightPrimer != NULL)
        primersR = load_primers(rightPrimer);

    if ((primersL == NULL) & (primersR == NULL)) {
        exit(EXIT_FAILURE);
    }

    int ro = trim_primers(infile, primersL, primersR);

    return PyLong_FromLong((long) ro);
}

