//
// Created by redwards on 8/4/20.
//
#include <Python.h>
#include "pyprimer-trimming.h"
#include "trimprimers.h"
#include "pyprinseq.h"

PyObject *
pyprimer_trimming(PyObject *self, PyObject *args) {
    // COMMAND LINE OPTIONS
    char *infile = NULL;
    char **primersL = NULL;
    char **primersR = NULL;

    char *leftPrimer = NULL;
    char *rightPrimer = NULL;

    if(!PyArg_ParseTuple(args, "szz", &infile, &leftPrimer, &rightPrimer)) {
        PyErr_SetString(PyExc_RuntimeError, "Could not parse the arguments to python_input");
        return NULL;
    }

    if (0) {
        fprintf(stderr, "infile: %s\n", infile);
        fprintf(stderr, "left primer file: %s\n", leftPrimer);
        fprintf(stderr, "right primer file: %s\n", rightPrimer);
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

