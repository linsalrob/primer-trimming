//
// Created by redwards on 8/4/20.
//
#include <Python.h>
#include "pyprinseq.h"

PyMethodDef PyPrinseqMethods[] = {
        {"primerpredict", primer_predictions, METH_VARARGS, "Python interface for ANSI-C primer predictions"},
        {"primertrimming", pyprimer_trimming, METH_VARARGS, "Python interface for ANSI-C primer trimming"},
        {NULL, NULL, 0, NULL}
};

struct PyModuleDef PyPrinseqModule = {
        PyModuleDef_HEAD_INIT,
        "pyprinseq",
        "Python interface for ANSI-C primer predictions",
        -1,
        PyPrinseqMethods
};
