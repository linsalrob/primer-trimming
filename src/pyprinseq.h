// Common function definitions for all the prinseq python modules. These are basically
// all the PyObject exports
// Created by redwards on 8/4/20.
//

#ifndef PYPRINSEQ_H
#define PYPRINSEQ_H

#include "pyprimer-predictions.h"
#include "pyprimer-trimming.h"


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



#endif //PYPRINSEQ_H
