//
// Created by redwards on 8/4/20.
//

#ifndef PRIMER_TRIMMING_PYPRIMER_PREDICTIONS_H
#define PRIMER_TRIMMING_PYPRIMER_PREDICTIONS_H


static PyObject * primer_predictions(PyObject *self, PyObject *args);

static PyMethodDef PyPrinseqMethods[] = {
        {"primerpredict", primer_predictions, METH_VARARGS, "Python interface for ANSI-C primer predictions"},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef PyPrinseqModule = {
        PyModuleDef_HEAD_INIT,
        "pyprinseq",
        "Python interface for ANSI-C primer predictions",
        -1,
        PyPrinseqMethods
};



#endif //PRIMER_TRIMMING_PYPRIMER_PREDICTIONS_H
