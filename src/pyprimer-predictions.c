/*
 * This is source code for integration with Python to make things easy!
 * Created by redwards on 8/4/20.
*/

#include <Python.h>
#include <stdbool.h>
#include "predictprimers.h"
#include "pyprimer-predictions.h"
#include "pyprinseq.h"


PyObject *
primer_predictions(PyObject *self, PyObject *args) {

    char * infile = NULL;
    int kmerlen = 0;
    double minpercent = 20.0;

    bool three_prime = false;
    bool fasta_output = false, print_kmer_counts = false, print_abundance = false;
    bool print_short = false, debug = false;

    /* Parse arguments */
    // , , &fasta_output, &print_kmer_counts, &print_abundance, &print_short, &debug
    if(!PyArg_ParseTuple(args, "sidb", &infile, &kmerlen, &minpercent, &three_prime)) {
            PyErr_SetString(PyExc_RuntimeError, "Could not parse the arguments to python_input");
        return NULL;
    }

    // for the results
    char **allprimers;
    allprimers = malloc(sizeof(*allprimers) * 1); // initializing with 1 member, but will realloc later
    int allprimerposition=0;

    int ro = predict_primers(infile, kmerlen, minpercent, fasta_output, three_prime, print_kmer_counts, print_abundance, print_short, debug, allprimers, &allprimerposition);
    if (ro != 0) {
        fprintf(stderr, "Error: Running the primer search returned %d\n", ro);
        return NULL;
    }

    // convert our list of primers to a python object
    PyObject *result = PyList_New(0);

    for (int i=0; i < allprimerposition; i++) {
        PyObject *item = Py_BuildValue("s", allprimers[i]);
        fprintf(stderr, "Appended %s\n", allprimers[i]);
        PyList_Append(result, item);
    }

    // CRITICAL: We need to reset the repeats before we
    // call this method again from within the same run!
    for (int i=0; i<allprimerposition; i++)
        free(allprimers[i]);

    free(allprimers);
    return result;
}


PyMODINIT_FUNC PyInit_PyPrinseq(void) {
    return PyModule_Create(&PyPrinseqModule);
}
