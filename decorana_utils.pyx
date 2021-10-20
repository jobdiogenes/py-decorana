from numpy cimport ndarray as ar

cdef extern from "pydecorana_utils.h":
    void c_cutup(double* x, int* ix,int* mi, int* mk)
    void c_eigy(double*, double*, double*, int*, int*, int*, double*,
		    int*, int*, int*, int*, int*, int*, int*, double*,
		    double*, double*, double*, double*, double*, double*,
		    double*, int*, int*, int*, double*, double*)
    void c_yxmult(double*, double*, int*, int*, int*, int*, int*,
		      int*, double*)

def cutup(double x, int ix, int mi,int mk):
    cdef ar[double, mode="c"] 