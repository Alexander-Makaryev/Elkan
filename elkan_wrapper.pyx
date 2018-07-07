cimport numpy as np
import numpy as np
cimport cython
from libc.stdlib cimport malloc, free
from libcpp cimport bool
import sys

cdef extern from "elkan.cpp":
    void elkan(int p, int n, double** x, int k, double** cStart, double threshold, int maxIterationsCount, double** c, int* pointsCenters, int* iterationsCount, bool* isEmptyCluster)

def execute(int p, int n, x, int k, cStart, threshold, int maxIterationsCount):
    cdef np.ndarray[double, ndim=1, mode="c"] temp
    c = np.empty([k, p], dtype=np.double)
    cdef double **c_c = <double **> malloc(k*sizeof(double *))
    for i in range(k):
        temp = c[i]
        c_c[i] = &temp[0]
    cdef np.ndarray[int, ndim=1, mode="c"] pointsCenters = np.empty(n, dtype=np.int)
    cdef int iterationsCount
    cdef bool isEmptyCluster

    
    cdef double **c_x = <double **> malloc(n*sizeof(double *))
    for i in range(n):
        temp = x[i]
        c_x[i] = &temp[0]
    
    cdef double **c_cStart = <double **> malloc(k*sizeof(double *))
    for i in range(k):
        temp = cStart[i]
        c_cStart[i] = &temp[0]
    
    elkan(p, n, c_x, k, c_cStart, threshold, maxIterationsCount, c_c, &pointsCenters[0], &iterationsCount, &isEmptyCluster)
    
    free(c_x)
    free(c_c)
    free(c_cStart)
    return (c, pointsCenters, iterationsCount, isEmptyCluster)
    
    