from scipy.io import mmread
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import norm
import numpy as np

"""
Simple test script to compare the results of my class to.
Used to check the construction, matrix-vector multiplication
as well as the norm computation.
"""

if __name__ == "__main__":
    # matrix = mmread("small_example.mtx")
    matrix = mmread("lnsp_131.mtx")

    print("CSR-MATRIX")
    csr = csr_matrix(matrix)
    # to_multiply = np.array([42, -17, 100, 0, 73, -5, 21, 8, -33, 55])
    # print(csr @ to_multiply)
    # print(csr)
    # print(csr.data)
    # print(csr.indices)
    # print(csr.indptr)

    # print("CSC-MATRIX")
    # csc = csc_matrix(matrix)
    # print(csc)
    # print(csc.data)
    # print(csc.indices)
    # print(csc.indptr)

    print("Norm")
    for ord in ["fro", np.inf, 1]:
        print("{} = {}".format(ord, norm(csr, ord=ord)))
    # print(csr.data)
    # print(np.sum(csr.data))
    # print(np.sqrt(np.sum(csr.data)))
