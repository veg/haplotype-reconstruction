import itertools as it

import numpy as np
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve


def simplex_projection(x):
    n = len(x)
    u = np.sort(x)[::-1]
    v = (np.cumsum(u)-1)/(np.arange(n)+1)
    K = np.argmax(v >= u)-1
    tau = v[K]
    return np.maximum(x-tau, np.zeros(n))


class QuadraticProx:
    def __init__(self, L, y, lambd):
        self.lambd = lambd
        LT = L.transpose()
        self.LTy = np.dot(LT, y)
        I = np.eye(L.shape[1])
        A = I + lambd*np.dot(LT, L)
        self.lu_factors = lu_factor(A)

    def apply(self, x):
        z = x + self.lambd*self.LTy
        return lu_solve(self.lu_factors, z)


def is_valid_superread(superread):
    is_source = superread['id'] == 'source'
    is_target = superread['id'] == 'target'
    if is_source or is_target:
        return False
    return not superread['discarded']


def build_L_and_y(superread_info, candidate_info):
    surviving_superread_indices = list(set(it.chain.from_iterable(candidate_info)))
    surviving_superread_indices.sort()
    reindex_map = {ssi: i for i, ssi in enumerate(surviving_superread_indices)}
    valid_superreads = [
        superread for superread in superread_info
        if is_valid_superread(superread)
    ]
    n = len(surviving_superread_indices)
    m = len(candidate_info)

    L = np.zeros((n, m))
    for i, row_info in enumerate(candidate_info):
        reindexed_row = [reindex_map[j] for j in row_info]
        L[reindexed_row, i] = 1

    y = np.array([
        row_info['frequency'] for row_info in valid_superreads
    ])

    return L, y


def perform_regression(superread_info, candidate_info, maxit=100000,
        tolerance=1e-6, lambd=.1):
    L, y = build_L_and_y(superread_info, candidate_info)
    quadratic_prox = QuadraticProx(L, y, lambd)

    y = np.ones(L.shape[1])/L.shape[1]
    x = simplex_projection(y)
    y = y + lambd*(quadratic_prox.apply(2*x - y) - x)
    for i in range(maxit):
        xold = np.array(x)
        x = simplex_projection(y)
        y = y + lambd*(quadratic_prox.apply(2*x - y) - x)
        relative_error = np.linalg.norm(x-xold)/np.linalg.norm(x) 
        if i % 1000 == 0:
            params = (i, relative_error)
            print('Iteration %d, relative error %.2e...' % params)
        if relative_error < tolerance:
            params = (relative_error, i)
            print('Reached relative error of %.5e at iteration %d!' % params)
            break
    actual_haplotypes = []
    for i, frequency in enumerate(x):
        if frequency > 0:
            actual_haplotypes.append(i)
            print('Candidate %2d: frequency %.3f' % (i, frequency))
    return actual_haplotypes

