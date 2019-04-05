import numpy as np


class NumericSequence:
    def __init__(self):
        self.iupac_nucleotides = list('ACGTRYSWKMBDHVN-')
        self.n_nucleotides = len(self.iupac_nucleotides)
        self.embedding = np.array([
            [1,   0,   0,    0,   0], # A
            [0,   1,   0,    0,   0], # C
            [0,   0,   1,    0,   0], # G
            [0,   0,   0,    1,   0], # T
            [.5,  0,  .5,    0,   0], # R
            [0,   .5,  0,    .5,  0], # Y
            [0,   .5,  .5,   0,   0], # S
            [.5,  0,   0,    .5,  0], # W
            [0,   0,   .5,   .5,  0], # K
            [.5,  .5,  0,    0,   0], # M
            [0,   1/3, 1/3,  1/3, 0], # B
            [1/3, 0,   1/3,  1/3, 0], # D
            [1/3, 1/3, 0,    0/3, 0], # H
            [1/3, 1/3, 1/3,  0,   0], # V
            [.25, .25, .25, .25,  0], # N
            [0,   0,   0,    0,   1], # -
        ])

