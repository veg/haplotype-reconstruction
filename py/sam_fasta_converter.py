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
        self.nuc2ind = {
            nuc: i for i, nuc in enumerate(self.iupac_nucleotides)
        }

    def get_numeric_representation(self, record):
        return np.array(
            [self.nuc2ind[char] for char in record.seq],
            dtype=np.int
        )

    def aligned_segment_to_fasta(self, aligned_segment):
        read_sequence = aligned_segment.query_alignment_sequence
        position = 0
        alignment = ''
        cigar_tuples = aligned_segment.cigartuples
        number_of_cigar_tuples = len(cigar_tuples)
        for i, cigar_tuple in enumerate(cigar_tuples):
            action = cigar_tuple[0]
            stride = cigar_tuple[1]
            match = action == 0
            insertion = action == 1
            deletion = action == 2
            if match or insertion:
                alignment += read_sequence[position: position + stride]
                position += stride
            elif deletion:
                if len(alignment) > 0 and i < number_of_cigar_tuples:
                    alignment += stride * '-'
        full_seq_np = np.array(
            list(alignment),
            dtype='<U1'
        )
        start = np.argmax(full_seq_np != '-')
        stop = len(full_seq_np) - np.argmax(full_seq_np[::-1] != '-')
        return full_seq_np[start:stop]

    def embed_numeric_fasta(numeric_fasta):
        n = numeric_fasta.shape[0]
        k = numeric_fasta.shape[1]
        embedded_fasta = np.reshape(
            embedding[numeric_fasta, :],
            newshape=(n, 5*k)
        )

