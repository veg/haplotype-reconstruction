import itertools as it

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import pysam


class ErrorCorrection:
    def __init__(self, pysam_alignment):
        self.pysam_alignment = pysam_alignment
        self.reference_length = pysam_alignment.header['SQ'][0]['LN']

    @staticmethod
    def read_count_data(read):
        sequence_length = np.array([
            cigar_tuple[1]
            for cigar_tuple in read.cigartuples
            if cigar_tuple[0] != 1
        ]).sum()
        first_position = read.reference_start
        last_position = first_position + sequence_length
        positions = np.arange(first_position, last_position)
        segments = []
        number_of_cigar_tuples = len(read.cigartuples)
        unaligned_sequence = read.query_alignment_sequence
        position = 0
        for i, cigar_tuple in enumerate(read.cigartuples):
            action = cigar_tuple[0]
            stride = cigar_tuple[1]
            match = action == 0
            insertion = action == 1
            deletion = action == 2
            if match:
                segments.append(unaligned_sequence[position: position + stride])
                position += stride
            elif insertion:
                position += stride
            elif deletion:
                if len(segments) > 0 and i < number_of_cigar_tuples:
                    segments.append(stride * '-')
        sequence = np.concatenate([list(segment) for segment in segments])
        return sequence, positions

    def nucleotide_counts(self):
        characters = ['A', 'C', 'G', 'T', '-']
        counts = np.zeros((self.reference_length, 5))
        for read in self.pysam_alignment.fetch():
            sequence, positions = self.read_count_data(read)
            for character_index, character in enumerate(characters):
                rows = positions[sequence == character]
                counts[rows, character_index] += 1
        return pd.DataFrame(counts, columns=characters)

