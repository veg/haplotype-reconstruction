import itertools as it

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import pysam


class ErrorCorrection:
    def __init__(self, pysam_alignment):
        self.pysam_alignment = pysam_alignment
        self.reference_length = pysam_alignment.header['SQ'][0]['LN']
        self.number_of_reads = 0
        for read in pysam_alignment.fetch():
            self.number_of_reads += 1

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

        df = pd.DataFrame(counts, columns=characters)
        zeros = lambda character: (df[character] == 0).astype(np.int)
        zero_cols = zeros('A') + zeros('C') + zeros('G') + zeros('T')
        df['interesting'] = zero_cols < 3
        return df
    
    def pairs_to_test(self, threshold):
        counts = self.nucleotide_counts()
        interesting_sites = np.array(counts.index[counts.interesting])
        reads = self.pysam_alignment.fetch()
        reference_starts = np.array(
            [read.reference_start for read in reads],
            dtype=np.int
        )
        reads = self.pysam_alignment.fetch()
        reference_stops = np.array(
            [read.reference_end for read in reads],
            dtype=np.int
        )
        interesting_starts = np.searchsorted(
            interesting_sites,
            reference_starts
        )
        interesting_stops = np.searchsorted(
            interesting_sites,
            reference_stops
        )
        all_pairs = it.chain.from_iterable((
            it.combinations(interesting_sites[start:stop], 2)
            for start, stop in zip (interesting_starts, interesting_stops)
        ))
        columns = ['site_i', 'site_j']
        all_counts = pd.DataFrame(
            all_pairs, columns=columns
        ).groupby(
            columns
        ).size().reset_index(name='count')

        result = all_counts[
            all_counts['count'] > 20
        ].sort_values(
            by=columns
        )[columns].reset_index(drop=True)

        return result

