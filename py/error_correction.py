import itertools as it
from multiprocessing import Pool

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import pysam

from .sam_fasta_converter import SAMFASTAConverter


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return it.zip_longest(*args, fillvalue=fillvalue)


def perform_partial_covariation_test(arguments):
    pysam_path, pairs, i = arguments
    pairs, pairs_for_max, pairs_for_min = it.tee(
        it.filterfalse(lambda x: x is None, pairs),
        3
    )
    print('   ...processing block %d...' % i)
    pysam_alignment = pysam.AlignmentFile(pysam_path, 'rb')
    pair_min = min([min(pair[0], pair[1]) for pair in pairs_for_min])
    pair_max = max([max(pair[0], pair[1]) for pair in pairs_for_max])
    sfc = SAMFASTAConverter()
    fasta_window = sfc.sam_window_to_fasta(pysam_alignment, pair_min, pair_max+1)
    results = []
    for col_i, col_j in pairs:
        idx_i = col_i - pair_min
        idx_j = col_j - pair_min
        i_char_counts = pd.Series(
            fasta_window[:, idx_i]
        ).value_counts().drop('~', errors='ignore')
        i_chars = i_char_counts.index[i_char_counts > 0]

        j_char_counts = pd.Series(
            fasta_window[:, idx_j]
        ).value_counts().drop('~', errors='ignore')
        j_chars = j_char_counts.index[j_char_counts > 0]

        content_i = fasta_window[:, idx_i] != '~'
        content_j = fasta_window[:, idx_j] != '~'
        valid = content_i & content_j
        if valid.sum() == 0:
            continue
        for i_char, j_char in it.product(i_chars, j_chars):
            equals_i = fasta_window[valid, idx_i] == i_char
            equals_j = fasta_window[valid, idx_j] == j_char
            X_11 = (equals_i & equals_j).sum()
            X_21 = (~equals_i & equals_j).sum()
            X_12 = (equals_i & ~equals_j).sum()
            X_22 = (~equals_i & ~equals_j).sum()
            
            table = [
                [X_11 , X_12],
                [X_21 , X_22]
            ]
            _, p_value = fisher_exact(table)
            results.append((col_i, col_j, i_char, j_char, p_value))
    pysam_alignment.close()
    print('   ...done block %d!' % i)
    return pd.DataFrame(
        results, columns=('col_i', 'col_j', 'i_char', 'j_char', 'p_value')
    )


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
        print('Obtaining nucleotide counts...')
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

    def get_pairs(self):
        counts = self.nucleotide_counts()
        interesting = counts.index[counts.interesting]
        max_read_length = max([
            read.infer_query_length()
            for read in self.pysam_alignment.fetch()
        ])
        pairs = list(filter(
            lambda pair: pair[1] - pair[0] <= max_read_length,
            it.combinations(interesting, 2)
        ))
        return pairs
    
    def perform_full_covariation_test(self, threshold=20, stride=10000,
            ncpu=24, block_size=250, fdr=.001):
        pairs = self.get_pairs()
        arguments = [
            (self.pysam_alignment.filename, group, i)
            for i, group in enumerate(grouper(pairs, block_size))
        ]
        n_pairs = len(pairs)
        n_blocks = len(arguments)
        print('Processing %d blocks of %d pairs...' % (n_blocks, n_pairs))
        pool = Pool(processes=ncpu)
        result_dfs = pool.map(perform_partial_covariation_test, arguments)
        result_df = pd.concat(result_dfs).sort_values(by='p_value')
        pool.close()
        print('...done!')
        print('Performing multiple testing correction...')
        m = len(result_df)
        result_df['bh'] = result_df['p_value'] <= fdr*np.arange(1, m+1)/m
        cutoff = (1-result_df['bh']).to_numpy().nonzero()[0][0]
        covarying_sites = np.unique(
            np.concatenate([
                result_df['col_i'].iloc[:cutoff],
                result_df['col_j'].iloc[:cutoff]
            ])
        )
        covarying_sites.sort()
        return covarying_sites

    def __del__(self):
        self.pysam_alignment.close()

