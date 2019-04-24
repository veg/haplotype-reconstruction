import unittest
import os
import itertools as it

import pysam
import pandas as pd
import numpy as np
from Bio import SeqIO
from joblib import Memory
from tqdm import tqdm

from py import ErrorCorrection
from py import partial_covariation_test
from .mock import MockPysamAlignedSegment


data_dir = os.path.join('tests', 'data')
memory = Memory('.joblibcache', verbose=0)


@memory.cache
def get_numpy_fasta():
    fasta_path = os.path.join(data_dir, 'sorted.fasta')
    fasta = SeqIO.parse(fasta_path, 'fasta')
    return np.array([
        list(str(record.seq)) for record in fasta],
        dtype='<U1'
    )


@memory.cache
def get_fasta_counts():
    numpy_fasta = get_numpy_fasta()
    counts = pd.DataFrame([
        pd.Series(numpy_fasta[:, i]).value_counts()
        for i in range(numpy_fasta.shape[1])
    ]).fillna(0)
    zeros = lambda character: (counts[character] == 0).astype(np.int)
    counts['zero_cols'] = zeros('A') + zeros('C') + zeros('G') + zeros('T')
    counts['interesting'] = counts['zero_cols'] < 3
    return counts


@memory.cache
def get_pairs_to_test(threshold):
    print('Generating matrix of pairs...')
    numpy_fasta = get_numpy_fasta()
    counts = get_fasta_counts()
    interesting_sites = counts.index[counts['interesting']]
    number_of_interesting_sites = len(interesting_sites)
    total_pairs = number_of_interesting_sites*(number_of_interesting_sites+1)/2
    pairs = []
    all_pairs = it.combinations(interesting_sites, 2)
    pbar = tqdm(total=total_pairs)
    for i, pair in enumerate(all_pairs):
        site_i, site_j = pair
        content_i = numpy_fasta[:, site_i] != '-'
        content_j = numpy_fasta[:, site_j] != '-'
        valid = content_i & content_j
        if valid.sum() > threshold:
            pairs.append((site_i, site_j))
        pbar.update()
    pbar.close()
    pairs_df = pd.DataFrame({
        'site_i': pd.Series([pair[0] for pair in pairs], dtype='int'),
        'site_j': pd.Series([pair[1] for pair in pairs], dtype='int')
    })
    return pairs_df


class TestErrorCorrection(unittest.TestCase):
    threshold = 20
    test_data = os.path.join(data_dir, 'sorted.bam')

    def test_single_read_count_data(self):
        segment = MockPysamAlignedSegment(
            'CTGATCGCTAACTA', # read 8
            [(0, 4), (1, 1), (0, 6), (2, 1), (0, 3)],
            2
        )
        desired_sequence = np.array(list('CTGACGCTAA-CTA'), dtype='<U1')
        desired_positions = np.arange(2, 16)
        sequence, positions = ErrorCorrection.read_count_data(segment)

        sequences_agree = (sequence == desired_sequence).all()
        self.assertTrue(sequences_agree)

        positions_agree = (positions == desired_positions).all()
        self.assertTrue(positions_agree)

    def test_real_nucleotide_counts(self):
        bam = pysam.AlignmentFile(self.test_data, 'rb')
        error_correction = ErrorCorrection(bam)
        bam_counts = error_correction.get_nucleotide_counts()
        desired_columns = ['A', 'C', 'G', 'T', 'interesting']

        fasta_counts = get_fasta_counts()

        fasta_subset = fasta_counts.loc[:, desired_columns]
        bam_subset = bam_counts.loc[:, desired_columns]
        fasta_equals_bam = fasta_subset == bam_subset
        self.assertTrue(fasta_equals_bam.all().all())

    def test_partial_covariation_test(self):
        bam = pysam.AlignmentFile(self.test_data, 'rb')
        error_correction = ErrorCorrection(bam)
        pairs = error_correction.get_pairs()
        partial_covariation_test((bam.filename, pairs[:1000], 1)) 

    def test_full_covariation_test(self):
        bam = pysam.AlignmentFile(self.test_data, 'rb')
        error_correction = ErrorCorrection(bam)
        covarying_sites = error_correction.full_covariation_test()
        print(covarying_sites+1)

    def test_write_corrected_reads(self):
        corrected_bam_filename = 'corrected.bam'
        bam = pysam.AlignmentFile(self.test_data, 'rb')
        error_correction = ErrorCorrection(bam)
        error_correction.write_corrected_reads(corrected_bam_filename)
        os.remove(corrected_bam_filename)

