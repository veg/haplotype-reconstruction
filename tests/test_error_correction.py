import unittest
import os

import pysam
import pandas as pd
import numpy as np
from Bio import SeqIO

from py import ErrorCorrection
from .mock import MockPysamAlignedSegment


class TestErrorCorrection(unittest.TestCase):
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
        data_dir = os.path.join('tests', 'data')
        bam_path = os.path.join(data_dir, 'sorted.bam')
        bam = pysam.AlignmentFile(bam_path, 'rb')
        error_correction = ErrorCorrection(bam)
        fasta_path = os.path.join(data_dir, 'sorted.fasta')
        fasta = SeqIO.parse(fasta_path, 'fasta')
        numpy_fasta = np.array([
            list(str(record.seq)) for record in fasta],
            dtype='<U1'
        )
        fasta_counts = pd.DataFrame([
            pd.Series(numpy_fasta[:, i]).value_counts()
            for i in range(numpy_fasta.shape[1])
        ]).fillna(0)
        bam_counts = error_correction.nucleotide_counts()
        nucleotide_columns = ['A', 'C', 'G', 'T']
        fasta_subset = fasta_counts.loc[:, nucleotide_columns]
        bam_subset = bam_counts.loc[:, nucleotide_columns]
        equality = fasta_subset == bam_subset
        self.assertTrue(equality.all().all())

