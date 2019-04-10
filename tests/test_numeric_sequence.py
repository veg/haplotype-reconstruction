import unittest

import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from py import NumericSequence


class MockPysamAlignedSegment:
    def __init__(self, alignment_sequence, cigartuples):
        self.query_alignment_sequence = alignment_sequence
        self.cigartuples = cigartuples

class TestNumericSequence(unittest.TestCase):

    def setUp(self):
        self.numeric_sequence = NumericSequence()
        self.first_segment = MockPysamAlignedSegment(
            'ACTCGAACTC',
            [(0, 1), (1, 1), (2, 1), (0, 2), (2, 1), (1, 1), (0, 5)]
        )
        self.second_segment = MockPysamAlignedSegment(
            'ACTCCTCGAA',
            [(0, 1), (1, 1), (2, 1), (0, 5), (2, 1), (0, 3)]
        )

    def test_get_numeric_representation(self):
        seq = Seq('ACGGCAT')
        record = SeqRecord(seq)
        result = self.numeric_sequence.get_numeric_representation(record)
        desired_result = np.array([0, 1, 2, 2, 1, 0, 3], dtype=np.int)

    def test_aligned_segment_to_fasta(self):
        first_fasta = self.numeric_sequence.aligned_segment_to_fasta(self.first_segment)
        desired_first_fasta = np.array(list('AC-TC-GAACTC'), dtype='<U1')
        self.assertTrue(np.all(first_fasta==desired_first_fasta))

        second_fasta = self.numeric_sequence.aligned_segment_to_fasta(self.second_segment)
        desired_second_fasta = np.array(list('AC-TCCTC-GAA'), dtype='<U1')
        self.assertTrue(np.all(second_fasta == desired_second_fasta))

    def test_embed_sam_in_window(self):
        fasta = self.numeric_sequence.embed_sam_in_window(
            [self.first_segment, self.second_segment],
            (0, 10)
        )
        desired_fasta = np.array([
            list('AC-TC-GAAC-TC'),
            list('AC-TCC-TC-GAA')
        ], dtype='<U1')
        self.assertTrue(np.all(fasta == desired_fasta))

