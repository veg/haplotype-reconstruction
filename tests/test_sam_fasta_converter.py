import unittest

import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from py import SAMFASTAConverter
from py import AlignedSegment


class MockPysamAlignedSegment:
    def __init__(self, alignment_sequence, cigartuples, reference_start=None):
        self.query_alignment_sequence = alignment_sequence
        self.cigartuples = cigartuples
        self.reference_start = reference_start


class TestSAMFASTAConverter(unittest.TestCase):

    def setUp(self):
        self.sam_fasta_converter = SAMFASTAConverter()

    def test_get_numeric_representation(self):
        seq = Seq('ACGGCAT')
        record = SeqRecord(seq)
        result = self.sam_fasta_converter.get_numeric_representation(record)
        desired_result = np.array([0, 1, 2, 2, 1, 0, 3], dtype=np.int)

    def test_single_aligned_segment_to_fasta(self):
        segment = MockPysamAlignedSegment(
            'ACTCCTCGAA',
            [(0, 1), (1, 1), (2, 1), (0, 5), (2, 1), (0, 3)]
        )
        fasta = self.sam_fasta_converter.single_aligned_segment_to_fasta(segment)
        desired_fasta = np.array(list('AC-TC-GAACTC'), dtype='<U1')
        self.assertTrue(np.all(fasta==desired_fasta))

    def test_initiate_conversion(self):
        mock_segments = [
            MockPysamAlignedSegment('ATCCGACGATTAC', [(0, 13)],
                0
            ),
            MockPysamAlignedSegment(
                'TCGCCGATAACGCT',
                [(0, 1), (2, 1), (0, 12)],
                2
            ),
            MockPysamAlignedSegment(
                'GACGATAACAGCTAAC',
                [(0, 9), (1, 1), (0, 6)],
                4
            )
        ]
        segments = [AlignedSegment(segment) for segment in mock_segments]
        for segment in segments:
            segment.initiate_conversion(2)
            if segment.active:
                self.assertEqual(segments[0].position_along_reference, 2)
        self.assertFalse(segments[2].active)

