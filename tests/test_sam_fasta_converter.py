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


class TestSingleSAMFASTAConverter(unittest.TestCase):

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


class TestMultipleSAMFASTAConverter(unittest.TestCase):

    def setUp(self):
        self.sam_fasta_converter = SAMFASTAConverter()
        self.mock_segments = [
            MockPysamAlignedSegment(
                'ATCTGACGATTAC', # read 1
                [(0, 13)],
                0
            ),
            MockPysamAlignedSegment(
                'CGCCGATAACGCT', # read 5
                [(0, 1), (2, 1), (0, 12)],
                2
            ),
            MockPysamAlignedSegment(
                'CTGATCGATAAGCTA', # read 8
                [(0, 4), (1, 1), (0, 6), (2, 1), (0, 4)],
                2
            ),
            MockPysamAlignedSegment(
                'CGCCGATAACTGCTA', # read 11
                [(0, 10), (1, 1), (0, 4)],
                3
            ),
            MockPysamAlignedSegment(
                'GACGATAACAGCTAAC', # read 14
                [(0, 9), (1, 1), (0, 6)],
                4
            )
        ]

    def test_initiate_conversion(self):
        segments = [AlignedSegment(segment) for segment in self.mock_segments]
        window_start = 2
        for i, segment in enumerate(segments):
            print('Initiating segment %d for conversion...' % i)
            segment.initiate_conversion(window_start)
            if segment.entered:
                self.assertEqual(segment.position_in_reference, window_start)
        self.assertFalse(segments[3].entered)
        self.assertFalse(segments[4].entered)

        self.assertEqual(segments[0].current_character, 'C')
        self.assertEqual(segments[1].current_character, 'C')
        self.assertEqual(segments[2].current_character, 'C')
        self.assertEqual(segments[3].current_character, '-')
        self.assertEqual(segments[4].current_character, '-')

    def test_single_test_case(self):
        sam_fasta_converter = SAMFASTAConverter()
        reference_length = 20
        window_start = 2
        sam_fasta_converter.initialize(self.mock_segments, reference_length, window_start)

        desired_fasta = [
            'CTGA-CGATTAC-----',
            'C-GC-CGATAAC-GCT-',
            'CTGATCGATAA--GCTA',
            '-CGC-CGATAACTGCTA',
            '--GA-CGATAACAGCTA'
        ]
        for j in range(4):
            insertion_result = sam_fasta_converter.handle_insertions()
            self.assertFalse(insertion_result)
            sam_fasta_converter.handle_deletions_and_matches()
            column = sam_fasta_converter.fasta[:, j]
            desired_column = [row[j] for row in desired_fasta]
            print(desired_column)
            triplets = zip(range(len(desired_column)), column, desired_column)
            for i, character, desired_character in triplets:
                print(character, desired_character)
                error_message = 'row %d, column %d' % (i, j)
                self.assertEqual(character, desired_character, error_message)

