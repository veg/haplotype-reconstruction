import unittest
import re
import os

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from py import SAMFASTAConverter
from py import AlignedSegment
from .mock import MockPysamAlignedSegment
from .mock import excel_information as excel
from .mock import create_sample_alignment
from .mock import excel_filepath


class TestSingleSAMFASTAConverter(unittest.TestCase):

    def setUp(self):
        self.sfc = SAMFASTAConverter()

    def test_get_numeric_representation(self):
        seq = Seq('ACGGCAT')
        record = SeqRecord(seq)
        result = self.sfc.get_numeric_representation(record)
        desired_result = np.array([0, 1, 2, 2, 1, 0, 3], dtype=np.int)

    def test_single_segment_to_fasta(self):
        segment = MockPysamAlignedSegment(
            'ACTCCTCGAA',
            [(0, 1), (1, 1), (2, 1), (0, 5), (2, 1), (0, 3)]
        )
        fasta = self.sfc.single_segment_to_fasta(segment)
        desired_fasta = np.array(list('AC-TCCTC-GAA'), dtype='<U1')
        self.assertTrue(np.all(fasta==desired_fasta), 'with insertion')

        fasta = self.sfc.single_segment_to_fasta(segment, False)
        desired_fasta = np.array(list('A-TCCTC-GAA'), dtype='<U1')
        self.assertTrue(np.all(fasta==desired_fasta), 'without insertion')


class TestMultipleSAMFASTAConverter(unittest.TestCase):

    def setUp(self):
        self.sam_fasta_converter = SAMFASTAConverter()
        self.mock_segments = [
            MockPysamAlignedSegment(
                'ATCTCACGATTGC', # read 1
                [(0, 13)],
                0
            ),
            MockPysamAlignedSegment(
                'CGCAGATAACCT', # read 5
                [(0, 1), (2, 1), (0, 11)],
                2
            ),
            MockPysamAlignedSegment(
                'CTGATCGCTAACTA', # read 8
                [(0, 4), (1, 1), (0, 6), (2, 1), (0, 3)],
                2
            ),
            MockPysamAlignedSegment(
                'CGCCGATGACTGCTA', # read 11
                [(0, 10), (1, 2), (0, 3)],
                3
            ),
            MockPysamAlignedSegment(
                'GACGAGAACACTAAC', # read 14
                [(0, 9), (1, 1), (0, 5)],
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
        self.assertEqual(segments[3].current_character, '~')
        self.assertEqual(segments[4].current_character, '~')

    def test_partial_single_case_with_insertions(self):
        sam_fasta_converter = SAMFASTAConverter()
        reference_length = 20
        window_start = 2
        window_end = 16
        sam_fasta_converter.initialize(
            self.mock_segments, reference_length, window_start, window_end
        )

        desired_fasta = [
            'CTCA-CGATTGC~~~~~',
            'C-GC-AGATAAC--CT~',
            'CTGATCGCTAA---CTA',
            '~CGC-CGATGACTGCTA',
            '~~GA-CGAGAACA-CTA'
        ]
        number_of_rows = len(desired_fasta)
        number_of_columns = len(desired_fasta[0])
        for j in range(number_of_columns):
            insertion_result = sam_fasta_converter.handle_insertions()
            if not insertion_result:
                sam_fasta_converter.handle_deletions_and_matches()
            column = sam_fasta_converter.fasta[:, j]
            desired_column = [row[j] for row in desired_fasta]
            triplets = zip(range(number_of_rows), column, desired_column)
            for i, character, desired_character in triplets:
                error_message = 'row %d, column %d' % (i, j)
                self.assertEqual(character, desired_character, error_message)

    def test_full_single_case_with_insertions(self):
        sample_alignment = create_sample_alignment()
        df = pd.read_excel(excel_filepath, header=None)
        desired_fasta = np.array(
            df.iloc[
                excel['read_row_start']: excel['read_row_end'],
                excel['read_window_start'] : excel['read_window_end']
            ],
            dtype='<U1'
        )
        desired_reference_position = np.array(df.iloc[
            excel['reference_position_row'],
            excel['read_window_start']: excel['read_window_end']
        ], dtype=np.int)

        sfc = SAMFASTAConverter()
        reference_length = 20
        window_start = 2
        window_end = 16
        fasta, reference_position = sfc.sam_window_to_fasta_with_insertions(
            sample_alignment.fetch(), reference_length,
            window_start, window_end
        )
        self.assertTrue(np.all(fasta == desired_fasta))
        self.assertTrue(np.all(
            reference_position == desired_reference_position
        ))

    def test_full_single_case_without_insertions(self):
        sample_alignment = create_sample_alignment()
        df = pd.read_excel(excel_filepath, header=None)
        window_start = 2
        window_end = 16
        columns = np.concatenate([
            np.arange(
               excel['read_col_start'],
               excel['insertion_columns'][0],
            ),
            np.arange(
               excel['insertion_columns'][0] + 1,
               excel['insertion_columns'][1]
            ),
            np.arange(
               excel['insertion_columns'][1] + 1,
               excel['insertion_columns'][2]
            ),
            np.arange(
               excel['insertion_columns'][2] + 1,
               excel['read_col_end'],
            )
        ])
        full_fasta = np.array(
            df.iloc[
                excel['read_row_start']: excel['read_row_end'],
                columns
            ],
            dtype='<U1'
        )
        desired_fasta = full_fasta[:, window_start: window_end]

        sfc = SAMFASTAConverter()
        reference_length = 20
        fasta = sfc.sam_window_to_fasta(
            sample_alignment, window_start, window_end
        )
        self.assertTrue(np.all(fasta == desired_fasta))

