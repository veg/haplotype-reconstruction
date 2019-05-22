import unittest
import re
import os

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from py import MappedReads
from py import AlignedSegment
from .mock import MockPysamAlignedSegment
from .mock import MockPysamAlignment
from .mock import excel_information as excel
from .mock import create_sample_alignment
from .mock import excel_filepath


class TestAlignedSegment(unittest.TestCase):

    def test_single_segment_to_fasta(self):
        mock_segment = MockPysamAlignedSegment(
            'ACTCCTCGAA',
            [(0, 1), (1, 1), (2, 1), (0, 5), (2, 1), (0, 3)]
        )
        segment = AlignedSegment(mock_segment)
        fasta = segment.to_fasta()
        desired_fasta = np.array(list('AC-TCCTC-GAA'), dtype='<U1')
        self.assertTrue(np.all(fasta==desired_fasta), 'with insertion')

        fasta = segment.to_fasta(False)
        desired_fasta = np.array(list('A-TCCTC-GAA'), dtype='<U1')
        self.assertTrue(np.all(fasta==desired_fasta), 'without insertion')


class TestSeveralMappedReads(unittest.TestCase):

    def setUp(self):
        self.mock_alignment = MockPysamAlignment([
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
        ])
        self.mapped_reads = MappedReads(self.mock_alignment)

    def test_initiate_conversion(self):
        segments = self.mapped_reads.reads
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
        reference_length = 20
        window_start = 2
        window_end = 16
        self.mapped_reads.initiate_conversion(
            reference_length, window_start, window_end
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
            insertion_result = self.mapped_reads.handle_insertions()
            if not insertion_result:
                self.mapped_reads.handle_deletions_and_matches()
            column = self.mapped_reads.fasta[:, j]
            desired_column = [row[j] for row in desired_fasta]
            triplets = zip(range(number_of_rows), column, desired_column)
            for i, character, desired_character in triplets:
                error_message = 'row %d, column %d' % (i, j)
                self.assertEqual(character, desired_character, error_message)


class TestFullAlignment(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sample_alignment = create_sample_alignment()
        cls.mapped_reads = MappedReads(cls.sample_alignment)

    def test_full_single_case_with_insertions(self):
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

        reference_length = 20
        window_start = 2
        window_end = 16
        fasta, reference_position = self.mapped_reads. \
            sam_window_to_fasta_with_insertions(
                reference_length, window_start, window_end
            )
        self.assertTrue(np.all(fasta == desired_fasta))
        self.assertTrue(np.all(
            reference_position == desired_reference_position
        ))

    def test_full_single_case_without_insertions(self):
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

        reference_length = 20
        fasta = self.mapped_reads.sam_window_to_fasta(
            window_start, window_end
        )
        self.assertTrue(np.all(fasta == desired_fasta))

