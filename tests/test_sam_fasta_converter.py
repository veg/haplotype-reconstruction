import unittest
import re
import os

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from py import SAMFASTAConverter
from py import AlignedSegment


class MockPysamAlignedSegment:
    def __init__(self, alignment_sequence, cigartuples, reference_start=None):
        self.query_alignment_sequence = alignment_sequence
        self.cigartuples = cigartuples
        self.reference_start = reference_start


def cigar_string_to_tuples(cigar_string):
    cigar_regex = re.compile('(\d+|[A-Z])');
    tuple_atoms = list(filter(None, cigar_regex.split(cigar_string)))
    lengths = [int(atom) for atom in tuple_atoms[::2]]
    operation_encoding = { 'M': 0, 'I': 1, 'D': 2 }
    operations = [operation_encoding[atom] for atom in tuple_atoms[1::2]]
    cigar_tuples = list(zip(operations, lengths))
    return cigar_tuples


class TestCigarStringToTuples(unittest.TestCase):
    def test_several_cigar_strings(self):
        cigar_strings = [
            '13M',
            '9M1D3M',
            '1M1D11M',
            '14M',
            '4M1I6M1D3M',
            '10M2I3M',
            '9M1I5M',
            '8M1D3M',
            '3M1D9M'
        ]
        desired_cigar_tuples = [
            [(0, 13)],
            [(0, 9), (2, 1), (0, 3)],
            [(0, 1), (2, 1), (0, 11)],
            [(0, 14)],
            [(0, 4), (1, 1), (0, 6), (2, 1), (0, 3)],
            [(0, 10), (1, 2), (0, 3)],
            [(0, 9), (1, 1), (0, 5)],
            [(0, 8), (2, 1), (0, 3)],
            [(0, 3), (2, 1), (0, 9)]
        ]
        cigar_pairs = zip(cigar_strings, desired_cigar_tuples)
        for cigar_string, desired_cigar_tuple in cigar_pairs:
            cigar_tuple = cigar_string_to_tuples(cigar_string)
            self.assertEqual(cigar_tuple, desired_cigar_tuple)


excel = {
    'read_row_start': 3,
    'read_row_end': 23,
    'read_col_start': 2,
    'read_col_end': 24,
    'query_col': 24,
    'reference_start_col': 25,
    'cigar_col': 27,
    'read_window_start': 4,
    'read_window_end': 21,
    'reference_position_row': 1
}


class TestValidSampleExcel(unittest.TestCase):
    def setUp(self):
        excel_filepath = os.path.join(os.getcwd(), 'tests', 'data', 'SAMtoFASTA.xlsx')
        self.df = pd.read_excel(excel_filepath, header=None)
        self.reads = self.df.iloc[
            excel['read_row_start']: excel['read_row_end'],
            excel['read_col_start']: excel['read_col_end']
        ]
        self.reads[self.reads == '~'] = np.nan
    
    def test_queries_match_alignment(self):
        queries = self.df.iloc[
            excel['read_row_start']: excel['read_row_end'],
            excel['query_col']
        ]
        for i, read in self.reads.iterrows():
            not_gap = read != '-'
            valid_sites = not_gap & read.notna()
            fasta_query = ''.join(read[valid_sites])
            query = queries[i]
            self.assertEqual(fasta_query, query, i)

    def test_reference_start_agreement(self):
        offset = 2
        reference_starts = self.df.iloc[
            excel['read_row_start']: excel['read_row_end'],
            excel['reference_start_col']
        ]
        for i, read in self.reads.iterrows():
            reference_start = reference_starts[i]
            fasta_reference_start = int(read.first_valid_index()) - offset
            if i >= 20: # apply correction for first gap
                fasta_reference_start -= 1
            self.assertEqual(reference_start, fasta_reference_start, i)


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
        desired_fasta = np.array(list('AC-TCCTC-GAA'), dtype='<U1')
        self.assertTrue(np.all(fasta==desired_fasta))


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

    def test_partial_single_case(self):
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

    def test_full_single_case(self):
        excel_filepath = os.path.join(os.getcwd(), 'tests', 'data', 'SAMtoFASTA.xlsx')
        df = pd.read_excel(excel_filepath, header=None)
        desired_fasta = np.array(
            df.iloc[
                excel['read_row_start']: excel['read_row_end'],
                excel['read_window_start'] : excel['read_window_end']
            ],
            dtype='<U1'
        )
        queries = df.loc[excel['read_row_start']:, excel['query_col']]
        reference_starts = df.loc[
            excel['read_row_start']:, excel['reference_start_col']
        ]
        cigar_strings = df.loc[excel['read_row_start']:, excel['cigar_col']]
        cigar_tuples = [
            cigar_string_to_tuples(cigar_string) for cigar_string in cigar_strings
        ]
        triplets = zip(queries, cigar_tuples, reference_starts)
        mock_segments = [
            MockPysamAlignedSegment(query, cigar_tuple, reference_start)
            for query, cigar_tuple, reference_start in triplets
        ]
        desired_reference_position = np.array(df.iloc[
            excel['reference_position_row'],
            excel['read_window_start']: excel['read_window_end']
        ], dtype=np.int)

        sam_fasta_converter = SAMFASTAConverter()
        reference_length = 20
        window_start = 2
        window_end = 16
        fasta, reference_position = sam_fasta_converter.sam_window_to_fasta(
            mock_segments, reference_length, window_start, window_end
        )
        self.assertTrue(np.all(fasta == desired_fasta))
        self.assertTrue(np.all(
            reference_position == desired_reference_position
        ))

