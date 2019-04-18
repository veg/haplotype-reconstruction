import unittest
import os

import numpy as np
import pandas as pd

from .mock import excel_information as excel
from .mock import cigar_string_to_tuples


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


