import os
import re

import pandas as pd

from py import BaseMappedReads
from py import AlignedSegment

excel_information = {
    'read_row_start': 3,
    'read_row_end': 23,
    'read_col_start': 2,
    'read_col_end': 24,
    'query_col': 24,
    'reference_start_col': 25,
    'reference_end_col': 26,
    'cigar_col': 27,
    'read_window_start': 4,
    'read_window_end': 21,
    'reference_position_row': 1,
    'insertion_columns': [8, 16, 17]
}

class MockPysamAlignedSegment:

    def __init__(self, alignment_sequence, cigartuples,
            reference_start=None, reference_end=None):
        self.query_alignment_sequence = alignment_sequence
        self.cigartuples = cigartuples
        self.reference_start = reference_start
        self.reference_end = reference_end
    
    @property
    def pysam_aligned_segment(self):
        return self


class MockMappedReads(BaseMappedReads):

    def __init__(self, mock_segments):
        self.mock_segments = mock_segments

    def fetch(self, start=0, stop=10000000):
        return self.mock_segments


def cigar_string_to_tuples(cigar_string):
    cigar_regex = re.compile('(\d+|[A-Z])');
    tuple_atoms = list(filter(None, cigar_regex.split(cigar_string)))
    lengths = [int(atom) for atom in tuple_atoms[::2]]
    operation_encoding = { 'M': 0, 'I': 1, 'D': 2 }
    operations = [operation_encoding[atom] for atom in tuple_atoms[1::2]]
    cigar_tuples = list(zip(operations, lengths))
    return cigar_tuples


excel_filepath = os.path.join(os.getcwd(), 'tests', 'data', 'SAMtoFASTA.xlsx')

def create_sample_alignment():
    df = pd.read_excel(excel_filepath, header=None)
    queries = df.loc[
        excel_information['read_row_start']:,
        excel_information['query_col']
    ]
    reference_starts = df.loc[
        excel_information['read_row_start']:,
        excel_information['reference_start_col']
    ]
    reference_ends = df.loc[
        excel_information['read_row_start']:,
        excel_information['reference_end_col']
    ]
    cigar_strings = df.loc[
        excel_information['read_row_start']:,
        excel_information['cigar_col']
    ]
    cigar_tuples = [
        cigar_string_to_tuples(cigar_string) for cigar_string in cigar_strings
    ]
    mock_segments = [
        AlignedSegment(MockPysamAlignedSegment(
            query, cigar_tuple, reference_start, reference_end
        ))
        for query, cigar_tuple, reference_start, reference_end
        in zip(queries, cigar_tuples, reference_starts, reference_ends)
    ]
    alignment = MockMappedReads(mock_segments)
    return alignment

