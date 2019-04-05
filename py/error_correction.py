from math import floor

import numpy as np

from .numeric_sequence import NumericSequence


class ErrorCorrection:
    "Correct errors in a BAM/SAM."

    def __init__(self, pysam_alignment, number_of_reads,
                 neighbors=20, stride=50, slack=5, number_of_processes=1):
        self.pysam_alignment = pysam_alignment
        self.numeric_sequence = NumericSequence()
        self.counts = np.zeros(
            (self.numeric_sequence.n_nucleotides, number_of_reads)
        )
        self.reference_length = pysam_alignment.header['SQ'][0]['LN']
        self.number_of_windows = floor(self.reference_length/stride)
        self.reads_in_window = [[] for _ in range(self.number_of_windows)]
        window_boundaries = np.arange(stride, self.reference_length, stride)
        window_boundaries[-1] = self.reference_length
        self.left_boundaries = window_boundaries - slack
        self.right_boundaries = window_boundaries + slack

    def get_window_range(self, read):
        left_index = np.searchsorted(
            self.left_boundaries,
            read.reference_start
        )
        right_index = np.searchsorted(
            self.right_boundaries,
            read.reference_end
        )
        return left_index, right_index

    def assign_read_to_windows(self, read, read_index):
        left_window_index, right_window_index = self.get_window_range(read)
        for window_index in range(left_window_index, right_window_index+1):
            self.reads_in_window[window_index].append(read_index)

    def assign_all_reads(self):
        print('Assigning reads to respective windows...')
        reads_and_indices = enumerate(self.pysam_alignment.fetch())
        for read_index, mapped_read in reads_and_indices:
            self.assign_read_to_windows(mapped_read, read_index)

