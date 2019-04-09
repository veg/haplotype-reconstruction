from math import floor
from multiprocessing import Pool

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.neighbors import KDTree

from .numeric_sequence import NumericSequence


def cluster_and_correct(packed_arguments):
    index, numeric_window, window_start = packed_arguments
    print('Starting on window %d...' % index)
    print('...finished window %d!' % index)
    return pd.DataFrame({
        'read_index': pd.Series([], dtype='int'),
        'position_index': pd.Series([], dtype='int'),
        'correction': pd.Series([], dtype='<U1')
    })


class ErrorCorrection:
    "Correct errors in a BAM/SAM."

    def __init__(self, pysam_alignment, number_of_reads,
                 neighbors=20, stride=50, slack=5, number_of_processes=1):
        self.pysam_alignment = pysam_alignment
        self.numeric_sequence = NumericSequence()
        self.reference_length = pysam_alignment.header['SQ'][0]['LN']
        self.number_of_windows = floor(self.reference_length/stride)
        self.reads_in_window = [[] for _ in range(self.number_of_windows)]
        self.window_boundaries = np.arange(0, self.reference_length, stride)
        self.window_boundaries[-1] = self.reference_length
        self.left_boundaries = self.window_boundaries[:-1] + slack
        self.right_boundaries = self.window_boundaries[1:] - slack
        self.number_of_processes = number_of_processes
        self.corrections = None

    def get_window_range(self, read):
        left_index = np.searchsorted(
            self.left_boundaries,
            read.reference_start
        )
        right_index = np.searchsorted(
            self.right_boundaries,
            read.reference_end
        ) - 1
        return left_index, right_index

    def assign_read_to_windows(self, read, read_index):
        left_window_index, right_window_index = self.get_window_range(read)
        for window_index in range(left_window_index, right_window_index + 1):
            self.reads_in_window[window_index].append(read_index)

    def assign_all_reads(self):
        print('Assigning reads to respective windows...')
        reads_and_indices = enumerate(self.pysam_alignment.fetch())
        for read_index, mapped_read in reads_and_indices:
            self.assign_read_to_windows(mapped_read, read_index)

    def get_numeric_window_block(self, lower_window, upper_window):
        numeric_windows = []
        for window_index in range(lower_window, upper_window):
            number_of_rows = len(self.reads_in_window[window_index])
            current_start = self.window_boundaries[window_index]
            current_stop = self.window_boundaries[window_index + 1]
            number_of_columns = current_stop - current_start
            shape = (number_of_rows, number_of_columns)
            numeric_windows.append(np.zeros(shape))
        current_indices = self.number_of_windows * [0]
        for read_index, read in enumerate(self.pysam_alignment.fetch()):
            numeric_row = self.numeric_sequence.fasta_from_sam(read)
            for window_index in range(lower_window, upper_window):
                if read_index in self.reads_in_window[window_index]:
                    current_indices[window_index] += 1
        return numeric_windows

    def get_corrections(self):
        lower_window = 0
        all_correction_dfs = [[] for _ in range(self.number_of_windows)]
        print('Processing %d blocks...' % self.number_of_windows) 
        pool = Pool(processes=self.number_of_processes)
        while lower_window < self.number_of_windows:
            print('executing while')
            upper_window = min(
                lower_window + self.number_of_processes,
                self.number_of_windows
            )
            numeric_windows = self.get_numeric_window_block(
                lower_window, upper_window
            )
            arguments = zip(
                range(lower_window, upper_window),
                numeric_windows,
                self.window_boundaries[lower_window: upper_window]
            )
            corrections = pool.map(cluster_and_correct, arguments)
            print('made it past map')
            all_correction_dfs[lower_window: upper_window] = corrections
            lower_window += self.number_of_processes
        print('Done all blocks!')
        pool.close()
        self.corrections = pd.concat(all_correction_dfs, ignore_index=True)
        self.corrections.sort_values(
            by=['read_index', 'position_index'], inplace=True
        )

    def apply_corrections(self):
        pass

    def perform_error_correction(self):
        self.assign_all_reads()
        self.get_corrections()
        self.apply_corrections()
        return self.pysam_alignment, self.corrections

