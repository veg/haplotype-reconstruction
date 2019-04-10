from math import floor
from multiprocessing import Pool

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.neighbors import KDTree

from .numeric_sequence import NumericSequence


def cluster_and_correct(packed_arguments):
    print('\tStarting on window %d...' % index)

    pysam_alignment, read_indices_from_window, boundaries, index = packed_arguments
    read_indices_from_window = np.array(read_indices_from_window, dtype=np.int)
    left_boundary, right_boundary = boundaries
    numeric_sequence = NumericSequence()
    nrows = len(read_indices_from_window)
    ncols = 5*(right_boundary - left_boundary)
    embedding = np.zeros((nrows, ncols))

    current_row = 0
    for i, read in enumerate(pysam_alignment.fetch()):
        read_index = np.searchsorted(read_indices_from_window, i)
        if read_indices_from_window[read_index] == i:
            numeric_sequence.fasta_from_sam(read)
            

    print('\t...finished window %d!' % index)
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

    def get_corrections(self):
        print('Processing %d windows...' % self.number_of_windows) 
        pool = Pool(processes=self.number_of_processes)
        arguments = zip(
            n*[self.pysam_alignment],
            self.reads_in_window
        )
        all_correction_dfs = pool.map(cluster_and_correct, arguments)
        print('...done all windows!')
        pool.close()
        self.corrections = pd.concat(all_correction_dfs, ignore_index=True)
        self.corrections.sort_values(
            by=['read_index', 'position_index'], inplace=True
        )

    def apply_corrections(self):
        pass

    def perform_error_correction(self):
        print('------------------------------------------------')
        print('--- ACME Viral Quasispecies Error Correction ---')
        print('------------------------------------------------\n')
        self.assign_all_reads()
        self.get_corrections()
        self.apply_corrections()
        return self.pysam_alignment, self.corrections


def error_correction_io(input_bam, output_bam, output_csv):
    pass

