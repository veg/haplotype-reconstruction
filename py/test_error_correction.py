import unittest
import os

import pysam

from .error_correction import ErrorCorrection


class MockPysamAlignedSegment:

    def __init__(self, reference_start, reference_end):
        self.reference_start = reference_start
        self.reference_end = reference_end


class MockPysamAlignment:

    def __init__(self, reference_length=100):
        self.header = { 'SQ': [ {'LN': reference_length } ] }
        self.mapped_reads = [
            MockPysamAlignedSegment(0, 31)
        ]
        self.number_of_reads = len(self.mapped_reads)

    def fetch(self):
        return self.mapped_reads


class TestErrorCorrection(unittest.TestCase):

    def setUp(self):
        self.mock_alignment = MockPysamAlignment()
        self.mock_error_correction = ErrorCorrection(
            self.mock_alignment, self.mock_alignment.number_of_reads,
            stride=30, slack=2
        )
        real_alignment_path = os.path.join(os.getcwd(), 'py', 'test-data', 'sorted.bam')

        self.real_alignment = pysam.AlignmentFile(real_alignment_path, 'rb')
        # obtained with: samtools idxstats py/test-data/sorted.bam
        real_number_of_reads = 42746 
        self.real_error_correction = ErrorCorrection(
            self.real_alignment, real_number_of_reads
        )

    def test_get_window_range(self):
        first_read = self.mock_alignment.mapped_reads[0]
        left_index, right_index = self.mock_error_correction.get_window_range(first_read)
        self.assertEqual(left_index, 0)
        self.assertEqual(right_index, 0)

    def test_assign_read_to_windows(self):
        first_read = self.mock_alignment.mapped_reads[0]
        self.mock_error_correction.assign_read_to_windows(first_read, 0)
        desired = [ [0], [], [] ]
        self.assertEqual(desired, self.mock_error_correction.reads_in_window)

    def test_assign_all_reads(self):
        self.real_error_correction.assign_all_reads()

