import unittest
import os

import pysam

from py import ErrorCorrection


class MockPysamAlignedSegment:

    def __init__(self, reference_start, reference_end):
        self.reference_start = reference_start
        self.reference_end = reference_end


class MockPysamAlignment:

    def __init__(self, reference_length=100):
        self.header = { 'SQ': [ {'LN': reference_length } ] }

        maps_snugly_in_first = MockPysamAlignedSegment(0, 30)
        underlaps_first_slightly = MockPysamAlignedSegment(0, 29)
        too_small_for_first = MockPysamAlignedSegment(2, 25)
        overlaps_first_and_second = MockPysamAlignedSegment(0, 61)
        overlaps_second_and_third = MockPysamAlignedSegment(22, 100)

        self.mapped_reads = [
            maps_snugly_in_first,
            underlaps_first_slightly,
            too_small_for_first,
            overlaps_first_and_second,
            overlaps_second_and_third 
        ]
        self.desired_window_ranges = [
            {'left': 0, 'right': 0},
            {'left': 0, 'right': 0},
            {'left': 0, 'right': -1},
            {'left': 0, 'right': 1},
            {'left': 1, 'right': 2}
        ]
        self.desired_reads_in_window = [
            [[0], [], []],
            [[0, 1], [], []],
            [[0, 1], [], []],
            [[0, 1, 3], [3], []],
            [[0, 1, 3], [3, 4], [4]]
        ]
        assert len(self.mapped_reads) == len(self.desired_window_ranges)
        assert len(self.mapped_reads) == len(self.desired_reads_in_window)
        self.number_of_reads = len(self.mapped_reads)

    def fetch(self):
        return self.mapped_reads


class TestArtificialErrorCorrection(unittest.TestCase):

    def setUp(self):
        self.mock_alignment = MockPysamAlignment()
        self.mock_error_correction = ErrorCorrection(
            self.mock_alignment, self.mock_alignment.number_of_reads,
            stride=30, slack=2
        )

    def test_get_window_range(self):
        window_range_data = zip(
            range(self.mock_alignment.number_of_reads),
            self.mock_alignment.mapped_reads,
            self.mock_alignment.desired_window_ranges
        )
        for index, read, desired in window_range_data:
            left_index, right_index = self.mock_error_correction.get_window_range(read)
            self.assertEqual(left_index, desired['left'], index)
            self.assertEqual(right_index, desired['right'], index)

    def test_assign_read_to_windows(self):
        assign_read_data = zip(
            range(self.mock_alignment.number_of_reads),
            self.mock_alignment.mapped_reads,
            self.mock_alignment.desired_reads_in_window
        )
        for index, read, desired_assignment in assign_read_data:
            self.mock_error_correction.assign_read_to_windows(read, index)
            self.assertEqual(
                desired_assignment,
                self.mock_error_correction.reads_in_window
            )


class TestRealErrorCorrection(unittest.TestCase):

    def setUp(self):
        alignment_path = os.path.join(os.getcwd(), 'tests', 'data', 'sorted.bam')
        self.alignment = pysam.AlignmentFile(alignment_path, 'rb')
        # obtained with: samtools idxstats py/test-data/sorted.bam
        number_of_reads = 42746 
        self.error_correction = ErrorCorrection(
            self.alignment, number_of_reads,
            number_of_processes=10
        )
    def test_full_pipeline(self):
        self.error_correction.assign_all_reads()
        self.error_correction.get_corrections()
