import unittest
import json

import numpy as np
import pysam

from py import SuperReadGraph


class TestSuperReadGraph(unittest.TestCase):
    bam_path = 'tests/data/corrected.bam'
    json_path = 'tests/data/covarying_sites.json'

    def test_obtain_superreads(self):
        alignment = pysam.AlignmentFile(self.bam_path)
        with open(self.json_path) as json_file:
            covarying_sites = np.array(json.load(json_file))
        superread_graph = SuperReadGraph(alignment, covarying_sites)
        superreads = superread_graph.obtain_superreads()
        for superread in superreads:
            print(superread)
        print('Total superreads: ', len(superreads))

    def test_create_superread_graph(self):
        alignment = pysam.AlignmentFile(self.bam_path)
        with open(self.json_path) as json_file:
            covarying_sites = np.array(json.load(json_file))
        superread_graph = SuperReadGraph(alignment, covarying_sites)
        superread_graph.create()

    def test_candidate_haplotypes(self):
        alignment = pysam.AlignmentFile(self.bam_path)
        with open(self.json_path) as json_file:
            covarying_sites = np.array(json.load(json_file))
        superread_graph = SuperReadGraph(alignment, covarying_sites)
        superread_graph.create()
        superread_graph.candidate_haplotypes()
