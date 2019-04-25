from collections import Counter

import networkx
import pandas as pd
import numpy as np

from .sam_fasta_converter import SAMFASTAConverter


class SuperReadGraph:
    def __init__(self, pysam_alignment=None, covarying_sites=None):
        self.pysam_alignment = pysam_alignment
        self.covarying_sites = covarying_sites

        self.superreads = []
        self.superread_graph = None

    def obtain_superreads(self):
        read_information = pd.DataFrame(
            [
                (read.reference_start, read.reference_end)
                for read in self.pysam_alignment.fetch()
            ],
            columns=['reference_start', 'reference_end']
        )
        read_information['covarying_start'] = np.searchsorted(
            self.covarying_sites, read_information['reference_start']
        )
        read_information['covarying_end'] = np.searchsorted(
            self.covarying_sites, read_information['reference_end']
        )
        read_groups = {}
        for i, read in enumerate(self.pysam_alignment.fetch()):
            row = read_information.iloc[i, :]
            key = (row['covarying_start'], row['covarying_end'])
            if key in read_groups:
                read_groups[key].append(read)
            else:
                read_groups[key] = [read]
        all_superreads = []
        sfc = SAMFASTAConverter()
        for covarying_boundaries, read_group in read_groups.items():
            superreads = {}
            for read in read_group:
                fasta = sfc.single_segment_to_fasta(read)
                covarying_sites_in_read = self.covarying_sites[
                    covarying_boundaries[0]:covarying_boundaries[1]
                ]
                corrected_sites = covarying_sites_in_read - read.reference_start
                value_at_covarying_sites = ''.join(fasta[corrected_sites])
                if value_at_covarying_sites in superreads:
                    superreads[value_at_covarying_sites] += 1
                else:
                    superreads[value_at_covarying_sites] = 1
            for vacs, weight in superreads.items():
                all_superreads.append({
                    'vacs': vacs,
                    'weight': weight,
                    'cv_start': covarying_boundaries[0],
                    'cv_end': covarying_boundaries[1],
                })
        return all_superreads

    def create_superread_graph(self):
        pass

    def transitive_reduction(self):
        pass

    def candidate_haplotypes(self):
        pass
    

