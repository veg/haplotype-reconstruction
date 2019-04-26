from collections import Counter

import networkx as nx
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
                    'cv_start': int(covarying_boundaries[0]),
                    'cv_end': int(covarying_boundaries[1]),
                })
        self.superreads = all_superreads
        return all_superreads

    def create_superread_graph(self):
        superreads = self.obtain_superreads()
        G = nx.DiGraph()
        G.add_node('source')
        G.add_node('target')
        n_cv = len(self.covarying_sites)
        for i, superread in enumerate(superreads):
            G.add_node(i, **superread)
            if superread['cv_start'] == 0:
                G.add_edge('source', i)
            if superread['cv_end'] == n_cv:
                G.add_edge(i, 'target')
        for i, superread_i in enumerate(superreads):
            for j, superread_j in enumerate(superreads):
                if i == j:
                    continue
                i_start = superread_i['cv_start']
                i_end = superread_i['cv_end']
                j_start = superread_j['cv_start']
                j_end = superread_j['cv_end']
                start_before_start = i_start <= j_start
                start_before_end = j_start <= i_end
                end_before_end = i_end <= j_end
                if start_before_start and start_before_end and end_before_end:
                    cv_start = max(i_start, j_start)
                    cv_end = min(i_end, j_end)
                    delta = cv_end - cv_start
                    i_start = cv_start - i_start
                    i_end = i_start + delta
                    j_start = cv_start - j_start
                    j_end = j_start + delta
                    i_sequence = superread_i['vacs'][i_start:i_end]
                    j_sequence = superread_j['vacs'][j_start:j_end]
                    if i_sequence == j_sequence:
                        G.add_edge(i, j)
        did_not_connect = [
            node 
            for node in G.nodes
            if len(list(G.pred[node])) == 0 and node != 'source'
        ]
        G.remove_nodes_from(did_not_connect)
        self.superread_graph = nx.algorithms.dag.transitive_reduction(G)
        self.superread_graph = G

    def candidate_haplotypes(self):
        pass
    

