from collections import Counter

import networkx as nx
import pandas as pd
import numpy as np
from Bio import SeqIO
from scipy.stats import rankdata

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
                    covarying_boundaries[0]: covarying_boundaries[1]
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

    def get_aligned_superreads(self, reference, covarying_sites):
        superreads = self.obtain_superreads()

    def create(self):
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
                i_cv_start = superread_i['cv_start']
                i_cv_end = superread_i['cv_end']
                j_cv_start = superread_j['cv_start']
                j_cv_end = superread_j['cv_end']
                start_before_start = i_cv_start <= j_cv_start
                start_before_end = j_cv_start <= i_cv_end
                end_before_end = i_cv_end <= j_cv_end
                if start_before_start and start_before_end and end_before_end:
                    cv_start = max(i_cv_start, j_cv_start)
                    cv_end = min(i_cv_end, j_cv_end)
                    delta = cv_end - cv_start
                    i_start = cv_start - i_cv_start
                    i_end = i_cv_start + delta
                    j_start = cv_start - j_cv_start
                    j_end = j_start + delta
                    i_sequence = superread_i['vacs'][i_start: i_end]
                    j_sequence = superread_j['vacs'][j_start: j_end]
                    if i_sequence == j_sequence:
                        G.add_edge(i, j)
        did_not_connect = [None]
        while len(did_not_connect) > 0:
            did_not_connect = [
                node
                for node in G.nodes
                if len(list(G.pred[node])) == 0 and node != 'source'
            ]
            G.remove_nodes_from(did_not_connect)
        self.superread_graph = nx.algorithms.dag.transitive_reduction(G)
        self.superread_graph = G

    def layout(self):
        pos = nx.spring_layout(self.superread_graph)
        node_data = pd.DataFrame([
            (node_name, position[0], position[1])
            for node_name, position in pos.items()
        ], columns=['name', 'x', 'y_spring']).set_index('name')
        not_source = node_data.index != 'source'
        not_target = node_data.index != 'target'
        actual_nodes = not_source & not_target
        y = rankdata(node_data.loc[actual_nodes, 'y_spring'])
        node_data.loc[actual_nodes, 'y'] = y
        node_data.loc['source', 'y'] = 0
        node_data.loc['target', 'y'] = len(y)
        return node_data

    def candidate_haplotypes(self):
        G = self.superread_graph
        G.nodes['source']['candidate_haplotypes'] = np.array([[]], dtype='<U1')
        G.nodes['source']['cv_end'] = 0
        G.nodes['target']['cv_start'] = len(self.covarying_sites)
        G.nodes['target']['vacs'] = ''
        reverse_post_order = list(nx.dfs_postorder_nodes(G, 'source'))[::-1][1:]
        i = 0
        for descendant in reverse_post_order:
            descendant_start = G.nodes[descendant]['cv_start']
            vacs = G.nodes[descendant]['vacs']
            extended_candidates = []
            for ancestor in G.pred[descendant].keys():
                ancestor_end = G.nodes[ancestor]['cv_end']
                vacs_start_index = ancestor_end - descendant_start
                vacs_np = np.array([list(vacs[vacs_start_index:])], dtype='<U1')
                ancestor_candidates = G.nodes[ancestor]['candidate_haplotypes']
                n_ancestor_candidates = ancestor_candidates.shape[0]
                current_extended_candidates = np.hstack([
                    ancestor_candidates,
                    np.repeat(vacs_np, n_ancestor_candidates, axis=0)
                ])
                extended_candidates.append(current_extended_candidates)
            candidate_haplotypes = np.unique(
                np.vstack(extended_candidates), axis=0
            )
            i += 1
            print('%d of %d' % (i, len(G.nodes)), candidate_haplotypes.shape[0])
            G.nodes[descendant]['candidate_haplotypes'] = candidate_haplotypes
        return G.nodes['target']['candidate_haplotypes']

