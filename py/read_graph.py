import itertools as it

import networkx as nx
import pandas as pd
import numpy as np
from Bio import SeqIO
from scipy.stats import rankdata

from .mapped_reads import MappedReads
from .utils import read_reference_start_and_end


class SuperReadGraph:
    def __init__(self, mapped_reads=None, covarying_sites=None, node_link_data=None):
        self.mapped_reads = mapped_reads
        self.covarying_sites = covarying_sites
        if node_link_data:
            self.superread_graph = nx.readwrite.json_graph.node_link_graph(
                node_link_data, multigraph=False
            )
        else:
            self.superread_graph = None

    def obtain_superreads(self, minimum_weight=3):
        read_information = read_reference_start_and_end(
            self.mapped_reads, self.covarying_sites
        )
        read_groups = {}
        for i, read in enumerate(self.mapped_reads.fetch()):
            row = read_information.iloc[i, :]
            key = (row['covarying_start'], row['covarying_end'])
            if key in read_groups:
                read_groups[key].append(read)
            else:
                read_groups[key] = [read]
        all_superreads = []
        superread_index = 0
        for covarying_boundaries, read_group in read_groups.items():
            superreads = {}
            for read in read_group:
                fasta = read.to_fasta()
                covarying_sites_in_read = self.covarying_sites[
                    covarying_boundaries[0]: covarying_boundaries[1]
                ]
                corrected_sites = covarying_sites_in_read - read.reference_start
                value_at_covarying_sites = ''.join(fasta[corrected_sites])
                if value_at_covarying_sites in superreads:
                    superreads[value_at_covarying_sites] += 1
                else:
                    superreads[value_at_covarying_sites] = 1
            admission = lambda pair: pair[1] >= minimum_weight
            admissible_superreads = list(filter(admission, superreads.items()))
            total_weight = sum([
                superread[1] for superread in admissible_superreads
            ])
            for vacs, weight in admissible_superreads:
                if weight >= minimum_weight:
                    all_superreads.append({
                        'index': superread_index,
                        'vacs': vacs,
                        'weight': weight,
                        'frequency': weight/total_weight,
                        'cv_start': int(covarying_boundaries[0]),
                        'cv_end': int(covarying_boundaries[1]),
                        'discarded': False
                    })
                    superread_index += 1
        self.superreads = all_superreads
        return all_superreads

    def get_aligned_superreads(self, reference, covarying_sites):
        superreads = self.obtain_superreads()

    def check_compatability(self, superread_i, superread_j, minimum_overlap=2):
        if superread_i['index'] == superread_j['index']:
            return (False, 0)
        i_cv_start = superread_i['cv_start']
        i_cv_end = superread_i['cv_end']
        j_cv_start = superread_j['cv_start']
        j_cv_end = superread_j['cv_end']
        start_before_start = i_cv_start <= j_cv_start
        start_before_end = j_cv_start < i_cv_end
        end_before_end = i_cv_end <= j_cv_end
        if start_before_start and start_before_end and end_before_end:
            cv_start = max(i_cv_start, j_cv_start)
            cv_end = min(i_cv_end, j_cv_end)
            delta = cv_end - cv_start
            i_start = cv_start - i_cv_start
            i_end = i_start + delta
            j_start = cv_start - j_cv_start
            j_end = j_start + delta
            i_sequence = superread_i['vacs'][i_start: i_end]
            j_sequence = superread_j['vacs'][j_start: j_end]
            agree_on_overlap = i_sequence == j_sequence
            overlap = len(i_sequence)
            long_enough = overlap >= minimum_overlap
            compatible = agree_on_overlap and long_enough
            return (compatible, overlap)
        return (False, 0)

    def create_full(self, **kwargs):
        superreads = self.obtain_superreads(**kwargs)
        G = nx.DiGraph()
        G.add_node('source')
        G.add_node('target')
        n_cv = len(self.covarying_sites)
        for superread in superreads:
            G.add_node(superread['index'], **superread)
            if superread['cv_start'] == 0:
                G.add_edge('source', superread['index'])
            if superread['cv_end'] == n_cv:
                G.add_edge(superread['index'], 'target')
        for superread_i in superreads:
            for superread_j in superreads:
                should_include_edge, overlap = self.check_compatability(
                    superread_i, superread_j
                )
                if should_include_edge:
                    G.add_edge(
                        superread_i['index'], superread_j['index'],
                        overlap=overlap
                    )
        self.superread_graph = G

    def node_pruner(self, node):
        G = self.superread_graph
        starts_late = len(list(G.pred[node])) == 0 and node != 'source'
        ends_early = len(list(G.succ[node])) == 0 and node != 'target'
        return starts_late or ends_early

    def reduce(self):
        G = self.superread_graph

        for node in G.nodes:
            if node == 'source' or node == 'target':
                continue
            successors = [succ for succ in G.succ[node] if succ != 'target']
            if len(successors) < 2:
                continue
            edges = [G.edges[node, succ] for succ in successors]
            overlaps = [edge['overlap'] for edge in edges]
            #min_overlap = min(overlaps)
            #min_index = overlaps.index(min_overlap)
            #G.remove_edge(node, successors[min_index])
            max_overlap = max(overlaps)
            max_index = overlaps.index(max_overlap)
            for i, succ in enumerate(successors):
                if i != max_index:
                    G.remove_edge(node, succ)

        did_not_connect = [None]
        while len(did_not_connect) > 0:
            did_not_connect = [
                node for node in G.nodes if self.node_pruner(node)
            ]
            print('Removing ', ' '.join([str(i) for i in did_not_connect]))
            G.remove_nodes_from(did_not_connect)
            for node in did_not_connect:
                try:
                    self.superreads[node]['discarded'] = True
                except TypeError:
                    import pdb; pdb.set_trace()
        self.superread_graph = nx.algorithms.dag.transitive_reduction(G)
        for node in G.nodes:
            self.superread_graph.nodes[node].update(G.nodes[node])
        for edge in self.superread_graph.edges:
            source, target = edge
            self.superread_graph.edges[source, target].update(G.edges[source, target])

    def create(self, **kwargs):
        self.create_full(**kwargs)
        self.reduce()

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
        number_of_covarying_sites = max(
            [G.nodes[node]['cv_end'] for node in G.nodes if 'cv_end' in G.nodes[node]]
        )
        G.nodes['source']['candidate_haplotypes'] = np.array([[]], dtype='<U1')
        G.nodes['source']['describing_superreads'] = [[]]
        G.nodes['source']['cv_end'] = 0
        G.nodes['target']['cv_start'] = number_of_covarying_sites
        G.nodes['target']['cv_end'] = number_of_covarying_sites
        G.nodes['target']['vacs'] = ''
        reverse_post_order = list(nx.dfs_postorder_nodes(G, 'source'))[::-1][1:]

        G.nodes['source']['npath'] = 1
        for node in reverse_post_order:
            number_of_current_paths = sum([
                G.nodes[predecessor]['npath']
                for predecessor in G.predecessors(node)
            ])
            G.nodes[node]['npath'] = number_of_current_paths
        print('Obtaining', G.nodes['target']['npath'], 'candidate haplotypes...')

        for descendant in reverse_post_order:
            descendant_start = G.nodes[descendant]['cv_start']
            vacs = G.nodes[descendant]['vacs']
            extended_candidates = []
            extended_describing_superreads = []
            for ancestor in G.pred[descendant].keys():
                ancestor_end = G.nodes[ancestor]['cv_end']
                vacs_start_index = ancestor_end - descendant_start
                vacs_np = np.array([list(vacs[vacs_start_index:])], dtype='<U1')
                ancestor_candidates = G.nodes[ancestor]['candidate_haplotypes']
                ancestor_superreads = G.nodes[ancestor]['describing_superreads']
                n_ancestor_candidates = ancestor_candidates.shape[0]
                current_extended_candidates = np.hstack([
                    ancestor_candidates,
                    np.repeat(vacs_np, n_ancestor_candidates, axis=0)
                ])
                superread_extension = [descendant] if descendant != 'target' else []
                current_extended_describing_superreads = [
                    describing_list + superread_extension
                    for describing_list in ancestor_superreads
                ]
                extended_candidates.append(current_extended_candidates)
                extended_describing_superreads.append(
                    current_extended_describing_superreads
                )
            candidate_haplotypes = np.vstack(extended_candidates)
            G.nodes[descendant]['candidate_haplotypes'] = candidate_haplotypes
            G.nodes[descendant]['describing_superreads'] = list(
                it.chain.from_iterable(extended_describing_superreads)
            )
        candidate_haplotypes = G.nodes['target']['candidate_haplotypes']
        describing_superreads = G.nodes['target']['describing_superreads']
        return candidate_haplotypes, describing_superreads

