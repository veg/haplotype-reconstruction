import json

import pysam
import numpy as np
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .mapped_reads import MappedReads
from .error_correction import ErrorCorrection
from .read_graph import SuperReadGraph
from .regression import perform_regression


def error_correction_io(input_bam, output_bam, output_json, output_consensus):
    alignment = pysam.AlignmentFile(input_bam, 'rb')
    error_correction = ErrorCorrection(alignment)
    error_correction.write_corrected_reads(output_bam)
    with open(output_json, 'w') as json_file:
        json.dump(
            [int(i) for i in error_correction.covarying_sites],
            json_file
        )
    consensus_record = error_correction.consensus()
    SeqIO.write(consensus_record, output_consensus, 'fasta')


def superreads_io(input_consensus, input_json, input_bam,
        output_full, output_cvs, output_json):
    with open(input_json) as json_file:
        covarying_sites = np.array(json.load(json_file), dtype=np.int)
    consensus = SeqIO.read(input_consensus, 'fasta')
    consensus_np = np.array(list(str(consensus.seq)))
    corrected_reads = MappedReads(input_bam, 'rb')

    superread_graph = SuperReadGraph(corrected_reads, covarying_sites)
    superread_info = superread_graph.obtain_superreads()
    superread_graph.create_full()
    layout = superread_graph.layout()

    superread_records = []
    for i, info in enumerate(superread_info):
        current_sequence = np.array(len(consensus.seq) * ['-'], dtype='<U1')
        start = covarying_sites[info['cv_start']]
        end = covarying_sites[info['cv_end']-1]
        current_sequence[start: end] = consensus_np[start: end]
        covarying_indices = covarying_sites[info['cv_start']: info['cv_end']]
        current_sequence[covarying_indices] = np.array(
            list(info['vacs']), dtype='<U1'
        )
        superread_records.append(SeqRecord(
            Seq(''.join(current_sequence)),
            id='superread-%d_weight-%d' % (i, info['weight']),
            description=''
        ))
    SeqIO.write(superread_records, output_full, 'fasta')

    superread_records = []
    G = superread_graph.superread_graph
    for node in G.nodes:
        if node != 'source' and node != 'target':
            info = G.nodes[node]
            current_sequence = np.array(len(covarying_sites) * ['-'], dtype='<U1')
            start = info['cv_start']
            end = info['cv_end']
            current_sequence[start: end] = np.array(
                list(info['vacs']), dtype='<U1'
            )
            superread_records.append(SeqRecord(
                Seq(''.join(current_sequence)),
                id='superread-%d_weight-%d' % (node, info['weight']),
                description=''
            ))
    SeqIO.write(superread_records, output_cvs, 'fasta')

    superread_graph.reduce()
    superread_json = {
        'layout': layout.to_dict('records'),
        'nodeLinkData': nx.node_link_data(superread_graph.superread_graph)
    }
    with open(output_json, 'w') as json_file:
        json.dump(superread_json, json_file, indent=2)


def candidates_io(input_consensus, input_graph_json, input_cvs_json,
        output_fasta, output_json):
    with open(input_graph_json) as json_file:
        superread_json = json.load(json_file)
    with open(input_cvs_json) as json_file:
        covarying_sites = np.array(json.load(json_file), dtype=np.int)
    consensus = SeqIO.read(input_consensus, 'fasta')
    consensus_np = np.array(list(str(consensus.seq)))
    superread_graph = SuperReadGraph(
        node_link_data=superread_json['nodeLinkData']
    )

    candidate_haplotypes, describing_superreads = superread_graph.candidate_haplotypes()
    superread_records = []
    for i in range(candidate_haplotypes.shape[0]):
        current_sequence = np.copy(consensus_np)
        current_sequence[covarying_sites] = candidate_haplotypes[i, :]
        superread_records.append(SeqRecord(
            Seq(''.join(current_sequence)),
            id='candidate_%d' % i,
            description=''
        ))
    SeqIO.write(superread_records, output_fasta, 'fasta')
    
    with open(output_json, 'w') as json_file:
        json.dump(describing_superreads, json_file, indent=2)


def regression_io(input_superreads, input_candidates, output_fasta):
    with open(input_superreads) as superread_file:
        superreads = json.load(superread_file)['nodeLinkData']['nodes']

    with open(input_candidates) as candidates_file:
        candidates = json.load(candidates_file)

    perform_regression(superreads, candidates)

    with open(output_fasta, 'w') as fasta_file:
        fasta_file.write('to do')

