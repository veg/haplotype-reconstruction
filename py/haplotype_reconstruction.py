import json

import pysam
import numpy as np
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .error_correction import ErrorCorrection
from .read_graph import SuperReadGraph


def error_correction_io(input_bam, output_bam, output_json):
    alignment = pysam.AlignmentFile(input_bam, 'rb')
    error_correction = ErrorCorrection(alignment)
    error_correction.write_corrected_reads(output_bam)
    with open(output_json, 'w') as json_file:
        json.dump(
            [int(i) for i in error_correction.covarying_sites],
            json_file
        )


def superreads_io(input_reference, input_json, input_bam,
        output_full, output_cvs, output_json):
    with open(input_json) as json_file:
        covarying_sites = np.array(json.load(json_file), dtype=np.int)
    reference = SeqIO.read(input_reference, 'fasta')
    reference_np = np.array(list(str(reference.seq)))
    corrected_reads = pysam.AlignmentFile(input_bam, 'rb')

    superread_graph = SuperReadGraph(corrected_reads, covarying_sites)
    superread_info = superread_graph.obtain_superreads()
    superread_graph.create()
    layout = superread_graph.layout()

    superread_records = []
    for i, info in enumerate(superread_info):
        current_sequence = np.array(len(reference.seq) * ['-'], dtype='<U1')
        start = covarying_sites[info['cv_start']]
        end = covarying_sites[info['cv_end']-1]
        current_sequence[start: end] = reference_np[start: end]
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

    superread_json = {
        'info': superread_info,
        'layout': layout.to_dict('records'),
        'links': nx.node_link_data(superread_graph.superread_graph)
    }
    with open(output_json, 'w') as json_file:
        json.dump(superread_json, json_file, indent=2)

