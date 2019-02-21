import json

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from sklearn.manifold import TSNE
from sklearn.mixture import BayesianGaussianMixture 
from sklearn.cluster import SpectralClustering


def embed_and_reduce_dimensions(records, ndim=2):
    nuc2ind = { '-': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4, 'Y': '2', 'R': 1, 'W': 1 }
    index = []
    rows = []
    for record in records:
        index.append(record.name)
        rows.append([nuc2ind[char] for char in record.seq])
    index = pd.Series(index)
    numeric_fasta = np.array(rows, dtype=np.int)
    embedding = np.array([
      [1, 0, 0, 0, 0],
      [0, 1, 0, 0, 0],
      [0, 0, 1, 0, 0],
      [0, 0, 0, 1, 0],
      [0, 0, 0, 0, 1],
    ])
    
    info_json = {
        'local_starts': [],
        'local_stops': []
    }
    #mean_read_length = np.ceil(np.sum(numeric_fasta != 0, axis=1).mean())
    mean_read_length = np.ceil(np.percentile(np.sum(numeric_fasta != 0, axis=1), 25))
    info_json['mean_read_length'] = mean_read_length
    overlap_fraction = .5
    overlap = np.ceil(overlap_fraction*mean_read_length)
    stride = int(mean_read_length - overlap)
    coverage = np.sum(numeric_fasta != 0, axis=0)
    sufficient_coverage = np.where(coverage >= 20)[0]
    start = sufficient_coverage[0]
    stop = sufficient_coverage[-1]
    local_start = int(start)
    local_stop = int(local_start + mean_read_length)

    reduced_dimension_dfs = []
    i = 0
    while local_stop < stop:
        info_json['local_starts'].append(local_start)
        info_json['local_stops'].append(local_stop)

        local_numeric_fasta = numeric_fasta[:, local_start:local_stop]
        desired_fraction = .8*(local_stop-local_start)
        local_coverage = np.sum(local_numeric_fasta != 0, axis=1)
        local_sufficient_coverage = local_coverage > desired_fraction
        local_numeric_fasta = local_numeric_fasta[local_sufficient_coverage, :]
        n = local_numeric_fasta.shape[0]
        k = local_numeric_fasta.shape[1]
        embedded_fasta = np.reshape(embedding[local_numeric_fasta, :], newshape=(n, 5*k))
        reduced_dimensions = TSNE(
            n_components=ndim, n_iter=5000, random_state=0
        ).fit_transform(embedded_fasta)
        columns = ['d%d' % (j+1) for j in range(ndim)]
        ids = pd.DataFrame({
            'sequence_id': index[local_sufficient_coverage],
            'block': n*[i]
        }).reset_index()
        rd_df = pd.DataFrame(reduced_dimensions, columns=columns).reset_index()
        df = pd.concat([ids, rd_df], axis=1).drop(columns='index')
        reduced_dimension_dfs.append(df)
        
        i += 1
        local_start += stride
        local_stop += stride

    df = pd.concat(reduced_dimension_dfs, axis=0)
    return df, info_json


def embed_and_reduce_dimensions_io(input_fasta, output_csv, output_json, ndim=2):
    records = SeqIO.parse(input_fasta, 'fasta')
    df, info_json = embed_and_reduce_dimensions(records, ndim)
    df.to_csv(output_csv, index_label='id')
    with open(output_json, 'w') as json_file:
        json.dump(info_json, json_file, indent=2)


def cluster_blocks(dr_df, ndim):
    n_blocks = dr_df['block'].max()+1
    dr_df['cluster'] = ''
    columns = ['d%d' % (j+1) for j in range(ndim)]
    for block in range(n_blocks):
        current_block = dr_df.loc[dr_df['block'] == block, columns]
        #clusters = BayesianGaussianMixture(n_components=5).fit(current_block).predict(current_block)
        #dr_df.loc[dr_df['block'] == block, 'cluster'] = ['c%d' % cluster for cluster in clusters]
        clusters = SpectralClustering(2).fit(current_block)
        dr_df.loc[dr_df['block'] == block, 'cluster'] = ['c%d' % cluster for cluster in clusters.labels_]
    return dr_df


def cluster_blocks_io(input_csv, output_csv, ndim):
    dr_df = pd.read_csv(input_csv)
    cluster_blocks(dr_df, ndim).to_csv(output_csv)
    

def obtain_consensus(cluster_df, reads, info_json):
    n_blocks = len(info_json['local_starts'])
    contigs = []
    for block in range(n_blocks):
        current_block = cluster_df.loc[cluster_df['block'] == block, :]
        clusters = list(set(current_block.loc[:, 'cluster']))
        for cluster in clusters:
            reads_in_cluster = current_block.loc[current_block['cluster'] == cluster, 'sequence_id']
            msa = MultipleSeqAlignment([reads[read] for read in reads_in_cluster])
            summary_info = AlignInfo.SummaryInfo(msa)
            seq = str(summary_info.gap_consensus()).replace('X', '-')
            record_id = 'block-%d_cluster-%s' % (block, cluster)
            contig = SeqRecord(Seq(seq), id=record_id, description='')
            contigs.append(contig)
    return contigs


def obtain_consensus_io(input_csv, input_fasta, input_json, output_fasta):
    cluster_df = pd.read_csv(input_csv)
    reads = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
    with open(input_json) as json_file:
        info_json = json.load(json_file)
    contigs = obtain_consensus(cluster_df, reads, info_json)
    SeqIO.write(contigs, output_fasta, 'fasta')

