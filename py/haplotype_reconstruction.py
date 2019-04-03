import json
from multiprocessing import Pool
from collections import OrderedDict
from collections import Counter

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
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import NearestNeighbors


def create_numeric_fasta(records):
    iupac_nucleotides = list('ACGTRYSWKMBDHVN-')
    nuc2ind = { nuc: i for i, nuc in enumerate(iupac_nucleotides) }
    index = []
    rows = []
    for record in records:
        index.append(record.name)
        rows.append([nuc2ind[char] for char in record.seq])
    index = pd.Series(index)
    numeric_fasta = np.array(rows, dtype=np.int)
    return index, numeric_fasta


def embed_and_reduce_dimensions(arguments):
    gap_index = 15
    numeric_fasta, local_start, local_stop, index, i = arguments
    embedding = np.array([
      [1,   0,   0,    0,   0], # A
      [0,   1,   0,    0,   0], # C
      [0,   0,   1,    0,   0], # G
      [0,   0,   0,    1,   0], # T
      [.5,  0,  .5,    0,   0], # R
      [0,   .5,  0,    .5,  0], # Y
      [0,   .5,  .5,   0,   0], # S
      [.5,  0,   0,    .5,  0], # W
      [0,   0,   .5,   .5,  0], # K
      [.5,  .5,  0,    0,   0], # M
      [0,   1/3, 1/3,  1/3, 0], # B
      [1/3, 0,   1/3,  1/3, 0], # D
      [1/3, 1/3, 0,    0/3, 0], # H
      [1/3, 1/3, 1/3,  0,   0], # V
      [.25, .25, .25, .25,  0], # N
      [0,   0,   0,    0,   1], # -
    ])
    local_numeric_fasta = numeric_fasta[:, local_start:local_stop]
    desired_fraction = .8*(local_stop-local_start)
    local_coverage = np.sum(local_numeric_fasta != gap_index, axis=1)
    local_sufficient_coverage = local_coverage > desired_fraction
    local_numeric_fasta = local_numeric_fasta[local_sufficient_coverage, :]
    n = local_numeric_fasta.shape[0]
    k = local_numeric_fasta.shape[1]
    display_string = 'Working on block %d, shape (%d, %d), covering %d to %d...'
    display_params = (i, n, k, local_start, local_stop)
    print(display_string % display_params)
    embedded_fasta = np.reshape(embedding[local_numeric_fasta, :], newshape=(n, 5*k))
    reduced_dimensions = TSNE(
        n_components=2, n_iter=5000, random_state=0
    ).fit_transform(embedded_fasta)
    columns = ['d%d' % (j+1) for j in range(2)]
    ids = pd.DataFrame({
        'sequence_id': index[local_sufficient_coverage],
        'block': n*[i]
    }).reset_index()
    rd_df = pd.DataFrame(reduced_dimensions, columns=columns).reset_index()
    df = pd.concat([ids, rd_df], axis=1).drop(columns='index')
    print('...done block %d!' % i)
    return df


def embed_and_reduce_all_dimensions(records, overlap_fraction=.3, ncpu=20):
    index, numeric_fasta = create_numeric_fasta(records)
    info_json = {
        'local_starts': [],
        'local_stops': []
    }
    gap_index = 15
    mean_read_length = np.ceil(np.percentile(np.sum(numeric_fasta != gap_index, axis=1), 25))
    print('Total length %d, mean read length %d ' %(numeric_fasta.shape[1], mean_read_length))
    info_json['mean_read_length'] = mean_read_length
    overlap = np.ceil(overlap_fraction*mean_read_length)
    stride = int(mean_read_length - overlap)
    coverage = np.sum(numeric_fasta != 0, axis=0)
    sufficient_coverage = np.where(coverage >= 20)[0]
    start = sufficient_coverage[0]
    stop = sufficient_coverage[-1]
    nblocks = int(np.ceil((stop-start)/stride))
    starts = np.linspace(start, stop-mean_read_length, nblocks).astype(np.int)
    stops = np.linspace(mean_read_length, stop, nblocks).astype(np.int)
    info_json['local_starts'] = [int(i) for i in starts]
    info_json['local_stops'] = [int(i) for i in stops]
    arguments = zip(nblocks*[numeric_fasta], starts, stops, nblocks*[index], range(nblocks))
    print('Processing %d blocks... pool opened.' % nblocks) 
    pool = Pool(processes=ncpu)
    reduced_dimension_dfs = pool.map(embed_and_reduce_dimensions, arguments)
    print('Done all blocks!')
    pool.close()
    df = pd.concat(reduced_dimension_dfs, axis=0)
    return df, info_json


def embed_and_reduce_all_dimensions_io(input_fasta, output_csv, output_json):
    records = SeqIO.parse(input_fasta, 'fasta')
    df, info_json = embed_and_reduce_all_dimensions(records)
    df.to_csv(output_csv, index_label='id')
    with open(output_json, 'w') as json_file:
        json.dump(info_json, json_file, indent=2)


def cluster_blocks(dr_df, ndim):
    n_blocks = dr_df['block'].max()+1
    dr_df['cluster'] = ''
    columns = ['d%d' % (j+1) for j in range(ndim)]
    for block in range(n_blocks):
        current_block = dr_df.loc[dr_df['block'] == block, columns]

        #clusters = BayesianGaussianMixture(n_components=2).fit(current_block).predict(current_block)
        #dr_df.loc[dr_df['block'] == block, 'cluster'] = ['c%d' % cluster for cluster in clusters]

        #clusters = SpectralClustering(2).fit(current_block)
        #dr_df.loc[dr_df['block'] == block, 'cluster'] = ['c%d' % cluster for cluster in clusters.labels_]

        cluster = AgglomerativeClustering(n_clusters=2, linkage='single')
        cluster.fit_predict(current_block)
        dr_df.loc[dr_df['block'] == block, 'cluster'] = ['c%d' % cluster for cluster in cluster.labels_]
    return dr_df


def cluster_blocks_io(input_csv, output_csv, ndim=2):
    dr_df = pd.read_csv(input_csv)
    cluster_blocks(dr_df, ndim).to_csv(output_csv)
    

def obtain_consensus(cluster_df, reads, info_json):
    n_blocks = len(info_json['local_starts'])
    contigs = []
    for block in range(n_blocks):
        current_block = cluster_df.loc[cluster_df['block'] == block, :]
        clusters = list(set(current_block.loc[:, 'cluster']))
        for cluster in clusters:
            print('Block %d, cluster %s...' % (block, cluster))
            labels_of_reads_in_cluster = current_block.loc[current_block['cluster'] == cluster, 'sequence_id']
            cluster_size = len(labels_of_reads_in_cluster)
            msa = np.array(
                [list(str(reads[label].seq)) for label in labels_of_reads_in_cluster],
                dtype='<U1'
            )
            for i in range(msa.shape[0]):
                  start = np.argmax(msa[i, :] != '-')
                  stop = msa.shape[1] - np.argmax(msa[i,:][::-1] != '-')
                  msa[i,:start] = '|'
                  msa[i,stop:] = '|'
            counts = np.array([
                np.sum(msa == '-', axis=0),
                np.sum(msa == 'A', axis=0),
                np.sum(msa == 'C', axis=0),
                np.sum(msa == 'G', axis=0),
                np.sum(msa == 'T', axis=0)
            ])
            consensus = np.array(list('-ACGT'), dtype='<U1')[np.argmax(counts, axis=0)]
            record_id = 'block-%d_cluster-%s_size-%d' % (block, cluster, cluster_size)
            contig = SeqRecord(Seq(''.join(consensus)), id=record_id, description='')
            contigs.append(contig)
    return contigs


def obtain_consensus_io(input_csv, input_fasta, input_json, output_fasta):
    cluster_df = pd.read_csv(input_csv)
    reads = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
    with open(input_json) as json_file:
        info_json = json.load(json_file)
    contigs = obtain_consensus(cluster_df, reads, info_json)
    SeqIO.write(contigs, output_fasta, 'fasta')


def superreads_to_haplotypes(superreads, max_haplotypes=10):
    blocks = OrderedDict()
    block_sizes = Counter()
    number_of_blocks = 0
    for superread in superreads:

        block = int(superread.name.split('_')[0].split('-')[1])
        size = int(superread.name.split('_')[2].split('-')[1])
        block_sizes[block] += size
        number_of_blocks = max(number_of_blocks, block+1)
        superread_np = np.array(list(str(superread.seq)), dtype='<U1')
        metadata = {
            'bp': superread,
            'np': superread_np,
            'size': size,
            'start': np.argmax(superread_np != '-'),
            'stop': len(superread_np) - np.argmax(superread_np[::-1] != '-') - 1
        }
        if block in blocks:
            blocks[block].append(metadata)
        else:
            blocks[block] = [metadata]

    no_empty_block = True
    records = []
    haplotype_length = len(blocks[0][0]['bp'].seq)
    total_size = 0
    num_haplotypes = 0
    while no_empty_block and num_haplotypes < max_haplotypes:
        first_block_sizes = [item['size'] for item in blocks[0]]
        max_size = max(first_block_sizes)
        max_index = first_block_sizes.index(max_size)
        preferences = [max_index]

        for block_index in range(1, number_of_blocks):
            last_preference_index = preferences[block_index-1]
            last_preference = blocks[block_index-1][last_preference_index]
            block = blocks[block_index]
            scores = np.zeros(len(block))
            stop = last_preference['stop']
            for i, superread in enumerate(block):
                start = superread['start']
                agreement = (superread['np'][start:stop] == last_preference['np'][start:stop]).sum()
                scores[i] = agreement
            max_score = np.argmax(scores)
            preferences.append(int(max_score))
        
        preferred_superreads = [blocks[i][preference] for i, preference in enumerate(preferences)]
        haplotype_size = np.floor(np.mean([sr['size'] for sr in preferred_superreads]))
        new_haplotype = np.array(haplotype_length*['-'], dtype='<U1')
        for i, preference in enumerate(preferences):
            preferred_superread = blocks[i][preference]
            start = preferred_superread['start']
            stop = preferred_superread['stop']
            new_haplotype[start:stop] = preferred_superread['np'][start:stop]
            blocks[i][preference]['size'] -= haplotype_size
            block_sizes[i] -= haplotype_size
            if block_sizes[i] <= 0:
                no_empty_block = False
        new_record = SeqRecord(
            Seq(''.join(new_haplotype)),
            id='haplotype-%d_size-%d' % (len(records)+1, haplotype_size),
            description=''
        )
        num_haplotypes += 1
        records.append(new_record)
        total_size += haplotype_size
    for record in records:        
        size = int(record.id.split('_')[1].split('-')[1])
        record.id += '_freq-%f' % (size/total_size)
    return records


def superreads_to_haplotypes_io(input_fasta, output_fasta):
    superreads = SeqIO.parse(input_fasta, 'fasta')
    haplotypes = superreads_to_haplotypes(superreads)
    SeqIO.write(haplotypes, output_fasta, 'fasta')


def evaluate(input_haplotypes, input_truth, output_json):
    haplotypes = SeqIO.parse(input_haplotypes, 'fasta')
    truth = SeqIO.parse(input_truth, 'fasta')
    haplotype_index, numeric_haplotypes = create_numeric_fasta(haplotypes)
    truth_index, numeric_truth = create_numeric_fasta(truth)
    full_numeric = np.vstack([numeric_truth, numeric_haplotypes])
    counts = np.array([
        np.sum(full_numeric == 15, axis=0),
        np.sum(full_numeric == 0, axis=0),
        np.sum(full_numeric == 1, axis=0),
        np.sum(full_numeric == 2, axis=0),
        np.sum(full_numeric == 3, axis=0)
    ])
    discordant = np.sum(counts == 0, axis=0) != 4
    get_headers_as_strings = lambda index: [str(header) for header in index]
    headers = get_headers_as_strings(truth_index)+get_headers_as_strings(haplotype_index)
    n_headers = len(headers)
    discordance_matrix = [n_headers*[0] for i in range(n_headers)]
    for i in range(n_headers):
        for j in range(n_headers):
            discordance = np.sum(full_numeric[i,:] != full_numeric[j,:])
            discordance_matrix[i][j] = int(discordance)
    output = {
        'number_of_discordant_sites': int(np.sum(discordant)),
        'headers': headers,
        'discordance_matrix': discordance_matrix
    }
    with open(output_json, 'w') as json_file:
        json.dump(output, json_file, indent=2)

