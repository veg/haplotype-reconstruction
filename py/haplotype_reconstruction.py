import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.manifold import TSNE


def embed_and_reduce_dimensions(records, blocks, ndim=2):
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
    
    if blocks:
        mean_read_length = np.sum(numeric_fasta != 0, axis=1).mean()
        coverage = np.sum(numeric_fasta != 0, axis=0)
        sufficient_coverage = np.where(coverage >= 20)[0]
        start = sufficient_coverage[0]
        stop = sufficient_coverage[-1]
        n_blocks = np.int(np.ceil((stop-start)/mean_read_length))
        intervals = np.ceil(np.linspace(start, stop, n_blocks+1)).astype(np.int)
        intervals[-1] += 1

        reduced_dimension_dfs = []
        for i in range(n_blocks):
            local_start = intervals[i]
            local_stop = intervals[i+1]
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
        return pd.concat(reduced_dimension_dfs, axis=0)
    else:
        n = numeric_fasta.shape[0]
        k = numeric_fasta.shape[1]
        embedded_fasta = np.reshape(embedding[numeric_fasta, :], newshape=(n, 5*k))
        reduced_dimensions = TSNE(
            n_components=ndim, n_iter=5000, random_state=0
        ).fit_transform(embedded_fasta)
        columns = ['d%d' % (j+1) for j in range(ndim)]
        return pd.DataFrame(
            reduced_dimensions,
            index=index,
            columns=columns
        )


def embed_and_reduce_dimensions_io(input_fasta, output_csv, blocks=True, ndim=2):
    records = SeqIO.parse(input_fasta, 'fasta')
    df = embed_and_reduce_dimensions(records, blocks, ndim)
    df.to_csv(output_csv, index_label='id')

