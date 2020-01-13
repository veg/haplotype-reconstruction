import json

from sklearn.manifold import SpectralEmbedding
import numpy as np
from Bio import SeqIO
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import pysam


characters = ['A', 'C', 'G', 'T', '-']


def single_read_count_data(read):
    sequence_length = np.array([
        cigar_tuple[1]
        for cigar_tuple in read.cigartuples
        if cigar_tuple[0] != 1
    ]).sum()
    first_position = read.reference_start
    last_position = first_position + sequence_length
    positions = np.arange(first_position, last_position)
    segments = []
    number_of_cigar_tuples = len(read.cigartuples)
    unaligned_sequence = read.query_alignment_sequence
    position = 0
    for i, cigar_tuple in enumerate(read.cigartuples):
        action = cigar_tuple[0]
        stride = cigar_tuple[1]
        match = action == 0
        insertion = action == 1
        deletion = action == 2
        if match:
            segments.append(
                unaligned_sequence[position: position + stride]
                )
            position += stride
        elif insertion:
            position += stride
        elif deletion:
            if len(segments) > 0 and i < number_of_cigar_tuples:
                segments.append(stride * '-')
    sequence = np.concatenate([list(segment) for segment in segments])
    return sequence, positions


def all_read_count_data(alignment):
    reference_length = alignment.header['SQ'][0]['LN']
    counts = np.zeros((reference_length, 5))
    for read in alignment.fetch():
        sequence, positions = single_read_count_data(read)
        for character_index, character in enumerate(characters):
            rows = positions[sequence == character]
            counts[rows, character_index] += 1
    return counts



def supplementary_info(row):
    result = row.loc[['A', 'C', 'G', 'T']] \
        .sort_values(ascending=False)
    result.index = ['c1', 'c2', 'c3', 'c4']
    return result.append(pd.Series(
        result.values/row['coverage'] if row['coverage'] else 0,
        index=['f1', 'f2', 'f3', 'f4']
    ))


def zeros(df, character):
    return (df[character] == 0).astype(np.int)


def site_table(alignment):
    counts = all_read_count_data(alignment)
    df = pd.DataFrame(counts, columns=characters)
    col_0s = zeros(df, 'A') + zeros(df, 'C') + zeros(df, 'G') + zeros(df, 'T')
    df['interesting'] = col_0s < 3
    df['nucleotide_max'] = df[['A', 'C', 'G', 'T']].max(axis=1)
    df['coverage'] = df[['A', 'C', 'G', 'T']].sum(axis=1)
    df['consensus'] = '-'
    for character in characters[:-1]:
        consensus_agreement = df['nucleotide_max'] == df[character]
        df.loc[consensus_agreement, 'consensus'] = character
    return pd.concat([df, df.apply(supplementary_info, axis=1)], axis=1)


def extract_label(query_name):
    return '-'.join(query_name.split('.')[:4]) if not '+' in query_name else 'AR'


def get_covarying_sites(alignment, threshold=.01, end_correction=10):
    nucleotide_counts = site_table(alignment)
    above_threshold = (
        nucleotide_counts
        .loc[:, ['f1', 'f2', 'f3', 'f4']] > threshold
    ).sum(axis=1)
    all_integers = np.arange(0, len(above_threshold))
    covarying_sites = all_integers[above_threshold > 1]
    number_of_sites = len(nucleotide_counts)
    after_head_correction = covarying_sites > end_correction
    final_site = number_of_sites - end_correction
    before_tail_correction = covarying_sites < final_site
    desired = after_head_correction & before_tail_correction
    return covarying_sites[desired]


def read_reference_start_and_end(alignment, site_boundaries):
    read_information = pd.DataFrame(
        [
            (read.reference_start, read.reference_end)
            for read in alignment.fetch()
        ],
        columns=['reference_start', 'reference_end']
    )
    read_information['covarying_start'] = np.searchsorted(
        site_boundaries, read_information['reference_start']
    )
    read_information['covarying_end'] = np.searchsorted(
        site_boundaries, read_information['reference_end']
    )
    return read_information


def admission(minimum_weight):
    def comparator(pair):
        return pair[1][0] >= minimum_weight
    return comparator


def obtain_superreads(alignment, covarying_sites, minimum_weight=3):
    read_information = read_reference_start_and_end(
        alignment, covarying_sites
    )
    read_groups = {}
    all_reads = list(alignment.fetch())
    for i, read in enumerate(all_reads):
        row = read_information.iloc[i, :]
        key = (row['covarying_start'], row['covarying_end'])
        if key in read_groups:
            read_groups[key].append(i)
        else:
            read_groups[key] = [i]
    all_superreads = []
    superread_index = 0
    for covarying_boundaries, read_group in read_groups.items():
        if covarying_boundaries[0] == covarying_boundaries[1]:
            continue
        superreads = {}
        for read_index in read_group:
            read = all_reads[read_index]
            label = extract_label(read.query_name)
            covarying_sites_in_read = covarying_sites[
                covarying_boundaries[0]: covarying_boundaries[1]
            ]
            value_at_covarying_sites = ''.join(
                [
                    read.query[triplet[0]].upper()
                    for triplet in read.get_aligned_pairs(True)
                    if triplet[1] in covarying_sites_in_read
                ]
            )
            has_ar = 1 if '+' in read.query_name else 0
            if value_at_covarying_sites in superreads:
                superreads[value_at_covarying_sites][0] += 1
                superreads[value_at_covarying_sites][1] += has_ar
                if not label in superreads[value_at_covarying_sites][2]:
                     superreads[value_at_covarying_sites][2][label] = 0
                superreads[value_at_covarying_sites][2][label] += 1
            else:
                superreads[value_at_covarying_sites] = [1, has_ar, {label: 1}]
        admissible_superreads = list(filter(admission(minimum_weight), superreads.items()))
        total_weight = sum([
            superread[1][0] for superread in admissible_superreads
        ])
        for vacs, weight in admissible_superreads:
            all_superreads.append({
                'index': superread_index,
                'vacs': vacs,
                'weight': weight[0],
                'frequency': weight[0]/total_weight,
                'ar': weight[1],
                'ar_frequency': weight[1]/weight[0],
                'cv_start': int(covarying_boundaries[0]),
                'cv_end': int(covarying_boundaries[1]),
                'composition': weight[2],
                'discarded': False
            })
            superread_index += 1
    return all_superreads


def superread_cv_filter(superreads, min_cv_start, max_cv_end):
    def cv_filter(sr):
        starts_after = sr['cv_start'] >= min_cv_start
        ends_before = sr['cv_end'] < max_cv_end
        return starts_after and ends_before
    return list(filter(cv_filter, superreads))


def get_score_matrix(superreads, min_cv_start, max_cv_end,
        minimum_agreement=0, power=1):
    n_sr = len(superreads)
    scores = []
    rows = []
    cols = []
    for i, superread_i in enumerate(superreads):
        for j, superread_j in enumerate(superreads):
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
                i_sequence = np.array(superread_i['vacs'][i_start: i_end], dtype='<U1')
                j_sequence = np.array(superread_j['vacs'][j_start: j_end], dtype='<U1')
                agreement = (i_sequence == j_sequence).sum()
                if agreement > minimum_agreement:
                    rows.append(i)
                    cols.append(j)
                    weight = min(superread_i['weight'], superread_j['weight'])
                    scores.append(agreement**power * weight)
    return scipy.sparse.csr_matrix((scores, (rows, cols)), shape=(n_sr, n_sr))


def perform_spectral_embedding(superreads, min_cv_start, max_cv_end):
    X = get_score_matrix(superreads, min_cv_start, max_cv_end)
    return SpectralEmbedding(
        n_components=2,
        random_state=0,
        affinity='precomputed'
    ).fit_transform(X.todense())


def get_labels(superreads):
    return [max(sr['composition'].items(), key=lambda x: x[1])[0] for sr in superreads]


def sc_covarying_sites_io(bam_path, json_path):
    alignment = pysam.AlignmentFile(bam_path, 'rb')
    covarying_sites = get_covarying_sites(alignment)
    covarying_sites_json = [int(site) for site in covarying_sites]
    with open(json_path, 'w') as json_file:
        json.dump(covarying_sites_json, json_file)


def sc_superread_io(bam_path, covarying_path, superread_path):
    with open(covarying_path) as json_file:
        covarying_sites = np.array(json.load(json_file), dtype=np.int)
    alignment = pysam.AlignmentFile(bam_path, 'rb')

    superreads = obtain_superreads(alignment, covarying_sites)
    with open(superread_path, 'w') as json_file:
        json.dump(superreads, json_file, indent=2)


def sc_embedding_io(superread_path, embedding_path, min_cv_start, max_cv_end):
    min_cv_start = int(min_cv_start)
    max_cv_end = int(max_cv_end)
    with open(superread_path) as json_file:
        superreads = superread_cv_filter(
            json.load(json_file),
            min_cv_start,
            max_cv_end
        )
    embedding = perform_spectral_embedding(superreads, min_cv_start, max_cv_end)
    df = pd.DataFrame({
        'x': embedding[:, 0],
        'y': embedding[:, 1],
        'label': get_labels(superreads)
    })
    df.to_csv(embedding_path)


def sc_srfasta_io(input_cvs, input_srdata, output_fasta):
    with open(input_cvs) as json_file:
        cvs = json.load(json_file)
    with open(input_srdata) as json_file:
        srdata = json.load(json_file)
    outfile = open(output_fasta, 'w')
    n_cvs = len(cvs)
    for sr in srdata:
        outfile.write('>superread-%d\n' % sr['index'])
        seq = sr['cv_start']*'-' + sr['vacs'] + (n_cvs - sr['cv_end'])*'-'
        outfile.write(seq + '\n')
    outfile.close()


def sc_truthcvs_io():
    truth = list(SeqIO.parse(input_fasta))
    for record in truth:
        record.seq = ''.join([record.seq[i] for i in cvs])
    SeqIO.write(truth, output_sr, 'fasta')
    SeqIO.write(truth, output_truth, 'fasta')
