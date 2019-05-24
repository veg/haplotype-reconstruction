import itertools as it
from multiprocessing import Pool

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import binom
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .mapped_reads import MappedReads
from .utils import read_reference_start_and_end


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return it.zip_longest(*args, fillvalue=fillvalue)


def partial_covariation_test(arguments):
    pysam_path, pairs, i = arguments
    pairs, pairs_for_max, pairs_for_min = it.tee(
        it.filterfalse(lambda x: x is None, pairs),
        3
    )
    print('   ...processing block %d...' % i)
    pair_min = min([min(pair[0], pair[1]) for pair in pairs_for_min])
    pair_max = max([max(pair[0], pair[1]) for pair in pairs_for_max])
    mr = MappedReads(pysam_path)
    fasta_window = mr.sam_window_to_fasta(pair_min, pair_max+1)
    results = []
    for col_i, col_j in pairs:
        idx_i = col_i - pair_min
        idx_j = col_j - pair_min
        i_char_counts = pd.Series(
            fasta_window[:, idx_i]
        ).value_counts().drop('~', errors='ignore')
        i_chars = i_char_counts.index[i_char_counts > 0]

        j_char_counts = pd.Series(
            fasta_window[:, idx_j]
        ).value_counts().drop('~', errors='ignore')
        j_chars = j_char_counts.index[j_char_counts > 0]

        content_i = fasta_window[:, idx_i] != '~'
        content_j = fasta_window[:, idx_j] != '~'
        valid = content_i & content_j
        if valid.sum() == 0:
            continue
        for i_char, j_char in it.product(i_chars, j_chars):
            equals_i = fasta_window[valid, idx_i] == i_char
            equals_j = fasta_window[valid, idx_j] == j_char
            X_11 = (equals_i & equals_j).sum()
            X_21 = (~equals_i & equals_j).sum()
            X_12 = (equals_i & ~equals_j).sum()
            X_22 = (~equals_i & ~equals_j).sum()
            
            table = [
                [X_11 , X_12],
                [X_21 , X_22]
            ]
            _, p_value = fisher_exact(table)
            results.append((col_i, col_j, i_char, j_char, p_value))
    mr.close()
    print('   ...done block %d!' % i)
    return pd.DataFrame(
        results, columns=('col_i', 'col_j', 'i_char', 'j_char', 'p_value')
    )


class ErrorCorrection:
    def __init__(self, pysam_alignment, all_fe_tests=None, error_threshold=1e-3):
        self.pysam_alignment = pysam_alignment
        self.reference_length = pysam_alignment.header['SQ'][0]['LN']
        self.number_of_reads = 0
        for read in pysam_alignment.fetch():
            self.number_of_reads += 1

        if all_fe_tests:
            self.all_fe_tests = pd.read_csv(all_fe_tests)
        else:
            self.all_fe_test = None
        self.covarying_sites = None
        self.pairs = None
        self.nucleotide_counts = None
        self.error_threshold = error_threshold

    @staticmethod
    def read_count_data(read):
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
                segments.append(unaligned_sequence[position: position + stride])
                position += stride
            elif insertion:
                position += stride
            elif deletion:
                if len(segments) > 0 and i < number_of_cigar_tuples:
                    segments.append(stride * '-')
        sequence = np.concatenate([list(segment) for segment in segments])
        return sequence, positions

    def get_nucleotide_counts(self):
        if not self.nucleotide_counts is None:
            return self.nucleotide_counts
        print('Calculating nucleotide counts...')
        characters = ['A', 'C', 'G', 'T', '-']
        counts = np.zeros((self.reference_length, 5))
        for read in self.pysam_alignment.fetch():
            sequence, positions = self.read_count_data(read)
            for character_index, character in enumerate(characters):
                rows = positions[sequence == character]
                counts[rows, character_index] += 1

        df = pd.DataFrame(counts, columns=characters)
        zeros = lambda character: (df[character] == 0).astype(np.int)
        zero_cols = zeros('A') + zeros('C') + zeros('G') + zeros('T')
        df['interesting'] = zero_cols < 3
        df['nucleotide_max'] = df[['A', 'C', 'G', 'T']].max(axis=1)
        df['coverage'] = df[['A', 'C', 'G', 'T']].sum(axis=1)
        df['consensus'] = '-'
        for character in characters[:-1]:
            df.loc[df['nucleotide_max'] == df[character], 'consensus'] = character
        self.nucleotide_counts = df
        return df

    def consensus(self):
        consensus_sequence = Seq(''.join(self.nucleotide_counts['consensus']))
        record = SeqRecord(consensus_sequence, id='consensus', description='')
        return record

    def get_pairs(self):
        if self.pairs:
            return self.pairs
        counts = self.get_nucleotide_counts()
        interesting = counts.index[counts.interesting]
        max_read_length = max([
            read.infer_query_length()
            for read in self.pysam_alignment.fetch()
        ])
        pairs = list(filter(
            lambda pair: pair[1] - pair[0] <= max_read_length,
            it.combinations(interesting, 2)
        ))
        self.pairs = pairs
        return pairs
    
    def full_covariation_test(self, threshold=20, stride=10000,
            ncpu=24, block_size=250):
        if self.covarying_sites:
            return self.covarying_sites
        pairs = self.get_pairs()
        arguments = [
            (self.pysam_alignment.filename, group, i)
            for i, group in enumerate(grouper(pairs, block_size))
        ]
        n_pairs = len(pairs)
        n_blocks = len(arguments)
        print('Processing %d blocks of %d pairs...' % (n_blocks, n_pairs))
        pool = Pool(processes=ncpu)
        result_dfs = pool.map(partial_covariation_test, arguments)
        pool.close()
        self.all_fe_tests = pd.concat(result_dfs).sort_values(by='p_value')
        print('...done!')

    def multiple_testing_correction(self, fdr=.001):
        print('Performing multiple testing correction...')
        m = len(self.all_fe_tests)
        self.all_fe_tests['bh'] = self.all_fe_tests['p_value'] <= fdr*np.arange(1, m+1)/m
        cutoff = (1-self.all_fe_tests['bh']).to_numpy().nonzero()[0][0]
        covarying_sites = np.unique(
            np.concatenate([
                self.all_fe_tests['col_i'].iloc[:cutoff],
                self.all_fe_tests['col_j'].iloc[:cutoff]
            ])
        )
        covarying_sites.sort()
        self.covarying_sites = covarying_sites
        self.nucleotide_counts.loc[:, 'covarying'] = False
        self.nucleotide_counts.loc[covarying_sites, 'covarying'] = True
        return covarying_sites

    def get_covarying_errors(self):
        nucleotide_counts = self.get_nucleotide_counts()
        summary = nucleotide_counts.loc[
            nucleotide_counts.covarying==False,
            ['nucleotide_max', 'coverage']
        ].sum()
        total_coverage = summary['coverage']
        total_consensus = summary['nucleotide_max']
        error_rate = np.abs(total_coverage - total_consensus) / total_consensus
        nucleotide_counts.loc[:, 'n_error'] = nucleotide_counts.loc[:, 'coverage'].apply(
            lambda count: binom.ppf(1-self.error_threshold, count, error_rate)
        )
        nucleotide_counts.loc[:, 'site'] = nucleotide_counts.index
        nucleotide_counts.loc[:, 'covarying'] = False
        nucleotide_counts.loc[self.covarying_sites, 'covarying'] = True
        site_counts = nucleotide_counts.loc[nucleotide_counts['covarying'], :] \
            .melt(id_vars=['n_error', 'site'], value_vars=['A', 'C', 'G', 'T'])
        covarying_values = site_counts['value']
        covarying_counts = site_counts['n_error']
        significant = (covarying_values <= covarying_counts) & (covarying_values > 0)
        covarying_errors = site_counts.loc[significant, :] \
            .sort_values(by='site') \
            .reset_index(drop=True)
        self.covarying_errors = covarying_errors
        return covarying_errors

    def get_covarying_correction(read):
        pass

    def get_all_covarying_corrections(self):
        covarying_errors = self.get_covarying_errors()
        self.read_information = read_reference_start_and_end(
            self.pysam_alignment, self.covarying_sites
        )
        corrections = {}
        for read in self.pysam_alignment.fetch():
            ce_start = (covarying_errors['site'] >= read.reference_start).idxmax()
            ce_end = (covarying_errors['site'] >= read.reference_end).idxmax()-1
            sequence, position = self.read_count_data(read)
            unshifted_indices = covarying_errors.loc[ce_start:ce_end, 'site']
            shifted_indices = unshifted_indices - read.reference_start
            cv_sequence = sequence[shifted_indices]
            cv_errors = covarying_errors.loc[ce_start:ce_end, 'variable']
            if (cv_sequence == cv_errors).any():
                corrections[read.query_name] = self.get_covarying_correction(read)

    def corrected_reads(self, **kwargs):
        nucleotide_counts = self.get_nucleotide_counts()
        self.full_covariation_test(**kwargs)
        covarying_sites = self.multiple_testing_correction()
        for read in self.pysam_alignment.fetch():
            sequence, _ = self.read_count_data(read)
            intraread_covarying_sites = covarying_sites[
                (covarying_sites >= read.reference_start) &
                (covarying_sites < read.reference_end)
            ]
            mask = np.ones(len(sequence), np.bool)
            mask[intraread_covarying_sites - read.reference_start] = False
            local_consensus = nucleotide_counts.consensus[
                read.reference_start: read.reference_end
            ]
            sequence[mask] = local_consensus[mask]
            
            corrected_read = pysam.AlignedSegment()
            corrected_read.query_name = read.query_name
            corrected_read.query_sequence = ''.join(sequence)
            corrected_read.flag = read.flag
            corrected_read.reference_id = 0
            corrected_read.reference_start = read.reference_start
            corrected_read.mapping_quality = read.mapping_quality
            corrected_read.cigar = [(0, len(sequence))]
            corrected_read.next_reference_id = read.next_reference_id
            corrected_read.next_reference_start = read.next_reference_start
            corrected_read.template_length = read.template_length
            corrected_read.query_qualities = pysam.qualitystring_to_array(
                len(sequence) * '<'
            )
            corrected_read.tags = read.tags
            yield corrected_read

    def write_corrected_reads(self, output_bam_filename):
        output_bam = pysam.AlignmentFile(
            output_bam_filename, 'wb', header=self.pysam_alignment.header
        )
        for read in self.corrected_reads():
            output_bam.write(read)
        output_bam.close()

    def __del__(self):
        self.pysam_alignment.close()

