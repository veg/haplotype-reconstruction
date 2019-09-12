import subprocess
import json
import os

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_orf(input_genome, output_genome, orf):
    orf = int(orf)
    record = SeqIO.read(input_genome, 'fasta')
    record.seq = record.seq[orf:]
    SeqIO.write(record, output_genome, 'fasta')


def backtranslate(input_nucleotide, input_protein, output_codon):
    nucleotides = SeqIO.parse(input_nucleotide, 'fasta')
    proteins = SeqIO.parse(input_protein, 'fasta')
    codons = []

    for protein_record, nucleotide_record in zip(proteins, nucleotides):
        i = 0
        codon_list = []
        for character in protein_record.seq:
            if character != '-':
                codon = str(nucleotide_record.seq[3*i:3*i+3])
                codon_list.append(codon)
                i += 1
            else:
                codon_list.append('---')
        codon_record = SeqRecord(
            Seq(''.join(codon_list)),
            id=protein_record.id,
            description=protein_record.description
        )
        codons.append(codon_record)

    SeqIO.write(codons, output_codon, 'fasta')


def select_simulated_gene(dataset, gene, output):
    aligned_filename = "output/simulation/%s/aligned_%s_orf-%d_codon.fasta"
    nucleotide_genome_filename = "output/simulation/%s/genome.fasta" % dataset
    nucleotide_genome = SeqIO.read(nucleotide_genome_filename, 'fasta')
    max_percent_identity = 0
    for i in range(3):
        non_gaps = 0
        matches = 0
        codon_list = []
        records = SeqIO.parse(aligned_filename % (dataset, gene, i), 'fasta')
        translated_genome = next(records)
        reference = next(records)
        genome_i = 0
        for j in range(len(reference)):
            if reference[j] != '-':
                non_gaps += 1
                codon = str(nucleotide_genome[3*genome_i+i:3*genome_i+i+3].seq)
                codon_list.append(codon)
                if reference[j] == translated_genome[j]:
                    matches += 1
            if translated_genome[j] != '-':
                genome_i += 1
        percent_identity = matches/non_gaps
        if percent_identity > max_percent_identity:
            max_percent_identity = percent_identity
            desired_codons = ''.join(codon_list)
    record = SeqRecord(
        Seq(desired_codons).ungap('-'),
        id=nucleotide_genome.id,
        description=gene
    )
    SeqIO.write(record, output, 'fasta')


def write_abayesqr_config(sam_filename, reference_filename, output):
    config_string = ("""filename of reference sequence (FASTA) : %s
filname of the aligned reads (sam format) : %s
paired-end (1 = true, 0 = false) : 0
SNV_thres : 0.01
reconstruction_start : 1
reconstruction_stop: 1300
min_mapping_qual : 20
min_read_length : 50
max_insert_length : 250
characteristic zone name : test
seq_err (assumed sequencing error rate(%%)) : 0.1
MEC improvement threshold : 0.0395 """ % (reference_filename, sam_filename))
    with open(output, 'w') as config_file:
        config_file.write(config_string)


def parse_abayesqr_output(input_text, output_fasta):
    with open(input_text) as input_file:
        lines = input_file.readlines()
    records = []
    for i, line in enumerate(lines):
        if i % 2 == 0:
            freq = float(line.split()[-1])
            number = int(i/2)+1
            header = 'haplotype-%d_freq-%f' % (number, freq)
        if i % 2 == 1:
            seq = Seq(line.strip())
            record = SeqRecord(seq, id=header, description='')
            records.append(record)
    SeqIO.write(records, output_fasta, 'fasta')


def pairwise_distance_csv(fasta_filename, csv_filename):
    records = list(SeqIO.parse(fasta_filename, 'fasta'))
    np_seqs = np.array(
        [list(str(record.seq)) for record in records],
        dtype='<U1'
    )
    first_records = []
    second_records = []
    distances = []
    for i in range(len(records)):
        for j in range(len(records)):
            first_records.append(records[i].id)
            second_records.append(records[j].id)
            distance = (np_seqs[i, :] != np_seqs[j, :]).sum()
            distances.append(distance)
    pd.DataFrame({
        'first_record': first_records,
        'second_record': second_records,
        'distance': distances,
    }).to_csv(csv_filename)


def add_subtype_information(input_csv, output_csv):
    df = pd.read_csv(input_csv)
    df['Subtype1'] = df['ID1'].apply(lambda row: row.split('.')[0])
    df['Subtype2'] = df['ID2'].apply(lambda row: row.split('.')[0])
    df.to_csv(output_csv, index=False)


def extract_5vm_truth(input_fasta, reference_path, output_path):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    aligned_sequences = []
    for sequence in sequences:
        output_dir = os.path.join("output", "FiveVirusMixIllumina_1")
        sequence_path = os.path.join(output_dir, sequence.id, "ref.fasta")
        alignment_path = os.path.join(output_dir, sequence.id, "aligned.fasta")
        SeqIO.write(sequence, alignment_path, "fasta")
        command = [
            "water", "-asequence", sequence_path, "-bsequence",
            reference_path, "-gapopen", "10.0", "-gapextend", ".5", "-aformat",
            "fasta", "-outfile", alignment_path
        ]
        subprocess.run(command)
        aligned_sequence = list(SeqIO.parse(alignment_path, "fasta"))[0]
        aligned_sequence.seq = aligned_sequence.seq.ungap('-')
        aligned_sequences.append(aligned_sequence)
    sequence_length = min([len(record.seq) for record in aligned_sequences])
    for record in aligned_sequences:
        record.seq = record.seq[:sequence_length]
    SeqIO.write(aligned_sequences, output_path, "fasta")


def covarying_truth(
        input_computed, input_actual, input_reference, output_json
        ):
    reference = SeqIO.read(input_reference, 'fasta')
    rl = len(reference.seq)
    with open(input_computed) as input_file:
        cvs = json.load(input_file)
    with open(input_actual) as input_file:
        true_cvs = json.load(input_file)

    tp = []
    fp = []
    tn = []
    fn = []

    for i in range(rl):
        if i in true_cvs and i in cvs:
            tp.append(i)
        elif i in true_cvs and i not in cvs:
            fn.append(i)
        elif i not in true_cvs and i in cvs:
            fp.append(i)
        elif i not in true_cvs and i not in cvs:
            tn.append(i)
    precision = len(tp)/(len(tp)+len(fp))
    recall = len(tp)/(len(tp)+len(fn))
    result = {
        'true_positives': tp,
        'true_negative': tn,
        'false_positives': fp,
        'false_negatives': fn,
        'precision': precision,
        'recall': recall
    }
    with open(output_json, 'w') as output_file:
        json.dump(result, output_file, indent=2)


def downsample_bam(input_bam_path, output_bam_path, downsample_amount):
    downsample_percentage = 1 - int(downsample_amount) / 100
    input_bam = pysam.AlignmentFile(input_bam_path, 'rb')
    number_of_reads = input_bam.count()
    downsample_number = np.ceil(downsample_percentage * number_of_reads) \
        .astype(np.int)
    np.random.seed(1)
    downsample_indices = np.random.choice(
        number_of_reads, downsample_number, replace=False
    )
    downsample_indices.sort()
    downsample_index = 0
    output_bam = pysam.AlignmentFile(
        output_bam_path, 'wb', header=input_bam.header
    )
    for i, read in enumerate(input_bam.fetch()):
        if i == downsample_indices[downsample_index]:
            output_bam.write(read)
            downsample_index += 1
        if downsample_index == len(downsample_indices):
            break
    output_bam.close()
    pysam.index(output_bam_path)
    input_bam.close()


def pluck_record(input_fasta_path, output_fasta_path, record):
    all_records = SeqIO.parse(input_fasta_path, 'fasta')
    desired_record = SeqIO.to_dict(all_records)[record]
    SeqIO.write(desired_record, output_fasta_path, 'fasta')


def single_mapping_dataset(bam_path, ref_path, output_path):
    bam = pysam.AlignmentFile(bam_path)
    ref = SeqIO.read(ref_path, 'fasta')
    percent_identity = np.zeros(bam.mapped, dtype=np.float)
    differences = np.zeros(bam.mapped, dtype=np.float)
    number_of_aligned_pairs = np.zeros(bam.mapped, dtype=np.float)
    for i, read in enumerate(bam.fetch()):
        aligned_pairs = read.get_aligned_pairs(matches_only=True)
        aligned_query = np.array([
            read.query[pair[0]] for pair in aligned_pairs
        ], dtype='<U1')
        aligned_reference = np.array([
            ref[pair[1]] for pair in aligned_pairs
        ], dtype='<U1')
        agreement = (aligned_query == aligned_reference).sum()
        number_of_aligned_pairs[i] = len(aligned_pairs)
        differences[i] = number_of_aligned_pairs[i] - agreement
        percent_identity[i] = agreement/number_of_aligned_pairs[i]

    quality = np.array([
        read.mapping_quality for read in bam.fetch()
    ], dtype=np.int)
    query_length = np.array([
        read.query_length for read in bam.fetch()
    ], dtype=np.int)
    result = pd.DataFrame({
        'mapping_quality': quality,
        'differences': differences,
        'number_of_aligned_pairs': number_of_aligned_pairs,
        'percent_identity': percent_identity,
        'query_length': query_length
    }, index=[read.query_name for read in bam.fetch()])
    result.to_csv(output_path, index_label='read_id')


def full_fvm_mapping_dataset(dataset_paths, output_csv_path):
    all_datasets = list(map(
        lambda path: pd.read_csv(path, index_col='read_id'),
        dataset_paths
    ))
    for dataset_path, dataset in zip(dataset_paths, all_datasets):
        dataset_name = dataset_path.split('/')[-2]
        dataset['reference'] = dataset_name
    pd.concat(all_datasets, axis=0, sort=False, ignore_index=True) \
        .to_csv(output_csv_path)


def true_covarying_kmers(input_fasta, input_json, output_csv, k):
    k = int(k)
    records = np.array([
        list(record.seq)
        for record in SeqIO.parse(input_fasta, 'fasta')
    ], dtype='<U1')
    data = {
        **{'index_%d' % i: [] for i in range(k)},
        **{'character_%d' % i: [] for i in range(k)}
    }
    with open(input_json) as json_file:
        covarying_sites = np.array(json.load(json_file), dtype=np.int)
    for i in range(len(covarying_sites) - k):
        covarying_indices = covarying_sites[i:i+k]
        covarying_kmers = set()
        for row_index in range(records.shape[0]):
            covarying_kmer = ''.join(records[row_index, covarying_indices])
            covarying_kmers.add(covarying_kmer)
        for covarying_kmer in list(covarying_kmers):
            for i in range(k):
                data['index_%d' % i].append(covarying_indices[i])
                data['character_%d' % i].append(covarying_kmer[i])
    pd.DataFrame(data).to_csv(output_csv, index=False)
