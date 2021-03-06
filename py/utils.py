import subprocess
import json
import os
import csv
from time import time

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
    search_term = 'quasispecies'
    for i in range(len(records)):
        for j in range(len(records)):
            if records[j].name[: len(search_term)] == search_term:
                continue
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


def extract_truth(
        input_fasta, reference_path, dataset, reference, output_path,
        output_json_path
        ):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    aligned_sequences = []
    output_dir = os.path.join("output", "truth", dataset)
    tmp_dir = os.path.join(
        output_dir, "truth-%s-%s-temp" % (dataset, reference)
    )
    os.mkdir(tmp_dir)
    for sequence in sequences:
        sequence_path = os.path.join(tmp_dir, "ref.fasta")
        alignment_path = os.path.join(tmp_dir, "aligned.fasta")
        SeqIO.write(sequence, sequence_path, "fasta")
        command = [
            "water", "-asequence", sequence_path, "-bsequence",
            reference_path, "-gapopen", "10.0", "-gapextend", ".5", "-aformat",
            "fasta", "-outfile", alignment_path
        ]
        subprocess.run(command)
        aligned_sequence = list(SeqIO.parse(alignment_path, "fasta"))[0]
        aligned_sequence.seq = aligned_sequence.seq.ungap('-')
        aligned_sequences.append(aligned_sequence)
        os.remove(sequence_path)
        os.remove(alignment_path)
    os.rmdir(tmp_dir)
    sequence_length = min([len(record.seq) for record in aligned_sequences])
    for record in aligned_sequences:
        record.seq = record.seq[:sequence_length]
    SeqIO.write(aligned_sequences, output_path, "fasta")

    pairwise_distances = []
    for i in range(len(aligned_sequences)):
        first_sequence = aligned_sequences[i]
        first_np = np.array(list(first_sequence.seq), dtype='<U1')
        for j in range(i+1, len(aligned_sequences)):
            second_sequence = aligned_sequences[j]
            second_np = np.array(list(second_sequence.seq), dtype='<U1')
            disagreement = int((first_np != second_np).sum())
            pairwise_distances.append({
                'sequenceA': first_sequence.name,
                'sequenceB': second_sequence.name,
                'disagreement': disagreement
            })
    with open(output_json_path, 'w') as json_file:
        json.dump(pairwise_distances, json_file, indent=2)

def covarying_truth(
        input_computed, input_truth, input_reference, output_json
        ):
    reference = SeqIO.read(input_reference, 'fasta')
    rl = len(reference.seq)
    with open(input_computed) as input_file:
        cvs = json.load(input_file)
    with open(input_truth) as input_file:
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


def restrict_fasta_to_cvs(input_fasta, input_cvs, output_fasta):
    with open(input_cvs) as json_file:
        cvs = json.load(json_file)
    records = list(SeqIO.parse(input_fasta, 'fasta'))
    for record in records:
        record.seq = Seq(''.join([record.seq[site] for site in cvs]))
    SeqIO.write(records, output_fasta, 'fasta')


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


def kmers_in_reads(input_bam, input_csv, output_csv, k):
    k = int(k)
    bam = pysam.AlignmentFile(input_bam)
    df = pd.read_csv(input_csv)
    df['support'] = np.zeros(len(df), dtype=np.int)
    for read in bam.fetch():
        starts_after = df.index_0 >= read.reference_start
        ends_before = df['index_%d' % (k-1)] <= read.reference_end
        relevent_kmers = df.loc[starts_after & ends_before, :]
        for i, row in relevent_kmers.iterrows():
            inds = list(row[['index_%d' % i for i in range(k)]])
            vacs = ''.join([
                read.query[pair[0]]
                for pair in read.get_aligned_pairs(matches_only=True)
                if pair[1] in inds
            ])
            kmer = ''.join(row[['character_%d' % i for i in range(k)]])
            if vacs == kmer:
                df.loc[i, 'support'] += 1
    df.to_csv(output_csv)


def result_json(distance_csv, output_json):
    df = pd.read_csv(distance_csv)
    not_quasispecies = df.first_record.apply(lambda x: x[:3] != 'qua')
    desired_records = list(set(df.first_record[not_quasispecies]))
    second_is_quasispecies = df.second_record.apply(lambda x: x[:3] == 'qua')
    results = {}
    for record in desired_records:
        if record[:3] != 'qua': continue
        first_is_desired = df.first_record == record
        best_match_index = df.loc[
            first_is_desired & second_is_quasispecies, 'distance'
        ].idxmin()
        results[record] = {
            'best_match': str(df.loc[best_match_index, 'second_record']),
            'distance': int(df.loc[best_match_index, 'distance']),
        }
    with open(output_json, 'w') as json_file:
        json.dump(results, json_file, indent=2)


def covarying_fasta(input_json, input_fasta, output_fasta, end_correction=10):
    with open(input_json) as json_file:
        covarying_sites = json.load(json_file)
    records = list(SeqIO.parse(input_fasta, 'fasta'))
    for record in records:
        last_site = len(record.seq) - end_correction
        record.seq = Seq(
            ''.join([
                record.seq[i]
                for i in covarying_sites
                if i > end_correction and i < last_site
            ])
        )
    SeqIO.write(records, output_fasta, 'fasta')


def report(input_files, output_csv, report_type):
    csvfile = open(output_csv, 'w')
    field_names = ['dataset', 'gene', 'worst_distance', 'report_type']
    writer = csv.DictWriter(csvfile, field_names)
    writer.writeheader()
    for file_path in input_files:
        with open(file_path) as json_file:
            result_data = json.load(json_file)
        dataset = file_path.split('/')[1]
        gene = file_path.split('/')[4]
        worst_distance = 0
        for key, value in result_data.items():
            if value['distance'] > worst_distance:
                worst_distance = value['distance']
        if report_type == 'reconstructing' and worst_distance > 5:
            raise Exception('A reconstruction dataset failed!', dataset)
        writer.writerow({
            'dataset': dataset,
            'gene': gene,
            'worst_distance': worst_distance,
            'report_type': report_type
        })
    csvfile.close()


def haplotyper_report(input_files, output_csv):
    csvfile = open(output_csv, 'w')
    field_names = ['dataset', 'worst_distance']
    writer = csv.DictWriter(csvfile, field_names)
    writer.writeheader()
    for file_path in input_files:
        file_path = file_path.split('.')[0] + '.csv'
        if not os.path.exists(file_path):
            continue
        with open(file_path) as json_file:
            result_data = json.load(json_file)
        for key, value in result_data.items():
            if value['distance'] > worst_distance:
                worst_distance = value['distance']
        writer.writerow({
            'dataset': file_path,
            'worst_distance': worst_distance,
        })
    csvfile.close()


def superread_agreement(input_superreads, input_fasta, input_json, output_csv):
    superreads = list(SeqIO.parse(input_superreads, 'fasta'))
    truth = list(SeqIO.parse(input_fasta, 'fasta'))
    with open(input_json) as json_file:
        sites = np.array(json.load(json_file), dtype=np.int)
    csvfile = open(output_csv, 'w')
    csvwriter = csv.DictWriter(
        csvfile, fieldnames=[
            'superread_id',
            'weight',
            'true_id',
            'smallest_diff',
            'smallest_recomb',
            'start',
            'stop'
        ]
    )
    csvwriter.writeheader()
    n_char = len(sites)
    for superread in superreads:
        smallest_diff = 1e6
        superread_id, weight = superread.name.split('_')
        weight = int(weight.split('-')[1])
        superread_np = np.array(list(superread.seq), dtype='<U1')[sites]
        start = (superread_np != '-').argmax()
        stop = ((np.arange(n_char) >= start) & (superread_np == '-')).argmax()
        smallest_recomb = 1e6
        for true_sequence_a in truth:
            true_a_np = np.array(list(true_sequence_a.seq), dtype='<U1')
            diff = (superread_np[start:stop] != true_a_np[start:stop]).sum()
            if diff < smallest_diff:
                smallest_diff = diff
                smallest_id = true_sequence_a.name
            for true_sequence_b in truth:
                true_b_np = np.array(list(true_sequence_b.seq), dtype='<U1')
                for i in range(start, stop):
                    first = true_a_np[start:i] != superread_np[start:i]
                    second = true_b_np[i:stop] != superread_np[i:stop]
                    recomb = first.sum() + second.sum()
                    if recomb < smallest_recomb:
                        smallest_recomb = recomb
        csvwriter.writerow({
            'superread_id': superread_id,
            'weight': weight,
            'true_id': smallest_id,
            'smallest_diff': smallest_diff,
            'smallest_recomb': smallest_recomb,
            'start': start,
            'stop': stop
        })
    csvfile.close()


def superread_scatter_data(superread_path, output_csv):
    with open(superread_path) as json_file:
        superreads = json.load(json_file)
    pd.DataFrame({
        'weight': [sr['weight'] for sr in superreads],
        'vacs_length': [len(sr['vacs']) for sr in superreads],
    }).to_csv(output_csv)


## HK getProportionAr_5.py
def quantify_ar(input_fasta, superread_path, output_json_path):
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Call files                                                              #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    ## Get all reference sequences
    refNames=[]
    refs=[]
    with open(input_fasta,'r') as myFile:
        for l in myFile:
            if l.startswith('>') and l.startswith('>superread')==False:
                refNames.append(l.rstrip().lstrip('>'))
                refs.append('')
            elif l.startswith('>superread'): break
            else: refs[-1]+=l.rstrip()

    ## Get json file where I can get: superread sequence, indices, and weights
    with open(superread_path,'r') as json_file:
        superread_json=json.load(json_file)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Obtain superreads with clear evidence of AR                             #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    ## Get recombined vs nonrecombined
    ar_total=set() # clear evidence of artificial recombination
    identical_total=set() # entire fragments exactly from ref
    misc_total=set() # either all nt mutated or all from other strain
    superreads=[]
    for ref in refs:
        ar_bin=[]
        identical_bin=[]
        misc_bin=[]
        ## Iterate over superread_json
        for superread_info in superread_json:
            ## Extract variables
            superreads.append(superread_info["index"])
            current_sr=superread_info["index"]
            myRead=superread_info["vacs"]
            first_ind=superread_info["cv_start"]
            last_ind_p1=superread_info["cv_end"]
            refFrag=ref[first_ind:last_ind_p1]
            if myRead==refFrag: identical_bin.append(current_sr)
            else:
                indSubseq_ends=set(['_'.join([str(j) for j in [0,i]]) for i in range(1,len(refFrag)+1)])
                indSubseq_ends=indSubseq_ends|set(['_'.join([str(j) for j in [i,len(refFrag)]]) for i in range(0,len(refFrag)+1)])
                indSubseq_ends=indSubseq_ends-{'_'.join([str(j) for j in [0,len(refFrag)]]),'_'.join([str(j) for j in [len(refFrag),len(refFrag)]])}
                indSubseq_ends=sorted([[int(j) for j in i.split('_')] for i in indSubseq_ends])
                isOverlap=False
                for ends in indSubseq_ends:
                    if refFrag[ends[0]:ends[1]]==myRead[ends[0]:ends[1]]:
                        isOverlap=True
                if isOverlap: ar_bin.append(current_sr)
                else: misc_bin.append(current_sr)
        ar_total=ar_total|set(ar_bin)
        identical_total=identical_total|set(identical_bin)
        misc_total=misc_total|set(misc_bin)
    ## Remove potential set intersects
    ar_total=ar_total-identical_total
    misc_total=misc_total-identical_total-ar_total

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Add weights to AR superreads                                            #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    superreads_count_weights=0
    for item in superread_json:
        superreads_count_weights+=item["weight"]
    ar_weighted=0
    for ar in ar_total:
        ar_weighted+=superread_json[ar]["weight"]

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Write json                                                              #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    to_write = [{ "fasta_name": input_fasta, \
                 "json_name": superread_path, \
                 "quantify_ar": ar_weighted/superreads_count_weights, \
                 "ar_weighted": ar_weighted, \
                 "superreads_tot_weighted": superreads_count_weights, \
                 "superreads_unweighted_length": len(superreads), \
                 "ar_unweighted_length": len(ar_total), \
                 "ar_indices": sorted(list(ar_total))
                 }]
    with open(output_json_path, 'w') as fp:
        json.dump(to_write, fp, indent=2)

## Merging multiple jsons (ex: multiple simulations);
## Input: Use one json with full json path and consolidate with similar ones
def consolidate_simjsons(input_json, output_csv_path):
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Get json files                                                          #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    input_json_split=input_json.split('/')
    diff_path_ind=[ind for ind in range(len(input_json_split)) if input_json_split[ind].startswith('sim')][0]
    dir_jsons=[]
    if input_json.startswith('sim'):
        dir_jsons=os.listdir()
    else:
        dir_jsons=os.listdir('/'.join(input_json_split[:diff_path_ind]))
    dir_jsons=[f for f in dir_jsons if f.startswith(input_json_split[diff_path_ind].split('_')[0])]
    dir_jsons=['/'.join(input_json_split[:diff_path_ind]+[f]+input_json_split[diff_path_ind+1:]) for f in dir_jsons]
    dir_jsons=[f for f in dir_jsons if os.path.isfile(f)]

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Write csv                                                               #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    with open(output_csv_path, 'w') as writeFile:
            writeFile.write(','.join(['path','true_ar','seed','quantified_ar'])+'\n')
            for f in dir_jsons:
                f_info=f.split('/')[diff_path_ind]
                ar=f_info.split('_')[1].split('-')[1]
                s=f_info.split('_')[2].split('-')[1]
                with open(f,'r') as json_file:
                    cur_json=json.load(json_file)
                writeFile.write(','.join([str(i) for i in [f,ar,s,cur_json[0]['quantify_ar']*100]])+'\n') # convert fraction to %
                
def read_statistics(input_bam, output_full, output_summary):
    mapping_qualities = []
    query_lengths = []
    for read in pysam.AlignmentFile(input_bam).fetch(until_eof=True):
        mapping_qualities.append(read.mapping_quality)
        query_lengths.append(read.qlen)
    full = pd.DataFrame({
        'mapping_qualities': mapping_qualities,
        'query_length': query_lengths
    })

    full.to_csv(output_full)
    mq_value_counts = full.mapping_qualities.value_counts().sort_index()
    mq_summary = pd.DataFrame({
        'value': mq_value_counts.index,
        'count': mq_value_counts.values,
        'stat': 'mapping_quality'
    })
    ql_value_counts = full.query_length.value_counts().sort_index()
    ql_summary = pd.DataFrame({
        'value': ql_value_counts.index,
        'count': ql_value_counts.values,
        'stat': 'query_length'
    })
    pd.concat([mq_summary, ql_summary], axis=0).to_csv(output_summary, index=False)


def pluck_superread_reads(input_json, superread_index, input_bam, output_bam):
    with open(input_json) as json_file:
        superreads = json.load(json_file)
    read_names = superreads[superread_index]['read_names']
    query_hash = { query: True for query in read_names }
    af_in = pysam.AlignmentFile(input_bam, 'rb')
    af_out = pysam.AlignmentFile(output_bam, 'wb', template=af_in)
    for read in af_in.fetch():
        if read.query_name in query_hash:
            af_out.write(read)
    af_in.close()
    af_out.close()


def extract_single_genbank(input_genbank):
    genbank_path = os.path.join('input', 'references', input_genbank)
    genbank_record = SeqIO.read(genbank_path, 'genbank')
    virus = input_genbank.split('.')[0]
    is_hiv = virus == 'hiv'
    feature = 'gene' if not is_hiv else 'note'
    genes = [
        feature
        for feature in genbank_record.features
        if feature.type == 'gene' or feature.type == 'CDS'
    ]
    for gene in genes:
        start = gene.location.start.position
        end = gene.location.end.position
        name = gene.qualifiers[feature][0].split()[0]
        seq = gene.extract(genbank_record.seq)
        if len(seq) % 3 == 0:
            filename = 'output/references/%s-%s.fasta' % (virus, name)
            fasta_record = SeqRecord(
                seq,
                id=name,
                description=''
            )
            SeqIO.write(fasta_record, filename, 'fasta')
    SeqIO.write(
        SeqRecord(genbank_record.seq, id=virus, description=''),
        'output/references/%s-wgs.fasta' % virus,
        'fasta'
    )
    return [(virus, gene.qualifiers[feature][0]) for gene in genes]


def extract_genbank_references(reference_csv):
    reference_directory = os.path.join('input', 'references')
    genbank_files = os.listdir(reference_directory)
    reference_list = []
    for genbank_file in genbank_files:
        single_reference_list = extract_single_genbank(genbank_file)
        reference_list.extend(single_reference_list)
    with open(reference_csv, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(('virus', 'gene'))
        for reference in reference_list:
            writer.writerow(reference)


def filter_bam(input_bam_file, output_bam_file, mapping_cutoff=1):
    input_bam = pysam.AlignmentFile(input_bam_file, "rb")
    output_bam = pysam.AlignmentFile(output_bam_file, "wb", template=input_bam)
    for read in input_bam.fetch(until_eof=True):
        if read.mapping_quality >= mapping_cutoff:
            output_bam.write(read)
    input_bam.close()
    output_bam.close()


def write_time(path):
    current_time = int(time())
    with open(path, 'w') as text_file:
        text_file.write('%d' % current_time)


def compute_total_time(start_path, stop_path, total_path):
    with open(start_path) as start_file, open(stop_path) as stop_file:
        start_value = int(start_file.read())
        stop_value = int(stop_file.read())
    total_value = stop_value - start_value
    with open(total_path, 'w') as total_file:
        total_file.write('%d' % total_value)
