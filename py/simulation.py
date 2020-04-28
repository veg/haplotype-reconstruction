import csv
import json
from itertools import tee

import numpy as np
from Bio import SeqIO
import pysam
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')


def extract_lanl_genome(lanl_input, lanl_id, fasta_output):
  records = SeqIO.to_dict(SeqIO.parse(lanl_input, 'fasta'))
  record = records[lanl_id]
  SeqIO.write(record, fasta_output, 'fasta')


def simulate_amplicon_dataset(dataset, gene, output_fastq, output_fasta):
  simulated_reads = []
  true_genes = []
  with open('simulations.json') as json_file:
    simulation_information = json.load(json_file)[dataset]
  for lanl_information in simulation_information:
    lanl_id = lanl_information['lanl_id']
    lanl_reads_filename = "output/lanl/%s/%s/reads.fastq" % (lanl_id, gene)
    lanl_reads = list(SeqIO.parse(lanl_reads_filename, 'fastq'))
    current_frequency = lanl_information['frequency']
    number_of_reads_to_extract = int(current_frequency * len(lanl_reads))
    simulated_reads.extend(lanl_reads[:number_of_reads_to_extract])

    lanl_gene_filename = "output/lanl/%s/%s/sequence.fasta" % (lanl_id, gene)
    lanl_gene = SeqIO.read(lanl_gene_filename, 'fasta')
    true_genes.append(lanl_gene)
  SeqIO.write(simulated_reads, output_fastq, 'fastq')
  SeqIO.write(true_genes, output_fasta, 'fasta')


def create_numeric_fasta(records):
    for_numeric, for_headers = tee(records, 2)
    mr = MappedReads()
    np_arrays = [
        mr.get_numeric_representation(record)
        for record in for_numeric
    ]
    numeric = np.vstack(np_arrays)
    headers = [record.id for record in for_headers]
    return headers, numeric

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


def covarying_sites(input_fasta, output_json):
    fasta_np = np.array([
        list(str(record.seq)) for record in SeqIO.parse(input_fasta, 'fasta')
    ], dtype='<U1')
    A_sum = np.sum(fasta_np == 'A', axis=0)
    C_sum = np.sum(fasta_np == 'C', axis=0)
    G_sum = np.sum(fasta_np == 'G', axis=0)
    T_sum = np.sum(fasta_np == 'T', axis=0)
    max_count = np.max(np.vstack([A_sum, C_sum, G_sum, T_sum]), axis=0)
    covarying_sites = np.arange(fasta_np.shape[1])[max_count < fasta_np.shape[0]]
    with open(output_json, 'w') as json_file:
        json.dump([int(cvs) for cvs in covarying_sites], json_file)


def get_sam_info(sam):
    sam_info = {}
    for i, read in enumerate(sam):
        location = (read.reference_start, read.reference_end)
        if location in sam_info:
            sam_info[location].append(i)
        else:
            sam_info[location] = [i]
    return sam_info


def get_reference_to_alignment_map(lanl_id, aligned_genomes):
    fasta_np = np.array(list(aligned_genomes[lanl_id].seq), dtype='<U1')
    indices = np.arange(len(fasta_np))[fasta_np != '-']
    return indices


def get_alignment_to_reference_map(lanl_id, aligned_genomes):
    fasta_np = np.array(list(aligned_genomes[lanl_id].seq), dtype='<U1')
    indices = np.cumsum(fasta_np != '-') - 1
    return indices


def get_mate(
        read, left_strain, right_strain, sams, sam_infos,
        r2a_maps, a2r_maps, stop=25
        ):
    sam = sams[right_strain]
    sam_info = sam_infos[right_strain]
    i = 0
    left_alignment_start = r2a_maps[left_strain][read.reference_start]
    left_alignment_end = r2a_maps[left_strain][read.reference_end-1]
    while True:
        for j in range(i+1):
            all_shifts = [(j, i-j), (j, j-i), (-j, i-j), (-j, j-i)]
            for left_shift, right_shift in all_shifts:
                right_alignment_start = left_alignment_start + left_shift
                right_alignment_end = left_alignment_end + right_shift
 
                starts_before = right_alignment_start < 0
                ends_after = right_alignment_end >= len(a2r_maps[right_strain])

                if starts_before or ends_after:
                    continue

                right_reference_start = a2r_maps[right_strain][right_alignment_start]
                right_reference_end = a2r_maps[right_strain][right_alignment_end]
                location = (right_reference_start, right_reference_end)
                if location in sam_info:
                    n = len(sam_info[location])
                    i = np.random.randint(n)
                    return sam_info[location][i]
        i += 1
        if i == stop:
            return None


def write_ar_dataset(
        lanl_ids, frequencies, ar, aligned_genomes, output_fastq, output_sam,
        number_of_reads
        ):
    should_recombine = len(lanl_ids) > 1
    sams = [
        list(pysam.AlignmentFile('output/lanl/%s/wgs.sam' % lanl_id, "r"))
        for lanl_id in lanl_ids
    ]
    sam_infos = [get_sam_info(sam) for sam in sams]
    number_of_ar_reads = np.ceil(ar*number_of_reads).astype(np.int)
    number_of_clean_reads = number_of_reads - number_of_ar_reads
    read_counts = (number_of_clean_reads, number_of_ar_reads)
    print('Simulating %d pure reads, %d recombined...' % read_counts)
    number_of_strains = len(frequencies)
    ar_left_strains = np.random.choice(
        number_of_strains, 2*number_of_ar_reads, p=frequencies
    )
    ar_left_indices = np.random.choice(
        number_of_clean_reads, 2*number_of_ar_reads, replace=False
    )
    ar_right_strains = np.zeros(2*number_of_ar_reads, dtype=np.int)
    if should_recombine:
        for i in range(number_of_strains):
            other_strains = [j for j in range(number_of_strains) if j != i]
            new_frequencies = np.array([frequencies[j] for j in other_strains])
            new_frequencies = new_frequencies/np.sum(new_frequencies)
            is_current_strain = ar_left_strains == i

            number_of_current_strain = is_current_strain.sum()
            ar_right_strains[ar_left_strains == i] = np.random.choice(
                other_strains, number_of_current_strain, p=new_frequencies
            )

        reference_to_alignment_maps = [
            get_reference_to_alignment_map(lanl_id, aligned_genomes)
            for lanl_id in lanl_ids
        ]

        alignment_to_reference_maps = [
            get_alignment_to_reference_map(lanl_id, aligned_genomes)
            for lanl_id in lanl_ids
        ]
    total_unsuitable = 0
    output_file = open(output_fastq, 'w')
    if should_recombine:
        i = 0
        for _ in range(number_of_ar_reads):
            found = False
            while not found:
                left_strain = ar_left_strains[i]
                left_read_index = ar_left_indices[i]
                left_read = sams[left_strain][left_read_index]
                right_strain = ar_right_strains[i]
                right_read_index = get_mate(
                    left_read, left_strain, right_strain, sams, sam_infos,
                    reference_to_alignment_maps, alignment_to_reference_maps
                )
                i += 1
                if right_read_index is not None:
                    found = True
                else:
                    total_unsuitable += 1
            right_read = sams[right_strain][right_read_index]

            left_aligned_pairs = left_read.get_aligned_pairs(matches_only=True)
            right_aligned_pairs = right_read.get_aligned_pairs(matches_only=True)
            left_a2r_map = alignment_to_reference_maps[left_strain]
            left_r2a_map = reference_to_alignment_maps[left_strain]
            right_a2r_map = alignment_to_reference_maps[right_strain]
            right_r2a_map = reference_to_alignment_maps[right_strain]

            recombination_lower = np.max([
                left_r2a_map[left_aligned_pairs[0][1]],
                right_r2a_map[right_aligned_pairs[0][1]]
            ])
            recombination_upper = np.min([
                left_r2a_map[left_aligned_pairs[-1][1]],
                right_r2a_map[right_aligned_pairs[-1][1]]
            ])
            recombination_site = np.random.randint(
                recombination_lower, recombination_upper
            )

            header = '@' + left_read.query_name + '+' + right_read.query_name + '\n'
            output_file.write(header)
            left_query = ''.join([
                left_read.query[pair[0]]
                for pair in left_aligned_pairs
                if pair[1] < left_a2r_map[recombination_site]
            ])
            right_query = ''.join([
                right_read.query[pair[0]]
                for pair in right_aligned_pairs
                if pair[1] >= right_a2r_map[recombination_site]
            ])
            query = left_query + right_query
            output_file.write(query + '\n')
            output_file.write('+\n')
            left_quality = ''.join([
                left_read.qual[pair[0]]
                for pair in left_aligned_pairs
                if pair[1] < left_a2r_map[recombination_site]
            ])
            right_quality = ''.join([
                right_read.qual[pair[0]]
                for pair in right_aligned_pairs
                if pair[1] >= right_a2r_map[recombination_site]
            ])
            quality = left_quality + right_quality
            output_file.write(quality + '\n')

    nonrecombined_strains = np.random.choice(
        number_of_strains, number_of_clean_reads, p=frequencies
    )
    nonrecombined_indices = np.random.choice(
        number_of_reads, number_of_clean_reads, replace=False
    )
    for_template = pysam.AlignmentFile('output/lanl/%s/wgs.sam' % lanl_ids[0], "r")
    sam_file = pysam.AlignmentFile(output_sam, 'w', template=for_template)
    for_template.close()
    for strain, index in zip(nonrecombined_strains, nonrecombined_indices):
        read = sams[strain][index]
        sam_file.write(read)
        output_file.write('@' + read.query_name + '\n')
        output_file.write(read.query + '\n')
        output_file.write('+\n')
        output_file.write(read.qual + '\n')
    print('Total unsuitable reads:', total_unsuitable)
    output_file.close()
    sam_file.close()

def simulation_truth(dataset, output_fasta):
    with open('simulations.json') as json_file:
        simulation_information = json.load(json_file)[dataset]
    lanl_ids = [lanl_info['lanl_id'] for lanl_info in simulation_information]
    true_genomes = []
    for lanl_id in lanl_ids:
        lanl_genome_filename = "output/lanl/%s/genome.fasta" % lanl_id
        lanl_genome = SeqIO.read(lanl_genome_filename, 'fasta')
        true_genomes.append(lanl_genome)
    SeqIO.write(true_genomes, output_fasta, 'fasta')


def simulate_wgs_dataset(
        dataset, ar, input_fasta, output_fastq, output_json, output_sam,
        seed=1, number_of_reads=300000
        ):
    np.random.seed(int(seed))
    ar = float(ar)/100
    with open('simulations.json') as json_file:
        simulation_information = json.load(json_file)[dataset]
    lanl_ids = [lanl_info['lanl_id'] for lanl_info in simulation_information]
    frequencies = np.array(
        [lanl_info['frequency'] for lanl_info in simulation_information]
    )
    aligned_genomes = SeqIO.to_dict(
        SeqIO.parse(input_fasta, 'fasta')
    )
    write_ar_dataset(
        lanl_ids, frequencies, ar, aligned_genomes, output_fastq, output_sam,
        number_of_reads
        )

    with open(output_json, 'w') as json_file:
        json.dump(
            evaluate_simulated_ar(lanl_ids, output_fastq), json_file, indent=2
        )


def evaluate_simulated_ar(lanl_ids, filename):
    ar_dataset = list(SeqIO.parse(filename, 'fastq'))
    recombined_reads = 0
    strain_counts = [0 for i in lanl_ids]
    total_reads = len(ar_dataset)
    info = {}
    for read in ar_dataset:
        if '+' in read.id:
            recombined_reads += 1
        else:
            for i, lanl_id in enumerate(lanl_ids):
                if lanl_id in read.id:
                    strain_counts[i] += 1
    non_recombined_reads = sum(strain_counts)
    info['totalReads'] = total_reads
    info['recombination'] = recombined_reads / total_reads
    for i, lanl_id in enumerate(lanl_ids):
        info[lanl_id] = strain_counts[i]/non_recombined_reads
    return info


def n_paths_boxplot(simulated_dataset, gene, output_filepath):
    template_string = "output/sim-%s_ar-%d_seed-%d/fastp/bowtie2/%s/acme/graph.json"
    input_files = []
    data_seeds = []
    data_ars = []
    data_npaths = []
    for seed in range(1, 11):
        for ar in [0, 5, 10, 15, 20]:
            json_filepath = template_string % (simulated_dataset, ar, seed, gene)
            with open(json_filepath) as json_file:
                graph = json.load(json_file)
            data_seeds.append(seed)
            data_ars.append(ar)
            data_npaths.append(graph['number_of_paths'])
    df = pd.DataFrame({
        "seed": data_seeds,
        "ARRate": data_ars,
        "LogNumberOfPaths": data_npaths,
        "LogNumberOfPaths": np.log10(data_npaths),
    })
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.set(font_scale=5)
    sns.boxplot(x='ARRate', y='LogNumberOfPaths', data=df, ax=ax)
    fig.savefig(output_filepath)


def superread_weight_distribution_data(superreads_filepath, csv_filepath):
    with open(superreads_filepath) as json_file:
        superreads = json.load(json_file)
    csv_file = open(csv_filepath, 'w')
    csv_writer = csv.DictWriter(csv_file, fieldnames=['weight', 'composition'])
    csv_writer.writeheader()
    for superread in superreads:
        composition = max(
            superread['composition'].items(),
            key=lambda i: i[1]
        )[0]
        csv_writer.writerow({
            'weight': superread['weight'],
            'composition': composition
        })
    csv_file.close()


def simulation_coverage(input_sam, output_csv, output_png):
    mapped_reads = pysam.AlignmentFile(input_sam)
    n_sites = mapped_reads.lengths[0]
    all_coverage = {}
    for read in mapped_reads.fetch():
        key = read.query_name.split('-')[0]
        coverage = np.zeros(n_sites)
        coverage[read.reference_start: read.reference_end] += 1
        if not key in all_coverage:
            all_coverage[key] = np.zeros(n_sites)
        all_coverage[key][read.reference_start: read.reference_end] += 1
    df = pd.concat([
        pd.DataFrame({
            'strain': key,
            'site': np.arange(n_sites),
            'coverage': value
        })
        for key, value in all_coverage.items()
    ], axis=0)
    df.to_csv(output_csv)
    fig, ax = plt.subplots(figsize=(30, 10))
    for strain in set(df.strain):
        ax.plot(
            df.loc[df.strain == strain, 'site'],
            df.loc[df.strain == strain, 'coverage']
        )
    fig.savefig(output_png)
