import json

import numpy as np
from Bio import SeqIO


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


def simulate_wgs_dataset(dataset, output_fastq, output_fasta):
  simulated_reads = []
  true_genomes = []
  with open('simulations.json') as json_file:
    simulation_information = json.load(json_file)[dataset]
  for lanl_information in simulation_information:
    lanl_id = lanl_information['lanl_id']
    lanl_reads_filename = "output/lanl/%s/wgs.fastq" % lanl_id
    lanl_reads = list(SeqIO.parse(lanl_reads_filename, 'fastq'))
    current_frequency = lanl_information['frequency']
    number_of_reads_to_extract = int(current_frequency * len(lanl_reads))
    simulated_reads.extend(lanl_reads[:number_of_reads_to_extract])

    lanl_genome_filename = "output/lanl/%s/genome.fasta" % lanl_id
    lanl_genome = SeqIO.read(lanl_genome_filename, 'fasta')
    true_genomes.append(lanl_genome)
  SeqIO.write(simulated_reads, output_fastq, 'fastq')
  SeqIO.write(true_genomes, output_fasta, 'fasta')


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

