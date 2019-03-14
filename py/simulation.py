import json

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

