import os

from py import embed_and_reduce_dimensions_io
from py import cluster_blocks_io
from py import obtain_consensus_io


ACCESSION_NUMBERS = ['ERS6610%d' % i for i in range(87, 94)]
SIMULATED_DATASETS = ["related_1", "related_2", "diverged_1", "diverged_2"]
REFERENCES = ["env_C2V5", "gag_p24", "pr", "rt"]
HYPHY_PATH = "/Users/stephenshank/Software/lib/hyphy"

##################
# RECONSTRUCTION #
##################

# Simulation

rule extract_lanl:
  input:
    "input/LANL-HIV.fasta"
  output:
    "output/related_1/genome.fasta",
    "output/related_2/genome.fasta",
    "output/diverged_1/genome.fasta",
    "output/diverged_2/genome.fasta"
  shell:
    "python py/extract_sequences.py"

rule extract_gene:
  input:
    reference="input/references/{gene}.fasta",
    genome="output/{simulated_dataset}/genome.fasta"
  output:
    sam="output/{simulated_dataset}/{gene}/sequence.sam",
    fasta="output/{simulated_dataset}/{gene}/sequence.fasta"
  shell:
    """
      bealign -r {input.reference} {input.genome} {output.sam}
      bam2msa {output.sam} {output.fasta}
    """

rule align_simulated:
  input:
    reference="input/references/{gene}.fasta",
    target=rules.extract_gene.output.fasta
  output:
    unaligned="output/{simulated_dataset}/{gene}/unaligned.fasta",
    aligned="output/{simulated_dataset}/{gene}/aligned.fasta"
  shell:
    """
      cat {input.target} {input.reference} > {output.unaligned}
      mafft --localpair {output.unaligned} > {output.aligned}
    """

rule amplicon_simulate_single:
  input:
    "output/{simulated_dataset}/{gene}/sequence.fasta"
  output:
    "output/{simulated_dataset}/{gene}/reads.fastq",
  shell:
    """
      art_illumina -ss HS25 -i {input} -l 120 -s 50 -c 15000 -o {output}
      mv {output}.fq {output}
    """

rule amplicon_simulate_mixed:
  input:
    "output/{mixed_dataset}_1/{gene}/reads.fastq",
    "output/{mixed_dataset}_2/{gene}/reads.fastq"
  output:
    "output/{mixed_dataset}_joint/{gene}/reads.fastq",
  shell:
    "cat {input[0]} {input[1]} > {output}"

rule wgs_simulate_single:
  input:
    "output/{simulated_dataset}/genome.fasta"
  output:
    "output/{simulated_dataset}/reads.fastq"
  shell:
    """
      art_illumina -ss HS25 -i {input} -l 120 -s 50 -c 15000 -o {output}
      mv {output}.fq {output}
    """

rule wgs_simulate_mixed:
  input:
    "output/{mixed_dataset}_1/reads.fastq",
    "output/{mixed_dataset}_2/reads.fastq"
  output:
    "output/{mixed_dataset}_joint/reads.fastq"
  shell:
    "cat {input[0]} {input[1]} > {output}"

# Situating other data

rule situate_intrahost_data:
  input:
    "input/evolution/{accession_number}.fastq"
  output:
    "output/{accession_number}/reads.fastq"
  shell:
    "cp {input} {output}"

# Quality control

rule qfilt:
  input:
    "output/{dataset}/reads.fastq"
  output:
    fasta="output/{dataset}/qfilt/qc.fasta",
    json="output/{dataset}/qfilt/qc.json",
  shell:
    "qfilt -Q {input} -q 20 -l 50 -P - -R 8 -j >> {output.fasta} 2>> {output.json}"


rule fastp:
  input:
    "output/{dataset}/reads.fastq"
  output:
    fastq="output/{dataset}/fastp/qc.fastq",
    json="output/{dataset}/fastp/qc.json",
    html="output/{dataset}/fastp/qc.html"
  shell:
    "fastp -A -q 10 -i {input} -o {output.fastq} -j {output.json} -h {output.html}"

rule trimmomatic:
  input:
    "output/{dataset}/reads.fastq"
  output:
    "output/{dataset}/trimmomatic/qc.fastq"
  shell:
    "trimmomatic SE {input} {output} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

# Read mapping

rule bealign:
  input:
    qc="output/{dataset}/qfilt/qc.fasta",
    reference="input/references/{gene}.fasta"
  output:
    bam="output/{dataset}/qfilt/bealign/{gene}/mapped.bam",
    discards="output/{dataset}/qfilt/bealign/{gene}/discards.fasta"
  shell:
    "bealign -r {input.reference} -e 0.5 -m HIV_BETWEEN_F -D {output.discards} -R {input.qc} {output.bam}"

rule bwa_reference_index:
  input:
    "input/references/{reference}.fasta"
  output:
    "output/references/{reference}.fasta"
  shell:
    """
      cp {input} {output}
      bwa index {output}
    """

rule bwa_map_reads:
  input:
    fastq="output/{dataset}/{qc}/qc.fastq",
    reference="output/references/{reference}.fasta"
  output:
    "output/{dataset}/{qc}/bwa/{reference}/mapped.bam"
  shell:
    "bwa mem {input.reference} {input.fastq} > {output}"

rule sort_and_index:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/mapped.bam"
  output:
    bam="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.fasta",
    index="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam.bai"
  shell:
    """
      samtools sort {input} > {output.bam}
      bam2msa {output.bam} {output.fasta}
      samtools index {output.bam}
    """

# Haplotype reconstruction (full pipelines)

rule regress_haplo_full:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam",
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam.bai",
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/final_haplo.fasta"
  script:
    "R/regress_haplo/full_pipeline.R"

# ACME haplotype reconstruction

rule mmvc:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.fasta"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/mmvc.json"
  shell:
    "mmvc -j {output} {input}"

rule embed_and_reduce_dimensions:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.fasta"
  output:
    csv="output/{dataset}/{qc}/{read_mapper}/{reference}/dr_{dim}d.csv",
    json="output/{dataset}/{qc}/{read_mapper}/{reference}/dr_{dim}d.json"
  params:
    dim=lambda w: int(w.dim)
  run:
    embed_and_reduce_dimensions_io(input[0], output.csv, output.json, params.dim)

rule cluster_blocks:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/dr_{dim}d.csv"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/clustered_{dim}d.csv"
  params:
    dim=lambda w: int(w.dim)
  run:
    cluster_blocks_io(input[0], output[0], params.dim)

rule obtain_consensus:
  input:
    csv="output/{dataset}/{qc}/{read_mapper}/{reference}/clustered_{dim}d.csv",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.fasta",
    json="output/{dataset}/{qc}/{read_mapper}/{reference}/dr_{dim}d.json"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/contigs_{dim}d.fasta",
  run:
    obtain_consensus_io(input.csv, input.fasta, input.json, output[0])

## Regress Haplo
#
#rule regress_haplo_bam_to_variant_calls:
#  input:
#    "output/{dataset}/{reference}/sorted.bam",
#    "output/{dataset}/{reference}/sorted.bam.bai"
#  output:
#    "output/{dataset}/{reference}/variant_calls.csv"
#  script:
#    "R/regress_haplo/bam_to_variant_calls.R"
#   
#rule regress_haplo_variant_calls_to_read_table:
#  input:
#    "output/{dataset}/{reference}/sorted.bam",
#    "output/{dataset}/{reference}/variant_calls.csv",
#  output:
#    "output/{dataset}/{reference}/read_table.csv"
#  script:
#    "R/regress_haplo/variant_calls_to_read_table.R"
#
#rule regress_haplo_read_table_to_loci:
#  input:
#    rules.regress_haplo_variant_calls_to_read_table.output[0]
#  output:
#    "output/{dataset}/{reference}/loci.csv"
#  script:
#    "R/regress_haplo/read_table_to_loci.R"
#
#rule regress_haplo_loci_to_haplotypes:
#  input:
#    rules.regress_haplo_read_table_to_loci.output[0]
#  output:
#    "output/{dataset}/{reference}/h.csv"
#  script:
#    "R/regress_haplo/loci_to_haplotypes.R"
#
#rule regress_haplo_haplotypes_to_parameters:
#  input:
#    rules.regress_haplo_loci_to_haplotypes.output[0]
#  output:
#    "output/{dataset}/{reference}/P.csv"
#  script:
#    "R/regress_haplo/haplotypes_to_parameters.R"
#
#rule regress_haplo_parameters_to_solutions:
#  input:
#    rules.regress_haplo_haplotypes_to_parameters.output[0]
#  output:
#    "output/{dataset}/{reference}/solutions.csv"
#  script:
#    "R/regress_haplo/parameters_to_solutions.R"
#
#rule regress_haplo_solutions_to_haplotypes:
#  input:
#    rules.regress_haplo_parameters_to_solutions.output[0]
#  output:
#    "output/{dataset}/{reference}/final_haplo.csv"
#  script:
#    "R/regress_haplo/solutions_to_haplotypes.R"
#
#rule regress_haplo_haplotypes_to_fasta:
#  input:
#    rules.regress_haplo_bam_to_variant_calls.input[0],
#    rules.regress_haplo_solutions_to_haplotypes.output[0]
#  output:
#    "output/{dataset}/{reference}/final_haplo.fasta"
#  script:
#    "R/regress_haplo/haplotypes_to_fasta.R"
#
#
##############
## EVOLUTION #
##############
#
#rule concatenate:
#  input:
#    expand("output/{dataset}/{{reference}}/full/final_haplo.fasta", dataset=ACCESSION_NUMBERS)
#  output:
#    "output/{reference}/unaligned.fasta"
#  params:
#    lambda wildcards: ' '.join(["output/%s/%s/haplotypes/final_haplo.fasta" % (accession, wildcards.reference) for accession in ACCESSION_NUMBERS])
#  shell:
#    "cat {params} > {output}"
#
#rule alignment:
#  input:
#    rules.concatenate.output[0]
#  output:
#    "output/{reference}/aligned.fasta"
#  shell:
#    "mafft {input} > {output}"
#
#rule recombination_screening:
#  input:
#    rules.alignment.output[0]
#  output:
#    gard_json="output/{reference}/GARD.json",
#    nexus="output/{reference}/seqs_and_trees.nex"
#  params:
#    gard_path="%s/TemplateBatchFiles/GARD.bf" % HYPHY_PATH,
#    gard_output=os.getcwd() + "/output/{reference}/aligned.GARD",
#    final_out=os.getcwd() + "/output/{reference}/aligned.GARD_finalout",
#    translate_gard_j=os.getcwd() + "/output/{reference}/aligned.GARD.json",
#    translated_json=os.getcwd() + "/output/{reference}/GARD.json",
#    lib_path=HYPHY_PATH,
#    alignment_path=os.getcwd() + "/output/{reference}/aligned.fasta"
#  shell:
#    """
#      mpirun -np 2 HYPHYMPI LIBPATH={params.lib_path} {params.gard_path} {params.alignment_path} '010010' None {params.gard_output}
#      translate-gard -i {params.gard_output} -j {params.translate_gard_j} -o {params.translated_json}
#      mv {params.final_out} {output.nexus}
#    """
#
#rule site_selection:
#  input:
#    rules.recombination_screening.output.nexus
#  output:
#    "output/{reference}/seqs_and_trees.nex.FUBAR.json"
#  params:
#    full_nexus_path=os.getcwd() + "/" + rules.recombination_screening.output.nexus,
#    fubar_path="%s/TemplateBatchFiles/SelectionAnalyses/FUBAR.bf" % HYPHY_PATH,
#    lib_path=HYPHY_PATH
#  shell:
#    "(echo 1; echo {params.full_nexus_path}; echo 20; echo 1; echo 5; echo 2000000; echo 1000000; echo 100; echo .5;) | HYPHYMP LIBPATH={params.lib_path} {params.fubar_path}"
#
#rule gene_selection:
#  input:
#    rules.recombination_screening.output.nexus
#  output:
#    "output/{reference}/seqs_and_trees.nex.BUSTED.json"
#  params:
#    full_nexus_path=os.getcwd() + "/" + rules.recombination_screening.output.nexus,
#    busted_path="%s/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf" % HYPHY_PATH,
#    lib_path=HYPHY_PATH
#  shell:
#    "(echo 1; echo {params.full_nexus_path}; echo 2;) | HYPHYMP LIBPATH={params.lib_path} {params.busted_path}"
#
#rule full_analysis:
#  input:
#    rules.site_selection.output[0],
#    rules.gene_selection.output[0]
#  output:
#    "output/{reference}/results.tar.gz"
#  shell:
#    "tar cvzf {output} {input[0]} {input[1]}"
#
