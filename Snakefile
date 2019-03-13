import os
import json

from py import extract_lanl_genome
from py import simulate_amplicon_dataset 
from py import simulate_wgs_dataset 
from py import write_abayesqr_config
from py import embed_and_reduce_dimensions_io
from py import cluster_blocks_io
from py import obtain_consensus_io
from py import superreads_to_haplotypes_io
from py import evaluate
from py import parse_abayesqr_output


with open('simulations.json') as simulation_file:
  SIMULATION_INFORMATION = json.load(simulation_file)
ACCESSION_NUMBERS = ['ERS6610%d' % i for i in range(87, 94)]
SIMULATED_DATASETS = ['simulation_' + dataset for dataset in SIMULATION_INFORMATION.keys()]
RECONSTRUCTION_DATASETS = [
  "93US141_100k_14-159320-1GN-0_S16_L001_R1_001",
  "PP1L_S45_L001_R1_001",
  "sergei1"
]
ALL_DATASETS = ACCESSION_NUMBERS + SIMULATED_DATASETS + RECONSTRUCTION_DATASETS
REFERENCES = ["env_C2V5", "gag_p24", "pr", "rt"]
HYPHY_PATH = "/Users/stephenshank/Software/lib/hyphy"

##################
# RECONSTRUCTION #
##################

# Simulation

rule extract_lanl_genome:
  input:
    "input/LANL-HIV.fasta"
  output:
    "output/lanl/{lanl_id}/genome.fasta"
  run:
    extract_lanl_genome(input[0], wildcards.lanl_id, output[0])

rule extract_gene:
  input:
    reference="input/references/{gene}.fasta",
    genome="output/lanl/{lanl_id}/genome.fasta"
  output:
    sam="output/lanl/{lanl_id}/{gene}/sequence.sam",
    fasta="output/lanl/{lanl_id}/{gene}/sequence.fasta"
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

rule amplicon_simulation:
  input:
    "output/lanl/{lanl_id}/{gene}/sequence.fasta"
  output:
    "output/lanl/{lanl_id}/{gene}/reads.fastq"
  shell:
    """
      art_illumina -ss HS25 -i {input} -l 120 -s 50 -c 15000 -o {output}
      mv {output}.fq {output}
    """

def amplicon_simulation_inputs(wildcards):
  dataset = SIMULATION_INFORMATION[wildcards.simulated_dataset]
  lanl_ids = [info['lanl_id'] for info in dataset]
  reads = ["output/lanl/%s/%s/reads.fastq" % (lanl_id, wildcards.gene) for lanl_id in lanl_ids]
  genes = ["output/lanl/%s/%s/sequence.fasta" % (lanl_id, wildcards.gene) for lanl_id in lanl_ids]
  return reads + genes

rule simulate_amplicon_dataset:
  input:
    amplicon_simulation_inputs
  output:
    fastq="output/simulation_{simulated_dataset}/{gene}/reads.fastq",
    fasta="output/simulation_{simulated_dataset}/{gene}/truth.fasta"
  run:
    simulate_amplicon_dataset(wildcards.simulated_dataset, wildcards.gene, output.fastq, output.fasta)

rule wgs_simulation:
  input:
    "output/lanl/{lanl_id}/genome.fasta"
  output:
    "output/lanl/{lanl_id}/wgs.fastq"
  shell:
    """
      art_illumina -ss HS25 -i {input} -l 120 -s 50 -c 150000 -o {output}
      mv {output}.fq {output}
    """

def wgs_simulation_inputs(wildcards):
  dataset = SIMULATION_INFORMATION[wildcards.simulated_dataset]
  lanl_ids = [info['lanl_id'] for info in dataset]
  reads = ["output/lanl/%s/wgs.fastq" % (lanl_id) for lanl_id in lanl_ids]
  genomes = ["output/lanl/%s/genome.fasta" % (lanl_id) for lanl_id in lanl_ids]
  return reads + genomes

rule simulate_wgs_dataset:
  input:
    wgs_simulation_inputs
  output:
    fastq="output/simulation_{simulated_dataset}/wgs.fastq",
    fasta="output/simulation_{simulated_dataset}/truth.fasta"
  run:
    simulate_wgs_dataset(wildcards.simulated_dataset, output.fastq, output.fasta)

# Situating other data

def situate_input(wildcards):
  dataset = wildcards.dataset
  is_evolution_dataset = dataset[:7] == 'ERS6610'
  if is_evolution_dataset:
    return "input/evolution/%s.fastq" % dataset
  is_simulated_dataset = 'simulation' in dataset
  if is_simulated_dataset:
    return "output/%s/wgs.fastq" % dataset
  return "input/reconstruction/%s.fastq" % dataset

rule situate_intrahost_data:
  input:
    situate_input
  output:
    "output/{dataset}/reads.fastq"
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
    "qfilt -Q {input} -q 15 -l 50 -P - -R 8 -j >> {output.fasta} 2>> {output.json}"

rule qfilt_454:
  input:
    fasta="input/reconstruction/3.GAC.454Reads.fna",
    qual="input/reconstruction/3.GAC.454Reads.qual"
  output:
    fasta="output/454/qfilt/qc.fasta",
    json="output/454/qfilt/qc.json"
  shell:
    "qfilt -F {input.fasta} {input.qual} -q 15 -l 50 -P - -R 8 -j >> {output.fasta} 2>> {output.json}"

rule fastp:
  input:
    "output/{dataset}/reads.fastq"
  output:
    fastq="output/{dataset}/fastp/qc.fastq",
    json="output/{dataset}/fastp/qc.json",
    html="output/{dataset}/fastp/qc.html"
  shell:
    "fastp -A -q 8 -i {input} -o {output.fastq} -j {output.json} -h {output.html}"

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
    sam="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.sam",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.fasta",
    index="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam.bai"
  shell:
    """
      samtools sort {input} > {output.bam}
      samtools view -h {output.bam} > {output.sam}
      bam2msa {output.bam} {output.fasta}
      samtools index {output.bam}
    """

rule all_acme_bams:
  input:
    expand(
      "output/{dataset}/qfilt/bealign/{reference}/sorted.bam",
      dataset=ALL_DATASETS,
      reference=REFERENCES
    )

# Haplotype reconstruction (full pipelines)

rule regress_haplo_full:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam",
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam.bai",
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/final_haplo.fasta"
  script:
    "R/regress_haplo/full_pipeline.R"

rule all_acme_regress_haplo:
  input:
    expand(
      "output/{dataset}/qfilt/bealign/{reference}/final_haplo.fasta",
      dataset=ALL_DATASETS,
      reference=REFERENCES
    )

rule haploclique:
  input:
    bam="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam",
    index="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam.bai"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/haplo_clique/result.fasta"
  params:
    move="output/{dataset}/{qc}/{read_mapper}/{reference}/haplo_clique/result.fasta.fasta"
  shell:
    """
      haploclique {input.bam} {output}
      mv {params.move} {output}
    """

rule quasirecomb:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/quasirecomb/output.fasta"
  params:
    basedir="output/{dataset}/{qc}/{read_mapper}/{reference}/quasirecomb/"
  shell:
    """
      java -jar ~/QuasiRecomb.jar -i {input}
      mv support/* {params.basedir}
      rmdir support
      mv piDist.txt {params.basedir}
      mv quasispecies.fasta {params.basedir}/output.fasta
    """

rule abayesqr:
  input:
    sam="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.sam",
    reference="input/references/{reference}.fasta"
  output:
    config="output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/config",
    freq="output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/test_Freq.txt",
    seq="output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/test_Seq.txt",
    viralseq="output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/test_ViralSeq.txt",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/haplotypes.fasta"
  run:
    write_abayesqr_config(input.sam, input.reference, output.config)
    shell("aBayesQR {output.config}")
    shell("mv test_Freq.txt {output.freq}")
    shell("mv test_Seq.txt {output.seq}")
    shell("mv test_ViralSeq.txt {output.viralseq}")
    parse_abayesqr_output(output.viralseq, output.fasta)

rule all_acme_abayesqr:
  input:
    expand(
      "output/{dataset}/qfilt/bealign/{reference}/abayesqr/test_Seq.txt",
      dataset=ALL_DATASETS,
      reference=REFERENCES
    )

# ACME haplotype reconstruction

rule mmvc:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.fasta"
  output:
    json="output/{dataset}/{qc}/{read_mapper}/{reference}/mmvc.json",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/mmvc.fasta"
  shell:
    "mmvc -j {output.json} -f {output.fasta} {input}"

rule embed_and_reduce_dimensions:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/mmvc.fasta"
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

rule cluster_blocks_image:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/clustered_2d.csv"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/clustered_2d.png"
  script:
    "R/cluster_plot.R"

rule obtain_consensus:
  input:
    csv="output/{dataset}/{qc}/{read_mapper}/{reference}/clustered_{dim}d.csv",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.fasta",
    json="output/{dataset}/{qc}/{read_mapper}/{reference}/dr_{dim}d.json"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/superreads_{dim}d.fasta",
  run:
    obtain_consensus_io(input.csv, input.fasta, input.json, output[0])

rule superreads_and_truth:
  input:
    superreads="output/{dataset}/{qc}/{read_mapper}/{reference}/superreads_{dim}d.fasta",
    truth="output/{dataset}/{reference}/truth.fasta"
  output:
    unaligned="output/{dataset}/{qc}/{read_mapper}/{reference}/truth_and_superreads_{dim}d_unaligned.fasta",
    aligned="output/{dataset}/{qc}/{read_mapper}/{reference}/truth_and_superreads_{dim}d.fasta"
  shell:
    """
      cat {input.truth} {input.superreads} > {output.unaligned}
      mafft {output.unaligned} > {output.aligned}
    """

rule readreduce:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/mmvc.fasta"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/haplosuperreads.fasta"
  shell:
    "readreduce -a resolve -l 30 -s 16 -o {output} {input}"

rule superreads_to_haplotypes:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/superreads_{dim}d.fasta"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/haplotypes_{dim}d.fasta"
  run:
    superreads_to_haplotypes_io(input[0], output[0])

rule haplotypes_and_truth:
  input:
    haplotypes="output/{dataset}/{qc}/{read_mapper}/{reference}/haplotypes_{dim}d.fasta",
    truth="output/{dataset}/{reference}/truth.fasta"
  output:
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/truth_and_haplotypes_{dim}d.fasta",
    json="output/{dataset}/{qc}/{read_mapper}/{reference}/truth_and_haplotypes_{dim}d.json"
  run:
    shell("cat {input.truth} {input.haplotypes} > {output.fasta}")
    evaluate(input.haplotypes, input.truth, output.json)

rule dashboard:
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/dashboard.html",
  params:
    path="output/{dataset}/{qc}/{read_mapper}/{reference}",
  shell:
    "npx webpack --output-path {params.path}"

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
