import os
import json

import pysam
from convex_qsr import *
from py import *


with open('simulations.json') as simulation_file:
  SIMULATION_INFORMATION = json.load(simulation_file)
ACCESSION_NUMBERS = ['ERS6610%d' % i for i in range(87, 94)]
SIMULATED_DATASETS = ['sim-' + dataset for dataset in SIMULATION_INFORMATION.keys()]
RECONSTRUCTION_DATASETS = [
  "93US141_100k_14-159320-1GN-0_S16_L001_R1_001",
  "PP1L_S45_L001_R1_001",
  "sergei1",
  "FiveVirusMixIllumina_1"
]
FVM_RECORDS = ['89_6', 'HXB2', 'JRCSF', 'NL43', 'YU2']
FVM_STRING = 'FiveVirusMix'
ALL_DATASETS = ACCESSION_NUMBERS + SIMULATED_DATASETS + RECONSTRUCTION_DATASETS
KNOWN_TRUTH = SIMULATED_DATASETS + ["FiveVirusMixIllumina_1"]
ALL_REFERENCES = ["env", "rev", "vif", "pol", "prrt", "rt", "pr", "gag", "int", "tat"] # + ["nef", "vpr"] 
REFERENCE_SUBSET = ["env", "pol", "gag"]
HYPHY_PATH = "/Users/stephenshank/Software/lib/hyphy"
HAPLOTYPERS = ["abayesqr", "regress_haplo", "quasirecomb", "acme"]
NUMBER_OF_READS = 300000
with open('input/sarscov2_accessions.txt') as sarscov2_accession_file:
  SARSCOV2_ACCESSIONS = [
    line.strip() for line in sarscov2_accession_file.readlines()
  ]
SARSCOV2_REFERENCES = [
  filename.split('.')[0]
  for filename in os.listdir('output/references')
  if filename[:4] == 'sars'
]
wildcard_constraints:
  dataset="[^/]+",
  simulated_dataset="[^/]+"

rule sarscov2:
  input:
    expand(
      "output/{dataset}/trimmomatic/bowtie2/{reference}/{haplotyper}/haplotypes.fasta",
      dataset=SARSCOV2_ACCESSIONS[:100],
      reference=SARSCOV2_REFERENCES,
      haplotyper=HAPLOTYPERS
    )

rule all_haplotypers:
  input:
    expand(
      "output/{dataset}/fastp/bowtie2/{reference}/{haplotyper}/haplotypes.fasta",
      dataset=ALL_DATASETS,
      reference=REFERENCE_SUBSET,
      haplotyper=HAPLOTYPERS
    ),
    expand(
      "output/{dataset}/reads_fastqc.html",
      dataset=ALL_DATASETS
    ),
    expand(
      "output/{dataset}/qfilt/bealign/{reference}/qualimapReport.html",
      dataset=ALL_DATASETS,
      reference=REFERENCE_SUBSET
    ),
    expand(
      "output/{dataset}/fastp/bowtie2/{reference}/qualimapReport.html",
      dataset=ALL_DATASETS,
      reference=REFERENCE_SUBSET
    )

rule all_acme_running:
  input:
    "output/sim-divergedFive_ar-10_seed-1/fastp/bowtie2/pol/acme/truth_and_haplotypes.json"
  output:
    "output/running_report.csv"
  run:
    report(input, output[0], 'running')

rule all_acme_reconstructing:
  input:
    "output/sim-divergedPair_ar-10_seed-1/fastp/bowtie2/pol/acme/truth_and_haplotypes.json",
    "output/sim-divergedTriplet_ar-10_seed-1/fastp/bowtie2/pol/acme/truth_and_haplotypes.json"
  output:
    "output/reconstruction_report.csv"
  run:
    report(input, output[0], 'reconstructing')

rule report:
  input:
    reconstructing=rules.all_acme_reconstructing.output[0],
    running=rules.all_acme_running.output[0]
  output:
    "output/report.csv"
  shell:
    """
      cat {input.reconstructing} > {output}
      tail -n +2 {input.running} >> {output}
    """

rule haplotyper_truth_report:
  input:
    expand(
      "output/{dataset}/{{qc}}/{{read_mapper}}/{reference}/{{haplotyper}}/truth_and_haplotypes.png",
      dataset=KNOWN_TRUTH,
      reference=REFERENCE_SUBSET
    )
  output:
    "output/reports/{haplotyper}-{qc}-{read_mapper}.csv"
  run:
    haplotyper_report(input, output[0])

rule all_bams:
  input:
    expand(
      "output/{dataset}/fastp/bowtie2/{reference}/sorted.bam",
      dataset=ALL_DATASETS,
      reference=REFERENCE_SUBSET,
      haplotyper=HAPLOTYPERS
    )

rule all_known_comparisons:
  input:
    expand(
      "output/{dataset}/fastp/bowtie2/{reference}/{haplotyper}/truth_and_haplotypes.png",
      dataset=KNOWN_TRUTH,
      reference=REFERENCE_SUBSET,
      haplotyper=HAPLOTYPERS
    )

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

rule extract_references:
   output:
     "output/references/reference_list.csv"
   run:
     extract_genbank_references(output[0])

rule extract_gene:
  input:
    rules.extract_references.output[0],
    genome=rules.extract_lanl_genome.output[0]
  output:
    sam="output/lanl/{lanl_id}/{gene}/sequence.sam",
    fasta="output/lanl/{lanl_id}/{gene}/sequence.fasta"
  conda:
    "envs/veg.yml"
  shell:
    """
      bealign -r {input.reference} {input.genome} {output.sam}
      bam2msa {output.sam} {output.fasta}
    """

rule simulation:
  input:
    rules.extract_lanl_genome.output[0]
  output:
    fastq="output/lanl/{lanl_id}/wgs.fastq",
    sam="output/lanl/{lanl_id}/wgs.sam",
    bam="output/lanl/{lanl_id}/wgs.bam",
    sorted_bam="output/lanl/{lanl_id}/sorted.bam",
    index="output/lanl/{lanl_id}/sorted.bam.bai",
    report="output/lanl/{lanl_id}/qualimapReport.html"
  params:
    out="output/lanl/{lanl_id}/wgs",
    directory="output/lanl/{lanl_id}"
  conda:
    "envs/ngs.yml"
  shell:
    """
      art_illumina -rs 1 -ss HS25 --samout -i {input} -l 140 -s 50 -c %d -o {params.out}
      mv {params.out}.fq {output.fastq}
      samtools view -Sb {output.sam} > {output.bam}
      samtools sort {output.bam} > {output.sorted_bam}
      samtools index {output.sorted_bam}
      qualimap bamqc -bam {output.sorted_bam} -outdir {params.directory}
    """ % NUMBER_OF_READS

def simulation_truth_input(wildcards):
  dataset = SIMULATION_INFORMATION[wildcards.simulated_dataset]
  lanl_ids = [info['lanl_id'] for info in dataset]
  genomes = ["output/lanl/%s/genome.fasta" % (lanl_id) for lanl_id in lanl_ids]
  return genomes

rule simulation_truth:
  input:
    simulation_truth_input
  output:
    "output/truth/sim-{simulated_dataset}/genomes.fasta"
  run:
    simulation_truth(wildcards.simulated_dataset, output[0])

rule simulation_truth_aligned:
  input:
    rules.simulation_truth.output[0],
  output:
    fasta="output/truth/sim-{simulated_dataset}/aligned.fasta",
    progress=os.getcwd()+"/output/truth/sim-{simulated_dataset}/progress.txt"
  shell:
    "mafft --progress {output.progress} {input} > {output.fasta}"

def wgs_simulation_inputs(wildcards):
  dataset = SIMULATION_INFORMATION[wildcards.simulated_dataset]
  lanl_ids = [info['lanl_id'] for info in dataset]
  reads = ["output/lanl/%s/wgs.fastq" % (lanl_id) for lanl_id in lanl_ids]
  genomes = ["output/lanl/%s/genome.fasta" % (lanl_id) for lanl_id in lanl_ids]
  return reads + genomes

rule simulate_wgs_dataset:
  input:
    wgs_simulation_inputs,
    fasta=rules.simulation_truth_aligned.output[0]
  output:
    fastq=temp("output/sim-{simulated_dataset}_ar-{ar}_seed-{seed}/wgs.fastq"),
    sam="output/sim-{simulated_dataset}_ar-{ar}_seed-{seed}/wgs.sam",
    json="output/sim-{simulated_dataset}_ar-{ar}_seed-{seed}/simulation_quality.json"
  run:
    simulate_wgs_dataset(
      wildcards.simulated_dataset, wildcards.ar, input.fasta,
      output.fastq, output.json, output.sam, wildcards.seed, NUMBER_OF_READS
    )

rule add_reference:
  input:
    truth=rules.simulation_truth_aligned.output[0],
    reference="output/references/{reference}.fasta"
  output:
    "output/truth/{simulated_dataset}/{reference}_added.fasta"
  shell:
    "mafft --add {input.reference} {input.truth} > {output}"

rule get_coordinates:
  input:
    truth=rules.simulation_truth_aligned.output[0],
    added=rules.add_reference.output[0]
  output:
    "output/truth/{simulated_dataset}/{reference}_coordinates.csv"
  run:
    get_coordinates(input.truth, input.added, output[0])

rule get_true_coverage:
  input:
    indicial_map=rules.get_coordinates.output[0],
    truth_sam=rules.simulate_wgs_dataset.output.sam
  output:
    "output/sim-{simulated_dataset}_ar-{ar}_seed-{seed}/{reference}_coverage.csv"
  run:
    get_true_coverage(input.indicial_map, input.truth_sam, output[0])

rule get_true_coverage_plot:
  input:
    rules.get_true_coverage.output[0]
  output:
    "output/sim-{simulated_dataset}_ar-{ar}_seed-{seed}/{reference}_coverage.png"
  script:
    "R/coverage_plot.R"

# Situating other data

rule sra_dataset:
  output:
    "output/sra/{sra_accession}.fastq"
  params:
    sra_dir="output/sra",
    first_pair="output/sra/{sra_accession}_1.fastq",
    second_pair="output/sra/{sra_accession}_2.fastq"
  shell:
    """
      prefetch -O output/sra {wildcards.sra_accession}
      fasterq-dump --outdir {params.sra_dir} {wildcards.sra_accession}
      if [ -f {params.first_pair} ]; then
        cat {params.first_pair} {params.second_pair} > {output}
      fi
    """

def situate_input(wildcards):
  dataset = wildcards.dataset
  is_evolution_dataset = dataset[:7] == 'ERS6610'
  is_amplicon_dataset = 'amplicon' in dataset
  is_simulated_dataset = 'sim-' in dataset
  is_sra_dataset = dataset[:2] == 'SR'

  if is_evolution_dataset:
    return "input/evolution/%s.fastq" % dataset
  if is_amplicon_dataset:
    return "output/%s/amplicon.fastq" % dataset
  if is_simulated_dataset:
    return "output/%s/wgs.fastq" % dataset
  if is_sra_dataset:
    return "output/sra/%s.fastq" % dataset
  return "input/reconstruction/%s.fastq" % dataset

rule situate_data:
  input:
    situate_input
  output:
    "output/{dataset}/reads.fastq"
  shell:
    "cp {input} {output}"

# Quality control

rule qfilt:
  input:
    rules.situate_data.output[0]
  output:
    fasta="output/{dataset}/qfilt/qc.fasta",
    json="output/{dataset}/qfilt/qc.json",
    html="output/{dataset}/reads_fastqc.html"
  params:
    dir="output/{dataset}"
  conda:
    "envs/veg.yml"
  shell:
    """
      qfilt -Q {input} -q 20 -l 50 -j >> {output.fasta} 2>> {output.json}
      fastqc {input} -o {params.dir}
    """

def fasta_454(wildcards):
  if wildcards.dataset == 'example_454':
    return "input/reconstruction/3.GAC.454Reads.fna"
  head = "input/compartmentalization/"
  tail = "/".join(wildcards.dataset.split('-'))
  return head + tail + "/reads.fasta"

def qual_454(wildcards):
  if wildcards.dataset == 'example_454':
    return "input/reconstruction/3.GAC.454Reads.qual"
  head = "input/compartmentalization/"
  tail = "/".join(wildcards.dataset.split('-'))
  return head + tail + "/scores.qual"

rule qfilt_454:
  input:
    fasta=fasta_454,
    qual=qual_454
  output:
    fasta="output/{dataset}/qfilt/qc.fasta",
    json="output/{dataset}/qfilt/qc.json"
  conda:
    "envs/veg.yml"
  shell:
    "qfilt -F {input.fasta} {input.qual} -q 20 -l 50 -j >> {output.fasta} 2>> {output.json}"

rule fastp:
  input:
    rules.situate_data.output[0]
  output:
    fastq="output/{dataset}/fastp/qc.fastq",
    fasta="output/{dataset}/fastp/qc.fasta",
    json="output/{dataset}/fastp/qc.json",
    html="output/{dataset}/fastp/qc.html"
  conda:
    "envs/ngs.yml"
  shell:
    """
      fastp -A -q 30 -i {input} -o {output.fastq} -j {output.json} -h {output.html} -n 50
      cat {output.fastq} | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > {output.fasta}
    """

rule trimmomatic:
  input:
    rules.situate_data.output[0]
  output:
    "output/{dataset}/trimmomatic/qc.fastq"
  conda:
    "envs/ngs.yml"
  shell:
    "trimmomatic SE {input} {output} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

# Read mapping

rule bealign:
  input:
    qc="output/{dataset}/{qc}/qc.fasta",
    reference="output/references/{reference}.fasta"
  output:
    bam="output/{dataset}/{qc}/bealign/{reference}/mapped.bam",
    discards="output/{dataset}/{qc}/bealign/{reference}/discards.fasta"
  conda:
    "envs/veg.yml"
  shell:
    "bealign -r {input.reference} -e 0.5 -m HIV_BETWEEN_F -D {output.discards} -R {input.qc} {output.bam}"

rule bwa:
  input:
    fastq="output/{dataset}/{qc}/qc.fastq",
    reference="output/references/{reference}.fasta"
  output:
    "output/{dataset}/{qc}/bwa/{reference}/mapped.bam"
  conda:
    "envs/ngs.yml"
  shell:
    """
      bwa index {input.reference}
      bwa mem {input.reference} {input.fastq} > {output}
    """

rule bowtie2_index:
  input:
    reference="output/references/{reference}.fasta"
  output:
    bt1="output/references/{reference}.1.bt2",
    bt2="output/references/{reference}.2.bt2",
    bt3="output/references/{reference}.3.bt2",
    bt4="output/references/{reference}.4.bt2",
    revbt1="output/references/{reference}.rev.1.bt2",
    revbt2="output/references/{reference}.rev.2.bt2"
  params:
    lambda wildcards: "output/references/%s" % wildcards.reference
  shell:
    "bowtie2-build {input.reference} {params}"

rule bowtie2_parameter_study:
  input:
    fastq="output/{dataset}/{qc}/qc.fastq",
    bt1=rules.bowtie2_index.output.bt1,
    bt2=rules.bowtie2_index.output.bt2,
    bt3=rules.bowtie2_index.output.bt3,
    bt4=rules.bowtie2_index.output.bt4,
    revbt1=rules.bowtie2_index.output.revbt1,
    revbt2=rules.bowtie2_index.output.revbt2
  output:
    sam="output/{dataset}/{qc}/bowtie2/{reference}/mapped_b-{b}_m-{m}_q-{q}.sam",
    bam="output/{dataset}/{qc}/bowtie2/{reference}/mapped_b-{b}_m-{m}_q-{q}.bam"
  params:
    index_loc=lambda wildcards: "output/references/%s" % wildcards.reference,
    score=lambda wildcards: "L,-%s,-%s" % (wildcards.b, wildcards.m)
  conda:
    "envs/ngs.yml"
  shell:
    """
      bowtie2 -x {params.index_loc} -U {input.fastq} -S {output.sam} \
        --score-min {params.score}
      samtools view -Sbq {wildcards.q} {output.sam} > {output.bam}
    """

rule bowtie2_alignment:
  input:
    "output/{dataset}/{qc}/bowtie2/{reference}/mapped_b-1_m-1_q-1.bam"
  output:
    "output/{dataset}/{qc}/bowtie2/{reference}/mapped.bam"
  shell:
    "cp {input} {output}"

rule minimap2:
  input:
    reads="output/{dataset}/{qc}/qc.fastq",
    reference="output/references/{reference}.fasta"
  output:
    "output/{dataset}/{qc}/minimap2/{reference}/mapped.bam"
  shell:
  	"minimap2 -a {input.reference} {input.reads} | samtools view -Sb > {output}"

rule read_statistics:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/mapped.bam"
  output:
    full="output/{dataset}/{qc}/{read_mapper}/{reference}/read-stats.csv",
    summary="output/{dataset}/{qc}/{read_mapper}/{reference}/read-stat-summary.csv"
  run:
    read_statistics(input[0], output.full, output.summary)

rule filter:
  input:
    bam="output/{dataset}/{qc}/{read_mapper}/{reference}/mapped.bam"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/filtered.bam"
  run:
    filter_bam(input.bam, output[0])

rule sort_and_index:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/mapped.bam"
  output:
    bam="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam",
    sam="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.sam",
    index="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam.bai"
  conda:
    "envs/ngs.yml"
  shell:
    """
      samtools sort {input} > {output.bam}
      samtools view -h {output.bam} > {output.sam}
      samtools index {output.bam}
    """

rule sorted_fasta:
  input:
    rules.sort_and_index.output.bam
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.fasta"
  conda:
    "envs/veg.yml"
  shell:
    "bam2msa {input} {output}"

rule deduplicated:
  input:
    rules.sort_and_index.output.bam
  output:
    bam="output/{dataset}/{qc}/{read_mapper}/{reference}/deduped.bam",
    metrics="output/{dataset}/{qc}/{read_mapper}/{reference}/picard_metrics.txt"
  shell:
    """
      picard MarkDuplicates \
        I={input} O={output.bam} M={output.metrics} \
        REMOVE_DUPLICATES=true ASSUME_SORTED=true
    """

rule realigned:
  input:
    bam=rules.deduplicated.output.bam,
    reference=rules.bowtie2_index.input.reference
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/realigned.bam"
  shell:
    "lofreq viterbi -f {input.reference} -o {output} {input.bam}"

rule call_variants:
  input:
    bam=rules.realigned.output[0],
    reference=rules.bowtie2_index.input.reference
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/variants.vcf"
  shell:
    "lofreq call -f {input.reference} -o {output} {input.bam}"

rule qualimap:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/qualimapReport.html"
  params:
    dir="output/{dataset}/{qc}/{read_mapper}/{reference}"
  conda:
    "envs/ngs.yml"
  shell:
    "qualimap bamqc -bam {input} -outdir {params.dir}"

rule insertion_plot:
  input:
    rules.sort_and_index.output.bam
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/insertion_plot.png"
  script:
    "R/insertion_plot.R"

# Haplotype reconstruction (full pipelines)

rule regress_haplo_full:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam",
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam.bai",
  output:
    temp("output/{dataset}/{qc}/{read_mapper}/{reference}/regress_haplo/final_haplo.fasta")
  script:
    "R/regress_haplo/full_pipeline.R"

rule regress_haplo_rightname:
  input:
    rules.regress_haplo_full.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/regress_haplo/haplotypes.fasta"
  shell:
    "mv {input} {output}"

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

rule quasirecomb_jar:
  output:
    "QuasiRecomb.jar"
  shell:
    "wget https://github.com/cbg-ethz/QuasiRecomb/releases/download/v1.2/QuasiRecomb.jar"

rule quasirecomb:
  input:
    jar=rules.quasirecomb_jar.output[0],
    bam="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/quasirecomb/haplotypes.fasta"
  params:
    basedir="output/{dataset}/{qc}/{read_mapper}/{reference}/quasirecomb"
  conda:
    "envs/quasirecomb.yml"
  shell:
    """
      java -jar QuasiRecomb.jar -conservative -o {params.basedir} -i {input.bam}
      mv {params.basedir}/quasispecies.fasta {params.basedir}/haplotypes.fasta
    """

rule savage:
  input:
    bam="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.bam",
    reference="input/references/{reference}.fasta"
  output:
    fastq="output/{dataset}/{qc}/{read_mapper}/{reference}/savage/reads.fastq",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/savage/haplotypes.fasta"
  params:
    outdir="output/{dataset}/{qc}/{read_mapper}/{reference}/savage",
    intermediate="output/{dataset}/{qc}/{read_mapper}/{reference}/savage/contigs_stage_c.fasta"
  shell:
    """
      bamToFastq -i {input.bam} -fq {output.fastq}
      savage -s {output.fastq} --ref `pwd`/{input.reference} --split 3 --num_threads 12 --outdir {params.outdir}
      mv {params.intermediate} {output.fasta}
    """

rule abayesqr_config:
  input:
    sam="output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.sam",
    reference="output/references/{reference}.fasta"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/config",
  run:
    write_abayesqr_config(input.sam, input.reference, output[0])

rule abayesqr:
  input:
    rules.abayesqr_config.output[0]
  output:
    freq="output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/test_Freq.txt",
    seq="output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/test_Seq.txt",
    viralseq="output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/test_ViralSeq.txt",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/abayesqr/haplotypes.fasta"
  run:
    shell("aBayesQR {input}")
    shell("mv test_Freq.txt {output.freq}")
    shell("mv test_Seq.txt {output.seq}")
    shell("mv test_ViralSeq.txt {output.viralseq}")
    parse_abayesqr_output(output.viralseq, output.fasta)

rule shorah:
  input:
    bam=rules.sort_and_index.output.bam,
    reference="output/references/{reference}.fasta"
  output:
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/shorah/haplotypes.fasta"
  params:
    workdir="output/{dataset}/{qc}/{read_mapper}/{reference}/shorah",
    bam="../sorted.bam",
    reference="../../../../../references/{reference}.fasta"
  conda:
    "envs/shorah.yml"
  shell:
    """
      cd {params.workdir}
      shorah.py -b {params.bam} -f {params.reference} -w 51
      mv sorted_global_haps.fasta haplotypes.fasta
    """

# VEG haplotype reconstruction

rule mmvc:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/sorted.fasta"
  output:
    json="output/{dataset}/{qc}/{read_mapper}/{reference}/mmvc.json",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/mmvc.fasta"
  conda:
    "envs/veg.yml"
  shell:
    "mmvc -j {output.json} -f {output.fasta} {input}"

rule readreduce:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/mmvc.fasta"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/veg/haplosuperreads.fasta"
  conda:
    "envs/veg.yml"
  shell:
    "readreduce -a resolve -l 30 -s 16 -o {output} {input}"

# ACME haplotype reconstruction

WEIGHT_FILTER = 5

def true_sequences_input(wildcards):
  if wildcards.known_dataset[:len(FVM_STRING)] == FVM_STRING:
    return "input/5VM.fasta"
  dataset = wildcards.known_dataset.split('_')[0]
  return "output/truth/%s/genomes.fasta" % dataset

rule true_sequences:
  input:
    wgs=true_sequences_input,
    reference="output/references/{reference}.fasta"
  output:
    fasta="output/truth/{known_dataset}/{reference}_gene.fasta",
    json="output/truth/{known_dataset}/{reference}_gene.json"
  run:
    extract_truth(input.wgs, input.reference, wildcards.known_dataset, wildcards.reference, output.fasta, output.json)

rule covarying_sites:
  input:
    rules.sort_and_index.output.bam
  output:
    json="output/{dataset}/{qc}/{read_mapper}/{reference}/acme/covarying_sites.json",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/acme/consensus.fasta",
    covariation="output/{dataset}/{qc}/{read_mapper}/{reference}/acme/covariation.csv"
  run:
    covarying_sites_io(input[0], output.json, output.fasta, output.covariation, threshold=.05)

rule true_covarying_sites:
  input:
    rules.true_sequences.output.fasta
  output:
    "output/truth/{known_dataset}/{reference}_covarying_sites.json"
  run:
    covarying_sites(input[0], output[0])

def covarying_truth_comparison_input(wildcards):
    dataset = wildcards.dataset
    if dataset[:len(FVM_STRING)] != FVM_STRING:
      known_dataset = dataset.split('_')[0]
    else:
      known_dataset = dataset
    computed_string = "output/%s/%s/%s/%s/acme/covarying_sites.json"
    computed_arguments = (dataset, wildcards.qc, wildcards.read_mapper, wildcards.reference)
    computed = computed_string % computed_arguments
    truth = "output/truth/%s/%s_covarying_sites.json" % (known_dataset, wildcards.reference)
    reference = "output/references/%s.fasta" % wildcards.reference
    return [computed, truth, reference]

rule covarying_truth_comparison:
  input:
    covarying_truth_comparison_input
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/covarying_truth.json"
  run:
    covarying_truth(input[0], input[1], input[2], output[0])

rule superreads:
  input:
    alignment=rules.sort_and_index.output.bam,
    covarying_sites=rules.covarying_sites.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/superreads.json",
  run:
    superread_json_io(input.alignment, input.covarying_sites, output[0])

rule superread_fasta:
  input:
    cvs=rules.covarying_sites.output[0],
    sr=rules.superreads.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/superreads.fasta"
  run:
    superread_fasta_io(
      input.cvs, input.sr, output[0], weight_filter=WEIGHT_FILTER
    )

rule superread_scatter_data:
  input:
    rules.superreads.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/superreads.csv"
  run:
    superread_scatter_data(input[0], output[0])

rule superread_scatter_plot:
  input:
    rules.superread_scatter_data.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/superread_scatter.png"
  script:
    "R/superread_scatterplot.R"


def inspect_superread_input(wildcards):
    dataset = wildcards.dataset
    if dataset[:len(FVM_STRING)] != FVM_STRING:
      known_dataset = dataset.split('_')[0]
    else:
      known_dataset = dataset
    bam_string = "output/%s/%s/%s/%s/sorted.bam"
    bam_arguments = (dataset, wildcards.qc, wildcards.read_mapper, wildcards.reference)
    bam = bam_string % bam_arguments
    truth = "output/truth/%s/%s_gene.fasta" % (known_dataset, wildcards.reference)
    sr_string = "output/%s/%s/%s/%s/acme/superreads.json"
    sr_arguments = (dataset, wildcards.qc, wildcards.read_mapper, wildcards.reference)
    superreads = sr_string % sr_arguments
    truth = "output/references/%s.fasta" % wildcards.reference
    return [bam, superreads, truth]


rule inspect_superread:
  input:
    inspect_superread_input
  output:
    bam="output/{dataset}/{qc}/{read_mapper}/{reference}/acme/sr-{sr}/sorted.bam",
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/acme/sr-{sr}/sorted.fasta",
    ref="output/{dataset}/{qc}/{read_mapper}/{reference}/acme/sr-{sr}/sorted-ref.fasta"
  run:
    pluck_superread_reads(input[1], int(wildcards.sr), input[0], output.bam)
    shell("bam2msa {output.bam} {output.fasta}")
    shell("cat {input[2]} {output.fasta} > {output.ref}")

rule inspect_superread_at_cvs:
  input:
    superread=rules.inspect_superread.output.ref,
    covarying_sites=rules.covarying_sites.output.json
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/sr-{sr}/sorted-ref-cvs.fasta"
  run:
    restrict_fasta_to_cvs(input.superread, input.covarying_sites, output[0])

def sorted_fasta_and_truth_input(wildcards):
  dataset = wildcards.dataset
  if dataset[:len(FVM_STRING)] != FVM_STRING:
    known_dataset = dataset.split('_')[0]
  else:
    known_dataset = dataset
  computed_string = "output/%s/%s/%s/%s/sorted.fasta"
  computed_arguments = (dataset, wildcards.qc, wildcards.read_mapper, wildcards.reference)
  computed = computed_string % computed_arguments
  truth = "output/truth/%s/%s_gene.fasta" % (known_dataset, wildcards.reference)
  return [truth, computed]


rule sorted_fasta_and_truth:
  input:
    sorted_fasta_and_truth_input
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/truth_and_sorted.fasta"
  shell:
    "cat {input[0]} {input[1]} > {output}"

def truth_at_cvs_input(wildcards):
  dataset = wildcards.dataset
  if dataset[:len(FVM_STRING)] != FVM_STRING:
    known_dataset = dataset.split('_')[0]
  else:
    known_dataset = dataset
  computed_string = "output/%s/%s/%s/%s/acme/covarying_sites.json"
  computed_arguments = (dataset, wildcards.qc, wildcards.read_mapper, wildcards.reference)
  computed = computed_string % computed_arguments
  truth = "output/truth/%s/%s_gene.fasta" % (known_dataset, wildcards.reference)
  return [truth, computed]

rule truth_at_cvs:
  input:
    truth_at_cvs_input
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/truth-cvs.fasta"
  run:
    restrict_fasta_to_cvs(input[0], input[1], output[0])

rule truth_and_superreads_cvs:
  input:
    truth=rules.truth_at_cvs.output[0],
    superreads=rules.superread_fasta.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/truth-and-superreads.fasta"
  shell:
    "cat {input.truth} {input.superreads} > {output}"

rule quantify_ar:
  input:
    fasta=rules.truth_at_cvs.output[0],
    json=rules.superreads.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/truth-sr-cvs.json"
  run:
    quantify_ar(input[0], input[1], output[0])

rule consolidate_sim_quantify_ar:
  input:
    json=rules.quantify_ar.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/similar_consolidated.csv"
  run:
    consolidate_simjsons(input[0], output[0])

rule superread_graph:
  input:
    rules.superreads.output[0],
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/graph-full.json"
  run:
    graph_io(input[0], output[0], 'full')

rule reduced_superread_graph:
  input:
    rules.superreads.output[0],
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/graph-reduced.json"
  run:
    graph_io(input[0], output[0], 'reduced', minimum_weight=WEIGHT_FILTER)

rule incremental_superread_graph:
  input:
    rules.superreads.output[0],
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/graph-incremental.json"
  run:
    graph_io(input[0], output[0], 'incremental')

rule candidates:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/graph-{graph_type}.json",
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/describing-{graph_type}.json"
  run:
    candidates_io(input[0], output[0])

rule regression:
  input:
    superreads=rules.superreads.output[0],
    describing=rules.candidates.output[0],
    consensus=rules.covarying_sites.output.fasta,
    covarying_sites=rules.covarying_sites.output.json
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/haplotypes-{graph_type}.fasta"
  run:
    regression_io(
      input.superreads, input.describing, input.consensus,
      input.covarying_sites, output[0]
    )

rule chosen_approach:
  input:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/haplotypes-reduced.fasta"
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/haplotypes.fasta"
  shell:
    "cp {input} {output}"

rule edge_list:
  input:
    rules.superreads.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/edge_list.csv"
  run:
    edge_list_io(input[0], output[0])

rule scaffold_describing:
  input:
    rules.superreads.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/scaffold_describing.json"
  run:
    scaffold_qsr_io(input[0], output[0])

rule scaffold_candidates:
  input:
    superreads=rules.superreads.output[0],
    description=rules.scaffold_describing.output[0]
  output:
    fasta="output/{dataset}/{qc}/{read_mapper}/{reference}/acme/scaffold_candidates.fasta",
    csv="output/{dataset}/{qc}/{read_mapper}/{reference}/acme/scaffold_candidates.csv"
  run:
    scaffold_candidates_io(
      input.superreads, input.description,
      output.fasta, output.csv
    )

# Simulation studies

rule simulation_coverage:
  input:
    rules.bowtie2_parameter_study.output.bam
  output:
    "output/{dataset}/{qc}/bowtie2/{reference}/coverage_b-{b}_m-{m}_q-{q}.csv"
  run:
    simulation_coverage(
      input[0], output[0],
      wildcards.b, wildcards.m, wildcards.q
    )

rule master_coverage_data:
  input:
    expand(
      "output/{{dataset}}/{{qc}}/bowtie2/{{reference}}/coverage_b-{b}_m-{m}_q-{q}.csv",
      b=[.6, 1.2, 1.8],
      m=[.6, 1.2, 1.8],
      q=[0, 20, 40]
    )
  output:
    "output/{dataset}/{qc}/bowtie2/{reference}/master-coverage.csv"
  shell:
    "cat {input} > {output}"

rule master_coverage_plot:
  input:
    rules.master_coverage_data.output[0]
  output:
    "output/{dataset}/{qc}/bowtie2/{reference}/master-coverage.png"
  script:
    "R/master_coverage_plot.R"

rule generic_coverage_data:
  input:
    rules.sort_and_index.output.bam
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/coverage.csv"
  run:
    simulation_coverage(input[0], output[0])

rule generic_coverage_plot:
  input:
    rules.generic_coverage_data.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/coverage.png"
  script:
    "R/coverage_plot.R"

def n_paths_boxplot_input(wildcards):
  template_string = "output/sim-%s_ar-%d_seed-%d/fastp/bowtie2/%s/acme/graph.json"
  input_files = []
  for seed in range(1, 11):
    for ar in [0, 5, 10, 15, 20]:
      parameters = (wildcards.simulated_dataset, ar, seed, wildcards.gene)
      input_files.append(template_string % parameters)
  return input_files

rule n_paths_boxplot:
  input:
    n_paths_boxplot_input
  output:
    "output/simulation/{simulated_dataset}/n_paths_boxplot_{gene}.png"
  run:
    n_paths_boxplot(wildcards.simulated_dataset, wildcards.gene, output[0])

rule superread_weight_distribution_data:
  input:
    rules.superreads.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/superread_distribution.csv"
  run:
    superread_weight_distribution_data(input[0], output[0])

rule superread_weight_distribution_plot:
  input:
    rules.superread_weight_distribution_data.output[0]
  output:
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/superread_histogram.png",
    "output/{dataset}/{qc}/{read_mapper}/{reference}/acme/superread_scatter_plot.png"
  script:
    "R/superread_distribution.R"

def reference_input(wildcards):
  format_string = "output/truth/%s/%s_gene.fasta"
  if wildcards.dataset[:3] == 'sim':
    dataset = wildcards.dataset.split('_')[0]
  else:
    dataset = wildcards.dataset
  parameters = (dataset, wildcards.reference)
  return format_string % parameters

rule haplotypes_and_truth:
  input:
    haplotypes="output/{dataset}/{qc}/{read_mapper}/{reference}/{haplotyper}/haplotypes.fasta",
    truth=reference_input
  output:
    unaligned="output/{dataset}/{qc}/{read_mapper}/{reference}/{haplotyper}/truth_and_haplotypes_unaligned.fasta",
    aligned="output/{dataset}/{qc}/{read_mapper}/{reference}/{haplotyper}/truth_and_haplotypes.fasta",
    progress=os.getcwd()+"output/{dataset}/{qc}/{read_mapper}/{reference}/{haplotyper}/progress.txt",
    csv="output/{dataset}/{qc}/{read_mapper}/{reference}/{haplotyper}/truth_and_haplotypes.csv",
    json="output/{dataset}/{qc}/{read_mapper}/{reference}/{haplotyper}/truth_and_haplotypes.json"
  run:
    shell("cat {input.haplotypes} {input.truth} > {output.unaligned}")
    shell("mafft --progress {output.progress} {output.unaligned} > {output.aligned}")
    pairwise_distance_csv(output.aligned, output.csv)
    result_json(output.csv, output.json)

rule haplotypes_and_truth_heatmap:
  input:
    rules.haplotypes_and_truth.output.csv
  output:
    png="output/{dataset}/{qc}/{read_mapper}/{reference}/{haplotyper}/truth_and_haplotypes.png"
  conda:
    "envs/R.yml"
  script:
    "R/truth_heatmap.R"

# Regress Haplo

rule regress_haplo_bam_to_variant_calls:
  input:
    "output/{dataset}/{reference}/sorted.bam",
    "output/{dataset}/{reference}/sorted.bam.bai"
  output:
    "output/{dataset}/{reference}/variant_calls.csv"
  script:
    "R/regress_haplo/bam_to_variant_calls.R"
   
rule regress_haplo_variant_calls_to_read_table:
  input:
    "output/{dataset}/{reference}/sorted.bam",
    "output/{dataset}/{reference}/variant_calls.csv",
  output:
    "output/{dataset}/{reference}/read_table.csv"
  script:
    "R/regress_haplo/variant_calls_to_read_table.R"

rule regress_haplo_read_table_to_loci:
  input:
    rules.regress_haplo_variant_calls_to_read_table.output[0]
  output:
    "output/{dataset}/{reference}/loci.csv"
  script:
    "R/regress_haplo/read_table_to_loci.R"

rule regress_haplo_loci_to_haplotypes:
  input:
    rules.regress_haplo_read_table_to_loci.output[0]
  output:
    "output/{dataset}/{reference}/h.csv"
  script:
    "R/regress_haplo/loci_to_haplotypes.R"

rule regress_haplo_haplotypes_to_parameters:
  input:
    rules.regress_haplo_loci_to_haplotypes.output[0]
  output:
    "output/{dataset}/{reference}/P.csv"
  script:
    "R/regress_haplo/haplotypes_to_parameters.R"

rule regress_haplo_parameters_to_solutions:
  input:
    rules.regress_haplo_haplotypes_to_parameters.output[0]
  output:
    "output/{dataset}/{reference}/solutions.csv"
  script:
    "R/regress_haplo/parameters_to_solutions.R"

rule regress_haplo_solutions_to_haplotypes:
  input:
    rules.regress_haplo_parameters_to_solutions.output[0]
  output:
    "output/{dataset}/{reference}/final_haplo.csv"
  script:
    "R/regress_haplo/solutions_to_haplotypes.R"

rule regress_haplo_haplotypes_to_fasta:
  input:
    rules.regress_haplo_bam_to_variant_calls.input[0],
    rules.regress_haplo_solutions_to_haplotypes.output[0]
  output:
    "output/{dataset}/{reference}/final_haplo.fasta"
  script:
    "R/regress_haplo/haplotypes_to_fasta.R"

#############
# EVOLUTION #
#############

rule concatenate:
  input:
    expand("output/{dataset}/qfilt/bealign/{{reference}}/{{haplotyper}}/haplotypes.fasta", dataset=ACCESSION_NUMBERS)
  output:
    "output/evolution/{qc}/{read_mapper}/{reference}/{haplotyper}/unaligned.fasta"
  params:
    lambda wildcards: ' '.join([
      "output/%s/qfilt/bealign/%s/%s/haplotypes.fasta" % 
      (accession, wildcards.reference, wildcards.haplotyper) for accession in ACCESSION_NUMBERS
    ])
  shell:
    "cat {params} > {output}"

rule alignment:
  input:
    rules.concatenate.output[0]
  output:
    "output/evolution/{qc}/{read_mapper}/{reference}/{haplotyper}/aligned.fasta"
  shell:
    "mafft {input} > {output}"

rule tree:
  input:
    rules.alignment.output[0]
  output:
    "output/evolution/{qc}/{read_mapper}/{reference}/{haplotyper}/tree.new"
  shell:
    "FastTree -nt {input} > {output}"
