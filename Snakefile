import os


# Evolution
ACCESSION_NUMBERS = ['ERS6610%d' % i for i in range(87, 94)]
SIMULATED_DATASETS = ["related_1", "related_2", "diverged_1", "diverged_2"]
REFERENCES = ["env_C2V5", "gag_p24", "pr", "rt"]
HYPHY_PATH = "/Users/stephenshank/Software/lib/hyphy"

rule reference_index:
  input:
    "input/references/{reference}.fasta"
  output:
    "output/references/{reference}.fasta"
  shell:
    """
      cp {input} {output}
      bwa index {output}
    """

rule extract_lanl:
  input:
    "input/LANL-HIV.fasta"
  output:
    "output/simulation/related_1/sequence.fasta",
    "output/simulation/related_2/sequence.fasta",
    "output/simulation/diverged_1/sequence.fasta",
    "output/simulation/diverged_2/sequence.fasta"
  shell:
    "python python/extract_sequences.py"

rule simulate_single:
  input:
    "output/simulation/{simulated_dataset}/sequence.fasta"
  output:
    "output/{simulated_dataset}/reads.fastq",
  shell:
    """
      art_illumina -ss HS25 -i {input} -l 120 -s 50 -c 15000 -o {output}
      mv {output}.fq {output}
    """

rule simulate_mixed:
  input:
    "output/{mixed_dataset}_1/reads.fastq",
    "output/{mixed_dataset}_2/reads.fastq"
  output:
    "output/{mixed_dataset}_joint/reads.fastq",
  shell:
    """
      cat {input[0]} {input[1]} > {output}
    """

rule situate_intrahost_data:
  input:
    "input/evolution/{accession_number}.fastq"
  output:
    "output/{accession_number}/reads.fastq"
  shell:
    "cp {input} {output}"

rule all_simulated_reads:
  input:
    expand(
      "output/{simulated_dataset}/reads.fastq",
      simulated_dataset=SIMULATED_DATASETS
    ),
    expand(
      "output/{mixed_dataset}_joint/reads.fastq",
      mixed_dataset=['related', 'diverged']
    )

rule all_intrahost_reads:
  input:
    expand(
      "output/{accession_number}/reads.fastq",
      accession_number=ACCESSION_NUMBERS
    )

rule quality_control:
  input:
    "output/{dataset}/reads.fastq",
    rules.all_simulated_reads.input,
    rules.all_intrahost_reads.input
  output:
    fastq="output/{dataset}/qc.fastq",
    json="output/{dataset}/fastp.json",
    html="output/{dataset}/fastp.html"
  shell:
    "fastp -A -q 10 -i {input} -o {output.fastq} -j {output.json} -h {output.html}"

rule map_reads:
  input:
    fastq=rules.quality_control.output[0],
    reference=rules.reference_index.output
  output:
    "output/{dataset}/{reference}/mapped.sam"
  shell:
    "bwa mem {input.reference} {input.fastq} > {output}"

rule sort_and_index:
  input:
    rules.map_reads.output
  output:
    "output/{dataset}/{reference}/sorted.bam",
    "output/{dataset}/{reference}/sorted.bam.bai"
  shell:
    """
      samtools sort {input} > {output[0]}
      samtools index {output[0]}
    """

rule regress_haplo_full_pipeline:
  input:
    rules.sort_and_index.output
  output:
    "output/{dataset}/{reference}/full/final_haplo.fasta"
  script:
    "R/regress_haplo/full_pipeline.R"

rule concatenate:
  input:
    expand("output/{dataset}/{{reference}}/full/final_haplo.fasta", dataset=ACCESSION_NUMBERS)
  output:
    "output/{reference}/unaligned.fasta"
  params:
    lambda wildcards: ' '.join(["output/%s/%s/haplotypes/final_haplo.fasta" % (accession, wildcards.reference) for accession in ACCESSION_NUMBERS])
  shell:
    "cat {params} > {output}"

rule alignment:
  input:
    rules.concatenate.output[0]
  output:
    "output/{reference}/aligned.fasta"
  shell:
    "mafft {input} > {output}"

rule recombination_screening:
  input:
    rules.alignment.output[0]
  output:
    gard_json="output/{reference}/GARD.json",
    nexus="output/{reference}/seqs_and_trees.nex"
  params:
    gard_path="%s/TemplateBatchFiles/GARD.bf" % HYPHY_PATH,
    gard_output=os.getcwd() + "/output/{reference}/aligned.GARD",
    final_out=os.getcwd() + "/output/{reference}/aligned.GARD_finalout",
    translate_gard_j=os.getcwd() + "/output/{reference}/aligned.GARD.json",
    translated_json=os.getcwd() + "/output/{reference}/GARD.json",
    lib_path=HYPHY_PATH,
    alignment_path=os.getcwd() + "/output/{reference}/aligned.fasta"
  shell:
    """
      mpirun -np 2 HYPHYMPI LIBPATH={params.lib_path} {params.gard_path} {params.alignment_path} '010010' None {params.gard_output}
      translate-gard -i {params.gard_output} -j {params.translate_gard_j} -o {params.translated_json}
      mv {params.final_out} {output.nexus}
    """

rule site_selection:
  input:
    rules.recombination_screening.output.nexus
  output:
    "output/{reference}/seqs_and_trees.nex.FUBAR.json"
  params:
    full_nexus_path=os.getcwd() + "/" + rules.recombination_screening.output.nexus,
    fubar_path="%s/TemplateBatchFiles/SelectionAnalyses/FUBAR.bf" % HYPHY_PATH,
    lib_path=HYPHY_PATH
  shell:
    "(echo 1; echo {params.full_nexus_path}; echo 20; echo 1; echo 5; echo 2000000; echo 1000000; echo 100; echo .5;) | HYPHYMP LIBPATH={params.lib_path} {params.fubar_path}"

rule gene_selection:
  input:
    rules.recombination_screening.output.nexus
  output:
    "output/{reference}/seqs_and_trees.nex.BUSTED.json"
  params:
    full_nexus_path=os.getcwd() + "/" + rules.recombination_screening.output.nexus,
    busted_path="%s/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf" % HYPHY_PATH,
    lib_path=HYPHY_PATH
  shell:
    "(echo 1; echo {params.full_nexus_path}; echo 2;) | HYPHYMP LIBPATH={params.lib_path} {params.busted_path}"

rule full_analysis:
  input:
    rules.site_selection.output[0],
    rules.gene_selection.output[0]
  output:
    "output/{reference}/results.tar.gz"
  shell:
    "tar cvzf {output} {input[0]} {input[1]}"

# Reconstruction
rule situate_regress_haplo_example:
  input:
    "input/reconstruction/regress_haplo.bam",
    "input/reconstruction/regress_haplo.bam.bai"
  output:
    "output/regress_haplo_example/mac239env/sorted.bam",
    "output/regress_haplo_example/mac239env/sorted.bam.bai"
  shell:
    """
      cp {input[0]} {output[0]}
      cp {input[1]} {output[1]}
    """

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

