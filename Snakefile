reads = ['ERS6610%d.fastq' for i in range(87, 94)]

rule reference_index:
  input:
    "input/{reference}.fasta"
  output:
    "output/references/{reference}.fasta"
  shell:
    """
      cp {input} {output}
      bwa index {output}
    """

rule map_reads:
  input:
    fastq="input/{accession}.fastq",
    reference=rules.reference_index.output
  output:
    "output/{accession}/{reference}/mapped.sam"
  shell:
    "bwa mem {input.reference} {input.fastq} > {output}"

rule sort_and_index:
  input:
    rules.map_reads.output
  output:
    "output/{accession}/{reference}/sorted.bam"
  shell:
    """
      samtools sort {input} > {output}
      samtools index {output}
    """

rule reconstruct_haplotypes:
  input:
    rules.sort_and_index.output
  output:
    "output/{accession}/{reference}/haplotypes/final_haplo.fasta"
  script:
    "invoke_regress_haplo.R"

