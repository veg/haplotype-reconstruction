# haplotype-reconstruction

Benchmarking of tools for viral haplotype reconstruction.

<b>WARNING</b>: under active construction! May not work as appears, but in such an event feel free to file an issue.

## Use

This repo can be used to test a variety of tools and settings to build pipelines for viral haplotype/quasispecies reconstruction on simulated, in vitro and in vivo data.

| Haplotyper | Implemented | Paper | Code |
| ------------- |:-------------:|:-----:|:-----:|
| QuasiRecomb  | Yes | [Link](https://www.liebertpub.com/doi/full/10.1089/cmb.2012.0232) | [Link](https://github.com/cbg-ethz/QuasiRecomb) |
| aBayesQR     | Yes | [Link](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_22) | [Link](https://github.com/SoyeonA/aBayesQR) |
| SAVAGE       | Yes | [Link](https://genome.cshlp.org/content/27/5/835.short) | [Link](https://bitbucket.org/jbaaijens/savage) |
| RegressHaplo | Yes | [Link](https://academic.oup.com/bioinformatics/article/33/16/2455/3100436) |   [Link](https://github.com/SLeviyang/RegressHaplo) |
| Haploclique  | Yes | [Link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003515) | [Link](https://github.com/armintoepfer/haploclique) | 
| SHORAH       | No | [Link](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-119) | [Link](https://github.com/cbg-ethz/shorah) | 
| PredictHaplo | No | [Link](https://ieeexplore.ieee.org/document/6661314) | [Link](https://bmda.dmi.unibas.ch/software.html) |
| QSdpR | No | [Link](https://www.sciencedirect.com/science/article/pii/S0888754317301568) | [Link](https://sourceforge.net/projects/qsdpr/) |
| TenSQR | No | [Link](https://academic.oup.com/bioinformatics/article/34/13/i23/5045739) | [Link](https://github.com/SoYeonA/TenSQR) |

A standard invocation will be of the form:

```
snakemake output/$DATASET/$QUALITY_CONTROL/$READ_MAPPER/$GENE/$HAPLOTYPER/haplotypes.fasta
```

where the variables can take one of the following values:

### Dataset
For more information, see the Data section below.

- <b>LANL simulations</b>: see keys of `simulations.json` for potential simulation names
- <b>Reconstruction</b>: see filenames in `reconstruction` directory of input data
- <b>Evolution</b>: see filenames in `evolution` directory of input data
- <b>Compartmentalization</b>: see entries in `compartmentalization.json` in root directory

### Quality Control
`qfilt`, `fastp`, `trimmomatic`

### Read mapper
`bealign`, `bowtie2`, `bwa`

### Gene
`env`, `gag`, `int`, `nef`, `pol`, `pr`, `prrt`, `rev`, `rt`, `tat`, `vif`, `vpr`

### Haplotyper
`quasirecomb`, `abayesqr`, `savage`, `regresshaplo`

Use of snakemake permits running on TORQUE, i.e.

```
snakemake --cluster 'qsub -o ./logs -e ./logs -V -d `pwd` -l nodes=1:ppn=$PPN' -j $JOBS -k target
```

## Data

### Input directory structure
<pre>
├── compartmentalization
│   ├── $PATIENT_ID/$DATE/$COMPARTMENT/$REPLICATE/reads.fasta
│   ├── $PATIENT_ID/$DATE/$COMPARTMENT/$REPLICATE/scores.qual
├── evolution
│   ├── ERS661087.fastq
│   ├── ERS661088.fastq
│   ├── ERS661089.fastq
│   ├── ERS661090.fastq
│   ├── ERS661091.fastq
│   ├── ERS661092.fastq
│   └── ERS661093.fastq
├── LANL-HIV-aligned.fasta
├── LANL-HIV.fasta
├── LANL-HIV.new
├── README.md
├── reconstruction
│   ├── 3.GAC.454Reads.fna
│   ├── 3.GAC.454Reads.qual
│   ├── 93US141_100k_14-159320-1GN-0_S16_L001_R1_001.fastq
│   ├── 93US141_100k_14-159320-1GN-0_S16_L001_R2_001.fastq
│   ├── BP_050100753.fasta
│   ├── BP_050100753.qual
│   ├── FiveVirusMixIllumina_1.fastq
│   ├── FiveVirusMixIllumina_2.fastq
│   ├── PP1L_S45_L001_R1_001.fastq
│   ├── PP1L_S45_L001_R2_001.fastq
│   ├── regress_haplo.bam
│   ├── regress_haplo.bam.bai
│   ├── sergei1.fastq
│   ├── sergei2.fastq
│   ├── SRR961514-Illumina.sra
│   ├── SRR961596-454.fastq
│   ├── SRR961596-454.sra
│   ├── SRR961669-PacBio.fastq
│   └── SRR961669-PacBio.sra
└── references
    ├── env.fasta
    ├── gag.fasta
    ├── int.fasta
    ├── nef.fasta
    ├── pol.fasta
    ├── pr.fasta
    ├── prrt.fasta
    ├── rev.fasta
    ├── rt.fasta
    ├── tat.fasta
    ├── vif.fasta
    ├── vpr.fasta
    └── vpu.fasta
</pre>

### Description

### LANL

```
LANL-HIV.fasta
LANL-HIV-aligned.fasta
LANL-HIV.new
```

HIV genomes from the [LANL database](https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html), as well as an alignment built with `mafft` and a tree built with `FastTree`. Used for simulation.

### Intrahost evolution

```
evolution/ERS6610*.fastq
```

NGS read data from a [study on HIV intra-host evolution](http://www.genetics.org/content/202/4/1449).

### 454 data

```
reconstruction/BP_050100753.fasta
reconstruction/BP_050100753.qual
```

ACME lab 454 data which shows a clear signal of segregating haplotypes.

### HIV benchmarking

```
reconstruction/SRR961514-Illumina.fastq
reconstruction/SRR961596-454.fastq
reconstruction/SRR961669-PacBio.fastq
```

A [gold standard dataset](https://github.com/cbg-ethz/5-virus-mix), consisting of mixed, known strains at known proportions.

### Sergei

```
reconstruction/sergei1.fastq
reconstruction/sergei2.fastq
```

A set of paired end reads given by Sergei.

### RegressHaplo test set

```
reconstruction/regress_haplo.bam
reconstruction/regress_haplo.bam.bai
```

Dataset that comes with the [RegressHaplo](https://github.com/SLeviyang/RegressHaplo) code.

### Reference genes

```
references/*.fasta
```

[HXB2](https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html) genes to be used as references when aligning reads.

## Installation

### Requirements

- Linux (tested on Ubuntu 18.04.1, CentOS 7)
- conda (tested on 4.6.10) with [standard BioConda channels](https://bioconda.github.io/#set-up-channels)

Further requirements listed in `environment.yml`.

### Install

```
conda env create -f environment.yml
conda activate haplotype-reconstruction
```

