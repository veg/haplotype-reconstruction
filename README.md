# haplotype-reconstruction

Benchmarking of tools for viral haplotype reconstruction.

## Installation

### Requirements

- conda
- virtualenv
- [art_illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)

### Install

Obtain a copy of `veg-haplo.tar.gz` and place in `data`.

Run:

```
bash install.sh
```

## Pipeline

Generate synthetic data with distinct haplotypes and assemble:

```
bash generate_and_assemble_synthetic.sh
```

