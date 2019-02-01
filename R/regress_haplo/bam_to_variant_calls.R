library(RegressHaplo)
split_path <- strsplit(snakemake@output[[1]], "/")[[1]]
directory <- paste(head(split_path, -1), collapse="/")
bam_to_variant_calls.pipeline(snakemake@input[[1]], directory)
