library(RegressHaplo)
split_path <- strsplit(snakemake@output[[1]], "/")[[1]]
directory <- paste(head(split_path, -1), collapse="/")
read_table_to_loci.pipeline(directory)
