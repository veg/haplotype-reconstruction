library(RegressHaplo)
split_path <- strsplit(snakemake@output[[1]], "/")[[1]]
directory <- paste(head(split_path, -1), collapse="/")
variant_calls_to_read_table.pipeline(snakemake@input[[1]], directory)
