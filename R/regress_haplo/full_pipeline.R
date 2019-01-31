library(RegressHaplo)
split_path <- strsplit(snakemake@output[[1]], "/")[[1]]
directory <- paste(head(split_path, -1), collapse="/")
full_pipeline(snakemake@input[[1]], directory, rho_vals=c(.1,1,5,10,20))
