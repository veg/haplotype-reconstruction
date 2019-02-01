library(RegressHaplo)
split_path <- strsplit(snakemake@output[[1]], "/")[[1]]
directory <- paste(head(split_path, -1), collapse="/")
parameters_to_solutions.pipeline(directory, num_trials=50, rho_vals=c(.1,1,5,10,20))
