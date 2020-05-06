library(RegressHaplo)

# write start time
fileConn <- file(snakemake@output[[2]])
time <- as.integer(as.POSIXct(Sys.time()))
writeLines(as.character(time), fileConn)
close(fileConn)

# run pipeline
split_path <- strsplit(snakemake@output[[1]], "/")[[1]]
directory <- paste(head(split_path, -1), collapse="/")
full_pipeline(snakemake@input[[1]], directory, rho_vals=c(.1,1,5,10,20))

# write stop time
fileConn <- file(snakemake@output[[3]])
time <- as.integer(as.POSIXct(Sys.time()))
writeLines(as.character(time), fileConn)
close(fileConn)
