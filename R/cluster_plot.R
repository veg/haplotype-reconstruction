library(tidyverse)
library(seqinr)

df <- read_csv(snakemake@input[[1]]);
df$true_cluster = ''
truth_fasta  <- read.fasta(snakemake@input[[2]])
for(i in 1:length(truth_fasta)) {
  cluster_label <- paste('c', i, sep='')
  current_header <- attr(truth_fasta[i], 'name')
  in_current_cluster <- grepl(current_header, df$sequence_id)
  df[in_current_cluster, 'true_cluster'] = cluster_label
}

ggplot(df)+geom_point(aes(d1,d2, color=true_cluster))+facet_wrap(~block, ncol=5)
ggsave(snakemake@output[[1]])

ggplot(df)+geom_point(aes(d1,d2, color=cluster))+facet_wrap(~block, ncol=5)
ggsave(snakemake@output[[2]])

ggplot(df)+geom_bar(aes(x=block, fill=cluster), position="dodge")
ggsave(snakemake@output[[3]])

