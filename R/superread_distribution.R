library(tidyverse)

df <- read_csv(snakemake@input[[1]]);

ggplot(df) +
  geom_histogram(mapping=aes(x=weight, fill=composition))
ggsave(snakemake@output[[1]], width=8, height=5)

ggplot(df) +
  geom_point(mapping=aes(x=weight, y=vacs_length, color=composition))
ggsave(snakemake@output[[2]], width=8, height=5)
