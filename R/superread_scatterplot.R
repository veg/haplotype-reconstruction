library(tidyverse)

df <- read_csv(snakemake@input[[1]])
ggplot(df, aes(weight, vacs_length)) +
  geom_point()
ggsave(snakemake@output[[1]], width=8, height=5)
