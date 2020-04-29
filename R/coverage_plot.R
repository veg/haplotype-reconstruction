library(tidyverse)

df <- read_csv(snakemake@input[[1]]);
ggplot(df) +
  geom_line(aes(x=site, y=coverage, color=strain)) +
  ylim(0, quantile(df$coverage, .99, na.rm=TRUE))

ggsave(snakemake@output[[1]], width=5, height=3)
