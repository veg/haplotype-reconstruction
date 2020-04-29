library(tidyverse)

df <- read_csv(snakemake@input[[1]]);
ggplot(df %>% filter(b==.6)) +
  geom_line(aes(x=site, y=coverage,color=strain)) +
  ylim(0, quantile(df$coverage, .99, na.rm=TRUE)) +
  facet_grid(m ~ q, labeller = label_both)

ggsave(snakemake@output[[1]], width=12, height=7)
