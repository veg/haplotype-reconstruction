library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)

bam <- scanBam(snakemake@input[[1]])[[1]]
df <- tibble(
  insertions=sapply(bam$cigar, function(c) {
    ifelse(is.na(c), 0, sum(explodeCigarOpLengths(c)[[1]][explodeCigarOps(c)[[1]] == "I"]))
  }),
  position=bam$pos
) %>% 
  drop_na() %>%
  group_by(position) %>%
  summarise(insertions=sum(insertions)) %>%
  mutate(insertionsp1 = insertions+1)
ggplot(df %>% dplyr::slice(1:20)) +
  geom_bar(aes(x=position, y=insertionsp1), stat='identity') + 
  #scale_y_log10() +
  ylab('log(insertions)')
ggsave(snakemake@output[[1]], width=8, height=5)

