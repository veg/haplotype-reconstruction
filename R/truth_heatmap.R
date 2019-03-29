library(tidyverse)


df <- read_csv(snakemake@input[[1]]);
dd_min <- 0
dd_max <- 5 
dd_mid = (dd_min + dd_max) / 2
ggplot(df) + geom_tile(aes(first_record, second_record, fill=distance)) +
  theme(
    axis.text.x = element_text(angle = 30, hjust=1),
    plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  scale_fill_gradient2(
    low="blue",
    high ="red",
    limits=c(dd_min, dd_max),
    midpoint=dd_mid
  ) +
  coord_equal()
ggsave(snakemake@output[[1]], width=8, height=5)

