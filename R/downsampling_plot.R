library(tidyverse)
library(jsonlite)

all_data <- lapply(snakemake@input, function(filename) { fromJSON(filename) })
plot_data <- tibble(
  downsampling=c(0, 20, 40, 60, 80),
  precision=as.numeric(lapply(all_data, function(datum) { datum$precision })),
  recall=as.numeric(lapply(all_data, function(datum) { datum$recall }))
) %>% gather(precision, recall, key="statistic", value="value")

plot_data

ggplot(plot_data) + 
  geom_bar(
    aes(x=downsampling, y=value, fill=statistic),
    stat="identity",
    position="dodge"
  ) + 
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))
ggsave(snakemake@output[[1]], width=8, height=8)
