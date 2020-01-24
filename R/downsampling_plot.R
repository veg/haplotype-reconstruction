library(tidyverse)
library(jsonlite)

downsample_values = 90:99
all_data <- lapply(snakemake@input, function(filename) { fromJSON(filename) })
plot_data <- tibble(
  downsampling=downsample_values,
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
  scale_x_continuous(breaks=downsample_values)
ggsave(snakemake@output[[1]], width=16, height=8)
