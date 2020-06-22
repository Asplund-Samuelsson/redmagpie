options(width=150)
library(tidyverse)

# Define infiles
supE_file = "results/Supplementary_Enrichment.tab"
supA_file = "results/Supplementary_ACE.tab"
supR_file = "results/Supplementary_Random_forest.tab"

# Load data
supE = read_tsv(supE_file)
supA = read_tsv(supA_file)
supR = read_tsv(supR_file)


supE = supE %>% filter(
  q < 0.05 &
  grepl("transport", Description, ignore.case=T) |
  grepl("transport", Name, ignore.case=T)
)

supA = supA %>% filter(
  Significant == 1 &
  grepl("transport", Description, ignore.case=T) |
  grepl("transport", Name, ignore.case=T)
)

supR = supR %>% filter(
  grepl("transport", Description, ignore.case=T) |
  grepl("transport", Name, ignore.case=T)
)

# Check number of significant features
supE %>% pull(Feature) %>% unique() %>% length()
supA %>% pull(Feature) %>% unique() %>% length()
supR %>% pull(Feature) %>% unique() %>% length()
c(
  supE %>% pull(Feature),
  supA %>% pull(Feature),
  supR %>% pull(Feature)
) %>% unique() %>% length()
