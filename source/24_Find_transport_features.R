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


# Check how many transporter features there are
bind_rows(
  select(supE, Feature, Name, Description),
  select(supA, Feature, Name, Description),
  select(supR, Feature, Name, Description)
) %>%
  distinct() %>%
  filter(
    grepl("transport|export|import|symport|antiport", Description, ignore.case=T) |
    grepl("transport|export|import|symport|antiport", Name, ignore.case=T) |
    grepl("ABC", Description) | grepl("ABC", Name)
  ) %>%
  pull(Feature) %>%
  unique() %>%
  length()


# Filter to significant

supE = supE %>% filter(
  q < 0.05 &
  grepl("transport|export|import|symport|antiport", Description, ignore.case=T) |
  grepl("transport|export|import|symport|antiport", Name, ignore.case=T) |
  grepl("ABC", Description) | grepl("ABC", Name)
)

supA = supA %>% filter(
  Significant == 1 &
  grepl("transport|export|import|symport|antiport", Description, ignore.case=T) |
  grepl("transport|export|import|symport|antiport", Name, ignore.case=T) |
  grepl("ABC", Description) | grepl("ABC", Name)
)

supR = supR %>% filter(
  grepl("transport|export|import|symport|antiport", Description, ignore.case=T) |
  grepl("transport|export|import|symport|antiport", Name, ignore.case=T) |
  grepl("ABC", Description) | grepl("ABC", Name)
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

# Load consensus table and check for transporters
supC = read_tsv("results/Supplementary_Consensus_rank.tab")

supC %>% filter(
  grepl("transport|export|import|symport|antiport", Description, ignore.case=T) |
  grepl("transport|export|import|symport|antiport", Name, ignore.case=T) |
  grepl("ABC", Description) | grepl("ABC", Name)
)

supC %>% filter(
  grepl("transport|export|import|symport|antiport", Description, ignore.case=T) |
  grepl("transport|export|import|symport|antiport", Name, ignore.case=T) |
  grepl("ABC", Description) | grepl("ABC", Name)
) %>% filter(Rank <= 1000) %>% pull(Feature) %>% unique() %>% length()
