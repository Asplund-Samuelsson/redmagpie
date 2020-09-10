options(width=150)
library(tidyverse)

# Load data for Table 1
tab1 = bind_rows(
  read_tsv("results/Table_Top_20_consensus_ECs.tab"),
  read_tsv("results/Table_Top_20_consensus_Pfams.tab")
)

# Format and save Table 1
write_tsv(
  tab1 %>%
    # Round values
    mutate(
      mean_Negative = round(mean_Negative, 2),
      mean_Positive = round(mean_Positive, 2),
      weighted_r = round(weighted_r, 2),
      Importance = round(Importance, 4)
    ) %>%
    # Change column names
    rename(
      `CBB+`  = mean_Positive,
      `CBB-`  = mean_Negative,
      `ACE r` = weighted_r,
      Imp. = Importance
    ) %>%
    # Reorder columns
    select(Rank, Feature, `CBB+`, `CBB-`, `ACE r`, Imp., Description),
  "results/Table_1_Top_consensus_ranks.tab"
)
