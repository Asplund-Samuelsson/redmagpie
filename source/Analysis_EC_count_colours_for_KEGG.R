options(width=150)
library(tidyverse)

# Define infiles
feim_file = "intermediate/EC_count_features.importance.tab.gz"
enri_file = "results/Supplementary_Enrichment.tab"
tran_file = "data/kegg_enzyme.old_new.tab"

# Load data
feim = read_tsv(feim_file)
enri = read_tsv(enri_file)
tran = read_tsv(tran_file, col_names = c("EC", "NewEC"))

# Clean up KEGG EC annotations
ecct = enri %>%
  rename(EC = Feature) %>%
  # Remove features that were filtered out
  filter(EC %in% feim$Feature) %>%
  # Clean up EC
  mutate(EC = str_replace(EC, "EC", ""))

ecct = ecct %>%
  # Select relevant data
  select(EC, mean_Negative, mean_Positive, q) %>%
  rename(Negative = mean_Negative, Positive = mean_Positive) %>%
  distinct() %>%
  # Rename EC if it has been updated
  left_join(tran) %>%
  mutate(
    OldEC = ifelse(is.na(NewEC), NA, EC),
    EC = ifelse(is.na(NewEC), EC, NewEC)
  ) %>%
  select(-NewEC)

ecct = ecct %>%
  # Calculate difference between Negative and Positive
  mutate(Difference = Positive - Negative) %>%
  # Sort by Difference
  arrange(-Difference)

# Add discrete colour to importances
dico = ecct %>%
  # Select relevant data
  select(Difference, EC, q) %>%
  # Add colours
  mutate(
    Colour = case_when(
      Difference > 0 & q < 0.05 ~ "#f1a340",
      Difference < 0 & q < 0.05 ~ "#998ec3",
      TRUE ~ "#f7f7f7"
    ),
    TextColour = ifelse(
      Colour != "#998ec3",
      "#000000", "#ffffff"
    ),
    KEGGColour = paste(Colour, TextColour, sep=",")
  )

# Write colours in format expected by KEGG
write_delim(
  select(dico, EC, KEGGColour),
  "results/EC_count_features.KEGG_difference_colours.discrete.txt",
  delim=" ", col_names=F
)
