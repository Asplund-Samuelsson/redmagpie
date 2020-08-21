options(width=150)
library(tidyverse)

# Define infiles
ecim_file = "intermediate/EC_count_features.importance.tab.gz"
ecct_file = "intermediate/EC_count_features.X.csv.gz"
ecnm_file = "intermediate/EC_count_features.feature_names.txt"
accn_file = "intermediate/EC_count_features.accession_ids.txt"
clss_file = "intermediate/EC_count_features.y.txt"

# Load data
ecim = read_tsv(ecim_file)
ecct = read_csv(ecct_file, col_names = scan(ecnm_file, character())) %>%
  mutate(
    Accession = scan(accn_file, character()),
    Class = scan(clss_file, integer())
  ) %>%
  gather(EC, Count, -Accession, -Class)

# Clean up KEGG EC annotations
ecct = ecct %>%
  # Calculate mean Count per Class and EC
  group_by(Class, EC) %>%
  summarise(Count = mean(Count)) %>%
  # Spread Class counts into columns
  ungroup() %>%
  mutate(Class = recode(Class, `1` = "Positive", `0` = "Negative")) %>%
  spread(Class, Count) %>%
  # Replace NA with zero
  mutate(
    Negative = replace_na(Negative, 0),
    Positive = replace_na(Positive, 0)
  ) %>%
  # Calculate difference between Negative and Positive
  mutate(Difference = Positive - Negative) %>%
  # Add importance of each feature
  inner_join(
    ecim %>%
      rename(EC = Feature) %>%
      group_by(EC) %>%
      summarise(Importance = mean(Importance))
  )

# Correlate absolute EC count Difference to Importance
cor(abs(ecct$Difference), ecct$Importance)
