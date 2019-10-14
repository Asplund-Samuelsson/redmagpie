options(width=150)
library(tidyverse)

# Define infiles
feim_file = "intermediate/EC_count_features.importance.tab.gz"
ecct_file = "intermediate/EC_count_features.X.csv.gz"
accn_file = "intermediate/EC_count_features.accession_ids.txt"
clss_file = "intermediate/EC_count_features.y.txt"
ecnm_file = "intermediate/EC_count_features.feature_names.txt"
tran_file = "data/kegg_enzyme.old_new.tab"
koec_file = "data/ko_ec.tab"

# Load data
feim = read_tsv(feim_file)
ecct = read_csv(ecct_file, col_names = scan(ecnm_file, character())) %>%
  mutate(
    Accession = scan(accn_file, character()),
    Class = scan(clss_file, integer())
  ) %>%
  gather(EC, Count, -Accession, -Class)
tran = read_tsv(tran_file, col_names = c("EC", "NewEC"))
koec = read_tsv(koec_file, col_names=c("KO", "EC"))

# Clean up KEGG EC annotations
ecct = ecct %>%
  # Remove features that were filtered out
  filter(EC %in% feim$Feature) %>%
  # Clean up EC
  mutate(EC = str_replace(EC, "EC", "")) %>%
  # Calculate mean Count per Class and EC
  group_by(Class, EC) %>%
  summarise(Count = mean(Count))

ecct = ecct %>%
  # Rename EC if it has been updated
  left_join(tran) %>%
  mutate(
    OldEC = ifelse(is.na(NewEC), NA, EC),
    EC = ifelse(is.na(NewEC), EC, NewEC)
  ) %>%
  select(-NewEC) %>%
  # Order by Count
  arrange(-Count)

ecct = ecct %>%
  # Add KO if available
  left_join(koec)

ecct = ecct %>%
  # Spread Class counts into columns
  ungroup() %>%
  mutate(Class = recode(Class, `1` = "Positive", `0` = "Negative")) %>%
  spread(Class, Count) %>%
  # Replace NA with zero
  mutate(
    Negative = replace_na(Negative, 0),
    Positive = replace_na(Positive, 0)
  )

ecct = ecct %>%
  # Calculate difference between Negative and Positive
  mutate(Difference = Positive - Negative) %>%
  # Sort by Difference
  arrange(-Difference)

# Calculate importance colours and write file for KEGG mapping
library(RColorBrewer)
library(grDevices)
pal = colorRamp(rev(brewer.pal(11, "PuOr")))

# Function to determine if colour is "bright"
# https://trendct.org/2016/01/22/how-to-choose-a-label-color-to-contrast-with-background/
is_bright = function(x) {(x[1]*299 + x[2]*587 + x[3]*114) / 1000 > 123}

# Add colour to importances
imco = ecct %>%
  # Remove if there is no KO assignment
  filter(!is.na(KO)) %>%
  # Select relevant data
  select(Difference, KO) %>%
  # Calculate mean for KOs with multiple ECs
  group_by(KO) %>%
  summarise(Difference = mean(Difference)) %>%
  ungroup() %>%
  # Normalize Difference to range 0 to 1 with 0.5 corresponding to 0
  mutate(
    xDifference =
    (Difference + max(abs(Difference)))/(max(Difference) + max(abs(Difference)))
  ) %>%
  # Add colours
  mutate(
    Colour = unlist(lapply(
      xDifference,
      function (a) {
        x = pal(a)
        rgb(x[1], x[2], x[3], maxColorValue=255)
      }
    )),
    TextColour = ifelse(
      unlist(lapply(xDifference, function (a) {is_bright(pal(a))})),
      "#000000", "#ffffff"
    ),
    KEGGColour = paste(Colour, TextColour, sep=",")
  )

# Write colours in format expected by KEGG
write_delim(
  select(imco, KO, KEGGColour),
  "results/EC_count_features.KEGG_difference_colours.txt",
  delim=" ", col_names=F
)
