options(width=150)
library(tidyverse)

# Define infiles
feco_file = "intermediate/EC_count_features.lr_coefficient.tab"
accn_file = "intermediate/EC_count_features.accession_ids.txt"
clss_file = "intermediate/EC_count_features.y.txt"
ecnm_file = "intermediate/EC_count_features.feature_names.txt"
tran_file = "data/kegg_enzyme.old_new.tab"
koec_file = "data/ko_ec.tab"

# Load data
feco = read_tsv(feco_file)
tran = read_tsv(tran_file, col_names = c("EC", "NewEC"))
koec = read_tsv(koec_file, col_names=c("KO", "EC"))

# Clean up KEGG EC annotations
feco = feco %>%
  # Clean up EC
  rename(EC = Feature) %>%
  mutate(EC = str_replace(EC, "EC", "")) %>%
  # Calculate mean Coefficient per EC
  group_by(EC) %>%
  summarise(Coefficient = mean(Coefficient)) %>%
  # Order by absolute Coefficient value
  arrange(-abs(Coefficient))

feco = feco %>%
  # Rename EC if it has been updated
  left_join(tran) %>%
  mutate(
    OldEC = ifelse(is.na(NewEC), NA, EC),
    EC = ifelse(is.na(NewEC), EC, NewEC)
  ) %>%
  select(-NewEC)

feco = feco %>%
  # Add KO if available
  left_join(koec)

# Calculate importance colours and write file for KEGG mapping
library(RColorBrewer)
library(grDevices)
pal = colorRamp(rev(brewer.pal(11, "PuOr")))

# Function to determine if colour is "bright"
# https://trendct.org/2016/01/22/how-to-choose-a-label-color-to-contrast-with-background/
is_bright = function(x) {(x[1]*299 + x[2]*587 + x[3]*114) / 1000 > 123}

# Add colour to importances
imco = feco %>%
  # Remove if there is no KO assignment
  filter(!is.na(KO)) %>%
  # Select relevant data
  select(Coefficient, KO) %>%
  # Calculate mean for KOs with multiple ECs
  group_by(KO) %>%
  summarise(Coefficient = mean(Coefficient)) %>%
  ungroup() %>%
  # Normalize Coefficient to range 0 to 1 with 0.5 corresponding to 0
  mutate(
    xCoefficient =
    (Coefficient + max(abs(Coefficient)))/(max(Coefficient) + max(abs(Coefficient)))
  ) %>%
  # Add colours
  mutate(
    Colour = unlist(lapply(
      xCoefficient,
      function (a) {
        x = pal(a)
        rgb(x[1], x[2], x[3], maxColorValue=255)
      }
    )),
    TextColour = ifelse(
      unlist(lapply(xCoefficient, function (a) {is_bright(pal(a))})),
      "#000000", "#ffffff"
    ),
    KEGGColour = paste(Colour, TextColour, sep=",")
  )

# Write colours in format expected by KEGG
write_delim(
  select(imco, KO, KEGGColour),
  "results/EC_count_features.KEGG_lr_coefficient_colours.txt",
  delim=" ", col_names=F
)
