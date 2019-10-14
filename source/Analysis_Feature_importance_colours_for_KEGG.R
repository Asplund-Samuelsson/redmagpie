options(width=150)
library(tidyverse)

# Define infiles
feim_file = "intermediate/EC_count_features.importance.tab.gz"
tran_file = "data/kegg_enzyme.old_new.tab"
koec_file = "data/ko_ec.tab"

# Load data
feim = read_tsv(feim_file)
tran = read_tsv(tran_file, col_names = c("EC", "NewEC"))
koec = read_tsv(koec_file, col_names=c("KO", "EC"))

# Clean up EC annotations
feim = feim %>%
  rename(EC = Feature) %>%
  # Clean up EC
  mutate(EC = str_replace(EC, "EC", "")) %>%
  # Calculate feature Importance mean
  group_by(EC) %>%
  summarise(Importance = mean(Importance), log10imp = log10(Importance))

feim = feim %>%
  # Rename EC if it has been updated
  left_join(tran) %>%
  mutate(
    OldEC = ifelse(is.na(NewEC), NA, EC),
    EC = ifelse(is.na(NewEC), EC, NewEC)
  ) %>%
  select(-NewEC) %>%
  # Order by importance
  arrange(-Importance)

feim = feim %>%
  # Add KO if available
  left_join(koec)

# Calculate importance colours and write file for KEGG mapping
library(viridis)
library(grDevices)
pal = colorRamp(viridis(64))

# Function to determine if colour is "bright"
# https://trendct.org/2016/01/22/how-to-choose-a-label-color-to-contrast-with-background/
is_bright = function(x) {(x[1]*299 + x[2]*587 + x[3]*114) / 1000 > 123}

# Add colour to importances
imco = feim %>%
  # Remove if there is no KO assignment or no importance
  filter(!is.na(KO), is.finite(log10imp)) %>%
  # Select relevant data
  select(log10imp, KO) %>%
  # Calculate mean for ECs with multiple KOs
  group_by(KO) %>%
  summarise(log10imp = mean(log10imp)) %>%
  ungroup() %>%
  # Normalize log10imp to range 0 to 1
  mutate(
    Importance = (log10imp - min(log10imp)) / (max(log10imp) - min(log10imp))
  ) %>%
  # Add colours
  mutate(
    Colour = unlist(lapply(
      Importance,
      function (a) {
        x = pal(a)
        rgb(x[1], x[2], x[3], maxColorValue=255)
      }
    )),
    TextColour = ifelse(
      unlist(lapply(Importance, function (a) {is_bright(pal(a))})),
      "#000000", "#ffffff"
    ),
    KEGGColour = paste(Colour, TextColour, sep=",")
  )

# Write colours in format expected by KEGG
write_delim(
  select(imco, KO, KEGGColour),
  "results/EC_count_features.KEGG_importance_colours.txt",
  delim=" ", col_names=F
)
