options(width=150)
library(tidyverse)

# Define infiles
feim_file = "intermediate/alignment_features.importance.tab.gz"
pred_file = "intermediate/alignment_features.prediction.tab.gz"
ecac_file = "intermediate/enzymes_to_align.tab"
tran_file = "data/kegg_enzyme.old_new.tab"
koec_file = "data/ko_ec.tab"

# Load data
feim = read_tsv(feim_file)
pred = read_tsv(pred_file)
ecac = read_tsv(ecac_file)
tran = read_tsv(tran_file, col_names = c("EC", "NewEC"))
koec = read_tsv(koec_file, col_names=c("KO", "EC"))

# Calculate accuracy
pred %>%
  group_by(Forest) %>%
  summarise(
    Accuracy = sum(Prediction == Class) / length(Class)
  ) %>%
  summarise(
    SD = sd(Accuracy),
    Accuracy = mean(Accuracy),
    CV = SD / Accuracy
  ) %>%
  select(Accuracy, SD, CV)

# Calculate average importance per feature
impo = feim %>%
  group_by(Feature) %>%
  summarise(
    SD = sd(Importance),
    Importance = mean(Importance),
    CV = SD / Importance
  ) %>%
  arrange(-Importance)

# Extract EC, Position, and Residue
str_split_i = function (x, p, i) {unlist(lapply(str_split(x, p), "[[", i))}

impo = impo %>%
  mutate(
    EC = str_split_i(Feature, "_", 2),
    Position = str_split_i(Feature, "_", 3),
    Residue = str_split_i(Feature, "_", 4),
    Residue = ifelse(Residue == "0", "-", Residue)
  )

# Summarise importances per Position and EC
imps = impo %>%
  mutate(Position = as.integer(Position)) %>%
  group_by(EC, Position) %>%
  summarise(Importance = sum(Importance)) %>%
  # Order EC by average Importance per Position
  ungroup() %>%
  mutate(
    EC = factor(
      EC,
      levels = group_by(., EC) %>%
        summarise(Ranking = sum(Importance) / length(Position)) %>%
        arrange(-Ranking) %>%
        pull(EC)
    )
  )

# Plot the importances
gp = ggplot(imps, aes(x=Position, y=Importance, fill=Importance))
gp = gp + geom_col()
gp = gp + facet_wrap(EC~., ncol=3, scales="free_x")
gp = gp + scale_fill_viridis_c()
gp = gp + theme_bw()
gp = gp + theme(
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black"),
  strip.background = element_blank(),
  panel.grid = element_blank()
)
ggsave(
  "results/alignment_features.importance_per_position.pdf",
  gp, w=768, h=432, units="mm"
)

# Count Features per EC for Features constituting minimal top including all ECs
for (i in 1:nrow(impo)) {
  if (length(unique(impo$EC[1:i])) == length(unique(impo$EC))) {break}
}

impt = head(impo, i)
impc = impt %>%
  group_by(EC) %>%
  summarise(Features = length(Feature)) %>%
  arrange(-Features)

# Create colours for KEGG based on top counts
impc = impc %>%
  # Rename EC if it has been updated
  left_join(tran) %>%
  mutate(
    OldEC = ifelse(is.na(NewEC), NA, EC),
    EC = ifelse(is.na(NewEC), EC, NewEC)
  ) %>%
  select(-NewEC)

impc = impc %>%
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
imco = impc %>%
  # Select relevant data
  select(Features, KO) %>%
  # Calculate mean for ECs with multiple KOs
  group_by(KO) %>%
  summarise(Features = log10(mean(Features))) %>%
  ungroup() %>%
  # Normalize Features count to range 0 to 1
  mutate(
    Intensity = (Features - min(Features)) / (max(Features) - min(Features))
  ) %>%
  # Add colours
  mutate(
    Colour = unlist(lapply(
      Intensity,
      function (a) {
        x = pal(a)
        rgb(x[1], x[2], x[3], maxColorValue=255)
      }
    )),
    TextColour = ifelse(
      unlist(lapply(Intensity, function (a) {is_bright(pal(a))})),
      "#000000", "#ffffff"
    ),
    KEGGColour = paste(Colour, TextColour, sep=",")
  )

# Write colours in format expected by KEGG
write_delim(
  select(imco, KO, KEGGColour),
  "results/alignment_features.KEGG_topcount_colours.txt",
  delim=" ", col_names=F
)
