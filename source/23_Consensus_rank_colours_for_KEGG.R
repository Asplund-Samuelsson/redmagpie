options(width=150)
library(tidyverse)

# Define infiles
cons_file = "results/Supplementary_Consensus_rank.tab"
tran_file = "data/kegg_enzyme.old_new.tab"

# Load data
cons = read_tsv(cons_file)
tran = read_tsv(tran_file, col_names = c("EC", "NewEC"))

# Clean up EC annotations
cons = cons %>%
  # Keep only ECs
  filter(Feature_Type == "DeepEC") %>%
  rename(EC = Feature) %>%
  # Clean up EC
  mutate(EC = str_replace(EC, "EC", ""))

cons = cons %>%
  # Rename EC if it has been updated
  left_join(tran) %>%
  mutate(
    OldEC = ifelse(is.na(NewEC), NA, EC),
    EC = ifelse(is.na(NewEC), EC, NewEC)
  ) %>%
  select(-NewEC) %>%
  # Order by rank
  arrange(Rank)

# Calculate rank colours and write file for KEGG mapping
library(viridis)
library(grDevices)
pal = colorRamp(rev(viridis(64)))

# Function to determine if colour is "bright"
# https://trendct.org/2016/01/22/how-to-choose-a-label-color-to-contrast-with-background/
is_bright = function(x) {(x[1]*299 + x[2]*587 + x[3]*114) / 1000 > 123}

# Add colour to ranks
raco = cons %>%
  # Calculate log of rank
  mutate(logrank = log(Rank)) %>%
  # Select relevant data
  select(logrank, EC) %>%
  # Normalize logrank to range 0 to 1
  mutate(
    Rank = (logrank - min(logrank)) / (max(logrank) - min(logrank))
  ) %>%
  # Add colours
  mutate(
    Colour = unlist(lapply(
      Rank,
      function (a) {
        x = pal(a)
        rgb(x[1], x[2], x[3], maxColorValue=255)
      }
    )),
    TextColour = ifelse(
      unlist(lapply(Rank, function (a) {is_bright(pal(a))})),
      "#000000", "#ffffff"
    ),
    KEGGColour = paste(Colour, TextColour, sep=",")
  )

# Write colours in format expected by KEGG
write_delim(
  select(raco, EC, KEGGColour),
  "results/EC_count_features.KEGG_consensus_rank_colours.txt",
  delim=" ", col_names=F
)

# Create a log-transformed rank plot
gp = ggplot(cons, aes(x=Rank, y=1, colour=log(Rank)))
gp = gp + geom_point()
gp = gp + theme_bw()
gp = gp + theme(
  axis.text=element_text(colour="black"),
  axis.ticks=element_line(colour="black")
)
gp = gp + scale_colour_viridis_c(direction=-1)
gp = gp + scale_x_log10()
gp = gp + annotation_logticks()

ggsave(
  "results/EC_count_features.KEGG_consensus_rank_colours.logscale.pdf",
  gp, w=8, h=3
)
