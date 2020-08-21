options(width=150)
library(tidyverse)
library(ggrepel)

# Define infiles
supe_file = "results/Supplementary_Enrichment.tab"
supa_file = "results/Supplementary_ACE.tab"
supr_file = "results/Supplementary_Random_forest.tab"
cons_file = "results/Supplementary_Consensus_rank.tab"
corr_file = "results/method_correlation.tab"
rubd_file = "results/Supplementary_Gene_proximity.tab"

# Load data
supe = read_tsv(supe_file)
supa = read_tsv(supa_file)
supr = read_tsv(supr_file)
cons = read_tsv(cons_file)
corr = read_tsv(corr_file)
rubd = read_tsv(rubd_file)

# Combine ranks
rnks = bind_rows(
  supe %>%
    select(Rank, Feature_Type, Feature) %>%
    distinct() %>%
    mutate(Method = "Enrichment"),
  supa %>%
    select(Rank, Feature_Type, Feature) %>%
    distinct() %>%
    mutate(Method = "ACE"),
  supr %>%
    select(Rank, Feature_Type, Feature) %>%
    distinct() %>%
    mutate(Method = "Random forest")
)

# Add consensus rank
rnks = rnks %>%
  inner_join(
    cons %>%
      select(Rank, Feature) %>%
      distinct() %>%
      rename(Consensus = Rank)
  ) %>%
  left_join(
    rubd %>%
      filter(Strand == "Same", cFeature == "Rubisco") %>%
      select(Feature, medD, Count) %>%
      distinct() %>%
      filter(Count > 200)
  ) %>%
  left_join(
    cons %>%
      select(Feature, Name) %>%
      mutate(
        Label = recode(
          Feature,
          "EC4.1.2.13" = "ALD",
          "PF00316.20" = "FBPase",
          "PF01116.20" = "ALD",
          "PF17866.1" = "CbbX",
          "EC3.1.3.11" = "FBPase",
          "EC2.2.1.6" = "ALS",
          "PF10531.9" = "SLBB",
          "PF04172.16" = "LrgB",
          "PF08406.10" = "CbbQ",
          "EC2.2.1.1" = "TKT",
          "EC5.1.3.1" = "RPE",
          "PF00834.19" = "RPE",
          "PF00936.19" = "BMC",
          "PF17042.5" = "NBD",
          "PF07005.11" = "SBD",
          "PF01058.22" = "NDH"
        ),
        Label = ifelse(Label == Feature, NA, Label)
      )
  ) %>%
  mutate(
    Method = factor(Method, levels=c("Enrichment", "ACE", "Random forest"))
  )

# Plot it
gp = ggplot(rnks, aes(x=Rank, y=Consensus, colour=medD, shape=Feature_Type))
gp = gp + geom_point(data=filter(rnks, is.na(medD)), size=1)
gp = gp + geom_point(data=filter(rnks, !is.na(medD)), size=1)
gp = gp + geom_text_repel(
  mapping=aes(label=Label), size=2.8, box.padding=0.4, segment.size=0.3,
  point.padding=0.01
)
gp = gp + theme_bw()
gp = gp + scale_colour_gradientn(
  na.value="#bcbddc",
  colours=c(
    "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404", "#662506"
  ),
  values = c(1,0.83,0.67,0.5,0.33,0.17,0),
  breaks=c(1,40,160,400,800),
  trans="sqrt"
)
gp = gp + facet_grid(~Method, scales="free_x", space="free")
gp = gp + annotation_logticks(
  size=0.3, short=unit(0.07, "cm"), mid=unit(0.14, "cm"), long=unit(0.21, "cm")
)
gp = gp + theme(
  strip.background = element_blank(),
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black", size=0.3),
  aspect.ratio=1,
  legend.position="bottom",
  legend.box="horizontal",
  legend.title = element_text(size=10),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(size=0.3)
)
gp = gp + labs(
  shape="", colour="Median number of genes to Rubisco (n > 200)",
  y="Consensus rank", x="Method rank"
)
gp = gp + scale_x_log10(limits=c(0.5, NA)) + scale_y_log10(limits=c(0.5, NA))

ggsave("results/method_feature_rank_comparison_2.pdf", gp, h=10, w=18, units="cm")
