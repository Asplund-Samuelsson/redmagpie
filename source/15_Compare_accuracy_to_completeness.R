options(width=150)
library(tidyverse)

# Define infiles
pred_file = "intermediate/EC_count_features.prediction.tab.gz"
meta_file = "data/gtdb_metadata.tab.gz"

# Load data
pred = read_tsv(pred_file)
meta = read_tsv(meta_file)

accu = pred %>%
  # Calculate accuracy
  group_by(Accession) %>%
  summarise(Accuracy = sum(Prediction == Class) / length(Prediction)) %>%
  # Add completeness
  inner_join(
    select(meta, accession, checkm_completeness) %>%
    rename(Accession = accession, Completeness = checkm_completeness)
  ) %>%
  mutate(Completeness = Completeness / 100) %>%
  # Add class
  inner_join(select(pred, Accession, Class) %>% distinct()) %>%
  mutate(Class = recode(Class, `1` = '+', `0` = '-'))

# Plot it
PuOr = c("#f1a340", "#998ec3")

gp = ggplot(accu, aes(x=Completeness, y=Accuracy, shape=Class, colour=Class))
gp = gp + geom_point(size=4, alpha=0.5)
gp = gp + scale_shape_identity()
gp = gp + theme_bw()
gp = gp + theme(
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black"),
  legend.position = "none",
  aspect.ratio = 1,
  panel.grid = element_blank()
)
gp = gp + scale_colour_manual(values=PuOr)

gp1 = gp

gp = ggplot(accu, aes(x=Completeness, fill=Class))
gp = gp + geom_density(alpha=0.5, colour=NA)
gp = gp + theme_bw()
gp = gp + theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none",
  panel.grid = element_blank()
)
gp = gp + scale_fill_manual(values=PuOr)

gp2 = gp

gp = ggplot(accu, aes(x=Accuracy, fill=Class))
gp = gp + geom_density(alpha=0.5, colour=NA)
gp = gp + theme_bw()
gp = gp + theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none",
  panel.grid = element_blank()
)
gp = gp + coord_flip()
gp = gp + scale_fill_manual(values=PuOr)

gp3 = gp

library(egg)

png(
  "results/EC_count_features.acc_vs_comp.png",
  h=120, w=120, units="mm", res=300
)
ggarrange(
  gp2, ggplot() + theme_minimal(), gp1, gp3,
  nrow=2, ncol=2, widths=c(4,1), heights=c(1,4)
)
garbage = dev.off()
