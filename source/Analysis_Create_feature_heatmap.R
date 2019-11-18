options(width=150)
library(tidyverse)

# Load metadata
exgn = read_tsv("intermediate/example_genomes.tab")
taxo = read_tsv("intermediate/accession_taxonomy.tab")
comp = read_tsv("data/gtdb_metadata.tab.gz") %>%
  select(accession, checkm_completeness) %>%
  rename(Accession = accession, Completeness = checkm_completeness) %>%
  filter(Accession %in% c(exgn$Genome, exgn$Relative)) %>%
  mutate(Completeness = Completeness / 100)

# Load features
dpec = read_csv(
  "intermediate/EC_count_features.X.csv.gz",
  col_names = scan(
    "intermediate/EC_count_features.feature_names.txt", character()
    )
  ) %>%
  # Add Accession IDs
  mutate(
    Accession = scan(
      "intermediate/EC_count_features.accession_ids.txt", character()
    )
  ) %>%
  # Gather into long format
  gather(Feature, Count, -Accession) %>%
  # Store origin of Data
  mutate(Data = "EC") %>%
  # Create a normalized Value per Feature
  group_by(Feature) %>%
  mutate(Value = (Count - min(Count))/(max(Count) - min(Count)))

pfam = read_csv(
    "intermediate/pfam_features.X.csv.gz",
    col_names = scan(
      "intermediate/pfam_features.feature_names.txt", character()
    )
  ) %>%
  # Add Accession IDs
  mutate(
    Accession = scan("intermediate/pfam_features.accession_ids.txt", character())
  ) %>%
  # Gather into long format
  gather(Feature, Count, -Accession) %>%
  mutate(Data = "Pfam") %>%
  # Create a normalized Value per Feature
  group_by(Feature) %>%
  mutate(Value = (Count - min(Count))/(max(Count) - min(Count)))

ecal = read_csv(
  "intermediate/alignment_features.X.csv.gz",
  col_names = scan(
    "intermediate/alignment_features.feature_names.txt", character()
    )
  ) %>%
  # Add Accession IDs
  mutate(
    Accession = scan(
      "intermediate/alignment_features.accession_ids.txt", character()
    )
  ) %>%
  # Gather into long format
  gather(Feature, Value, -Accession) %>%
  # Store origin of Data
  mutate(Data = "Alignment")

# Load prediction and importance data
alim = read_tsv("intermediate/alignment_features.importance.tab.gz")
alpr = read_tsv("intermediate/alignment_features.prediction.tab.gz")
ecim = read_tsv("intermediate/EC_count_features.importance.tab.gz")
ecpr = read_tsv("intermediate/EC_count_features.prediction.tab.gz")
pfim = read_tsv("intermediate/pfam_features.importance.tab.gz")
pfpr = read_tsv("intermediate/pfam_features.prediction.tab.gz")

# Calculate average Importance, select top 300, and order Feature
alim = alim %>%
  group_by(Feature) %>%
  summarise(Importance = mean(Importance)) %>%
  top_n(200, Importance) %>%
  arrange(Importance)
ecal = ecal %>%
  inner_join(alim)

ecim = ecim %>%
  group_by(Feature) %>%
  summarise(Importance = mean(Importance)) %>%
  top_n(200, Importance) %>%
  arrange(Importance)
dpec = dpec %>%
  inner_join(ecim) %>%
  ungroup()

pfim = pfim %>%
  group_by(Feature) %>%
  summarise(Importance = mean(Importance)) %>%
  top_n(200, Importance) %>%
  arrange(Importance)
pfam = pfam %>%
  inner_join(pfim) %>%
  ungroup()

# Calculate average Accuracy for each Accession, and order Accession
alac = alpr %>%
  group_by(Accession) %>%
  summarise(Accuracy = sum(Prediction == Class) / length(Class)) %>%
  mutate(Data = "Alignment")

ecac = ecpr %>%
  group_by(Accession) %>%
  summarise(Accuracy = sum(Prediction == Class) / length(Class)) %>%
  mutate(Data = "EC")

pfac = pfpr %>%
  group_by(Accession) %>%
  summarise(Accuracy = sum(Prediction == Class) / length(Class)) %>%
  mutate(Data = "Pfam")

accu = bind_rows(alac, ecac, pfac)
avac = accu %>%
  group_by(Accession) %>%
  summarise(Accuracy = mean(Accuracy)) %>%
  inner_join(select(taxo, Accession, Group)) %>%
  inner_join(comp) %>%
  mutate(
    Class = ifelse(Accession %in% exgn$Genome, "CBB-positive", "CBB-negative"),
    SortValue = ifelse(Class == "CBB-positive", Accuracy, -Accuracy),
    SortValue2 = ifelse(Class == "CBB-positive", Completeness, -Completeness)
  ) %>%
  arrange(desc(Class), -SortValue, Group, -SortValue2)
accu = accu %>%
  mutate(Accession = factor(Accession, levels=avac$Accession))
comp = comp %>%
  mutate(Accession = factor(Accession, levels=avac$Accession))

# Combine Values and order Accession according to Class, Accuracy, and Group
fval = bind_rows(select(dpec, -Count), select(pfam, -Count), ecal) %>%
  inner_join(select(avac, -Accuracy, -SortValue)) %>%
  mutate(
    Feature = factor(
      Feature,
      levels=c(pfim$Feature, ecim$Feature, alim$Feature)
    ),
    Accession = factor(Accession, levels=avac$Accession),
    Data = factor(Data, levels = c("Pfam", "EC", "Alignment"))
  )

# Create average Accuracy plot
gp = ggplot(accu, aes(x=Accession, y=Accuracy, shape=Data))
gp = gp + geom_point(alpha=0.2, size=0.8)
gp = gp + annotate(
  geom="text", x=10, y=0.5, label="Average accuracy",
  colour="black", size = 14, hjust=0
)
gp = gp + scale_y_continuous(breaks=c(0,0.5,1))
gp = gp + theme_bw()
gp = gp + theme(
  axis.text.y = element_text(size=12),
  axis.text.x = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  legend.direction="horizontal",
  legend.position=c(0.16,0.5),
  legend.background = element_blank(),
  panel.border = element_blank(),
  legend.text = element_text(size=18),
  legend.title = element_blank()
)
gpac = gp

# Create Completeness plot
gp = ggplot(comp, aes(x=Accession, y=Completeness))
gp = gp + geom_point(alpha=0.2, size=0.8)
gp = gp + annotate(
  geom="text", x=10, y=0.7, label="Completeness",
  colour="black", size = 14, hjust=0
)
gp = gp + scale_y_continuous(breaks=c(0.4, 0.7, 1.0))
gp = gp + theme_bw()
gp = gp + theme(
  axis.text.y = element_text(size=12),
  axis.text.x = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  legend.text = element_text(size=18),
  legend.title = element_blank()
)
gpcm = gp

# Create phylogenetic Group plot
grtb = select(fval, Accession, Group) %>% distinct()

# cltb = grtb %>%
#   # Select top 6 groups
#   group_by(Group) %>%
#   summarise(Count = length(Group)) %>%
#   arrange(-Count) %>%
#   top_n(6, Count) %>%
#   # Remove Cyanobacteria
#   filter(Group != "Cyanobacteria") %>%
#   # Add colour for remaining top 5 groups
#   mutate(
#     Colour = rev(c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3"))
#   ) %>%
#   # Add colour for Cyanobacteria and Other
#   select(-Count) %>%
#   bind_rows(
#     tibble(
#       Group = c("Cyanobacteria", "Other"),
#       Colour = c("#35978f", "#f5f5f5")
#     )
#   )

cltb = grtb %>%
  # Select top 5 groups
  group_by(Group) %>%
  summarise(Count = length(Group)) %>%
  top_n(5, Count) %>%
  arrange(-Count) %>%
  # Add colour
  mutate(
   Colour = rev(c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3"))
  ) %>%
  # Add colour for Other
  select(-Count) %>%
  bind_rows(tibble(Group = "Other", Colour = "#f5f5f5"))

grtb = grtb %>%
  # Rename Group if rare
  mutate(Group = ifelse(Group %in% cltb$Group, Group, "Other")) %>%
  # Add Colour of Group
  inner_join(cltb)

gp = ggplot(grtb, aes(x=Accession, fill=Colour, y=1))
gp = gp + geom_tile()
gp = gp + scale_fill_identity(
  labels = cltb$Group, breaks=cltb$Colour, guide="legend"
)
gp = gp + annotate(
  geom="text", x=10, y=1, label="Organism", colour="white", size = 14, hjust=0
)
gp = gp + theme_bw()
gp = gp + guides(fill = guide_legend(nrow=1))
gp = gp + theme(
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  legend.position = "top",
  legend.direction = "horizontal",
  panel.border = element_blank(),
  legend.text = element_text(size=18),
  legend.title = element_blank()
)
gpgr = gp

# Create Class plot
gp = ggplot(
  select(fval, Accession, Class) %>% distinct(),
  aes(x=Accession, fill=Class, y=1)
)
gp = gp + geom_tile()
gp = gp + scale_fill_manual(values=c("white", "black"), guide=F)
gp = gp + theme_bw()
gp = gp + annotate(
  geom="text", x=10, y=1, label="CBB-positive",
  colour="white", size = 14, hjust=0
)
gp = gp + annotate(
  geom="text", x=2021, y=1, label="CBB-negative",
  colour="black", size = 14, hjust=0
)
gp = gp + theme(
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_blank()
)
gpcl = gp

# Create Importance plot
gp = ggplot(
  select(fval, Feature, Data, Importance) %>% distinct(),
  aes(y=Importance, x=Feature)
)
gp = gp + geom_col()
gp = gp + theme_bw()
gp = gp + coord_flip()
gp = gp + facet_grid(Data~., scales="free_y")
gp = gp + scale_y_reverse()
gp = gp + theme(
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  strip.background = element_blank(),
  strip.text = element_blank(),
  panel.border = element_blank()
)
gpim = gp

# Create Value plot
gp = ggplot(fval, aes(x=Accession, y=Feature, fill=Value))
gp = gp + geom_tile()
gp = gp + theme_bw()
gp = gp + facet_grid(Data~., scales="free_y")
gp = gp + scale_fill_viridis_c()
gp = gp + theme(
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  strip.background = element_blank(),
  panel.border = element_blank(),
  strip.text = element_text(size=18),
  legend.title = element_blank(),
  legend.text = element_text(size=18)
)
gpvl = gp

# Arrange all plots
library(egg)

gp0 = ggplot() + theme_minimal()

png(
  "results/Feature_heatmap.png",
  h=700, w=1200, units="mm", res=300
)
ggarrange(
  gp0, gp0, gp0, gp0,  gpgr,
  gp0, gp0, gp0, gp0,  gpcl,
  gp0, gp0, gp0, gp0,  gpac,
  gp0, gp0, gp0, gp0,  gpcm,
  gp0, gp0, gp0, gpim, gpvl,
  ncol=5, widths=c(0,0,0,1,30), heights=c(1,1,1,1,30)
)
garbage = dev.off()
