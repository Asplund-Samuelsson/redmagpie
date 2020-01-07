options(width=150)
library(tidyverse)

# Define infiles
taxo_file = "intermediate/accession_taxonomy.tab"
exgn_file = "intermediate/example_genomes.tab"

# Load data
taxo = read_tsv(taxo_file)
exgn = read_tsv(exgn_file)

# Reformat example genomes table
exgn = exgn %>%
  select(-Distance) %>%
  gather(CBB, Accession, -Domain) %>%
  mutate(CBB = recode(CBB, "Genome" = "Positive", "Relative" = "Negative"))

# Determine top Orders
min_count = 30

orgs = taxo %>%
  # Select only example genomes, add domain and CBB status
  inner_join(exgn) %>%
  # Count genomes in each Order
  group_by(Group, Order) %>%
  # Use Group instead of Order if min_count is not achieved
  mutate(
    Organism0 = ifelse(
      length(Accession) >= min_count,
      Order,
      ifelse(
        Group %in% filter(., length(Accession) >= min_count)$Group,
        paste("Other", Group),
        Group
      )
    )
  ) %>%
  # Count genomes for each Organism
  group_by(Organism0) %>%
  # Set to "Other" unless min_count is achieved
  mutate(
    Organism = ifelse(
      length(Accession) >= min_count,
      Organism0,
      "Other"
    )
  )

# Count the number of genomes in each Organism set
orgc = orgs %>%
  group_by(CBB, Domain, Organism) %>%
  summarise(Count = length(Accession)) %>%
  # Add back Group
  left_join(
    orgs %>%
      ungroup() %>%
      select(Group, Organism) %>%
      mutate(Group = ifelse(Organism == "Other", Organism, Group)) %>%
      distinct()
    ) %>%
  group_by(Domain, Group, Organism) %>%
  summarise(Count = sum(Count)) %>%
  arrange(Domain, Group, -Count) %>%
  mutate(Index = 1:length(Organism))

# Decide colours
library(RColorBrewer)

clrs = tibble(
  Index = c(rep(rev(1:6), 5), 1),
  Group = c(
    rep("Halobacterota", 6),
    rep("Gammaproteobacteria", 6),
    rep("Actinobacteriota", 6),
    rep("Alphaproteobacteria", 6),
    rep("Firmicutes", 6),
    "Other"
  ),
  Colour = c(
    # c("#CCECE6", "#99D8C9", "#66C2A4", "#41AE76", "#238B45", "#005824")
    brewer.pal(7, "BuGn")[2:7],
    # c("#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177")
    brewer.pal(7, "RdPu")[2:7],
    # c("#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#4A1486")
    brewer.pal(7, "Purples")[2:7],
    # c("#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")
    brewer.pal(7, "Oranges")[2:7],
    # c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594")
    brewer.pal(7, "Blues")[2:7],
    "#F7F7F7"
  )
)

# Add colours
orgc = inner_join(orgc, clrs)
orgs = inner_join(orgs, select(ungroup(orgc), Organism, Colour) %>% distinct())

# Save Accession table with Organism and Colour definitions
write_tsv(orgs, "intermediate/accession_organism_colours.tab")

# Create plotting table
orgp = orgs %>%
  mutate(Group = ifelse(Organism == "Other", Organism, Group)) %>%
  group_by(CBB, Domain, Group, Organism, Colour) %>%
  summarise(Genomes = length(Accession)) %>%
  ungroup() %>%
  mutate(Organism = factor(
    Organism, levels=(
      select(., Group, Organism, Genomes) %>%
        group_by(Group, Organism) %>%
        summarise(Genomes = sum(Genomes)) %>%
        arrange(Genomes) %>%
        pull(Organism)
      )
    )
  ) %>%
  spread(CBB, Genomes) %>%
  gather(CBB, Genomes, -Domain, -Group, -Organism, -Colour) %>%
  mutate(
    Genomes = replace_na(Genomes, 0),
    Other = ifelse(Organism == "Other", "#c7c7c7", "#ffffff"),
    Label = ifelse(CBB == "Positive", "+", "-"),
    Label = ifelse(Genomes == 0, "", Label)
  )

gp = ggplot(
  orgp,
  aes(x=Organism, y=Genomes, group=CBB, fill=Colour, colour=Other, label=Label)
)
gp = gp + geom_col(position="dodge")
gp = gp + geom_text(
  colour="black", position=position_dodge(width=1), hjust=-0.5
)
gp = gp + theme_bw()
gp = gp + facet_grid(Domain+Group~., scales="free_y", space="free_y")
gp = gp + scale_fill_identity()
gp = gp + scale_colour_identity()
gp = gp + theme(
  strip.background=element_blank(),
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black"),
  strip.text.y = element_text(angle=0, hjust=0)
)
gp = gp + coord_flip()

ggsave(
  "results/taxonomic_distribution.pdf", gp, width=30, height=12, units="cm"
)
