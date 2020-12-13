options(width=120)
library(tidyverse)

# Define infiles
cyan_file = "intermediate/cyano_relatives.tab"
othr_file = "intermediate/example_genomes.tab"
taxo_file = "intermediate/accession_taxonomy.tab"

# Load data
cyan = read_tsv(cyan_file)
othr = read_tsv(othr_file)
taxo = read_tsv(taxo_file)

# Combine data and add group of Genome and Relative
comp = cyan %>%
  mutate(Domain = "Cyanobacteria") %>%
  bind_rows(othr) %>%
  left_join(
    taxo %>%
      select(Accession, Group) %>%
      rename(Genome = Accession, Genome_organism = Group)
  ) %>%
  left_join(
    taxo %>%
      select(Accession, Group) %>%
      rename(Relative = Accession, Relative_organism = Group)
  )

# Compute and compare average distances
diff = comp %>%
  mutate(Type = ifelse(Domain == "Cyanobacteria", "Cyanobacteria", "Other")) %>%
  group_by(Type) %>%
  summarise(
    Distances = list(Distance),
    Distance = median(Distance)
  )

# Perform Wilcox test
diff$p = wilcox.test(
  unlist(filter(diff, Type == "Cyanobacteria")$Distances),
  unlist(filter(diff, Type == "Other")$Distances)
)$p.value

# Save final result
write_tsv(
  diff %>% select(-Distances) %>% spread(Type, Distance),
  "results/cyano_vs_others_relative_dist.tab"
)
