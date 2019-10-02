options(width=150)
library(tidyverse)

# Define infiles
bacdist_file = "/ssd/common/db/gtdb/r89/bac120.dist"
arcdist_file = "/ssd/common/db/gtdb/r89/ar122.dist"
meta_file = "data/gtdb_metadata.tab.gz"
cbbg_file = "data/rpp_genomes.txt"

# Load data
bacdist = read_delim(bacdist_file, " ", skip=1, col_names=F)
arcdist = read_delim(arcdist_file, " ", skip=1, col_names=F)
meta = read_tsv(meta_file)
cbbg = scan(cbbg_file, character())

# Alter column names for distance matrices
colnames(arcdist) = c("Genome", arcdist$X1)
colnames(bacdist) = c("Genome", bacdist$X1)

# Remove Genomes that are CBB-negative
arcdist = arcdist %>% filter(Genome %in% cbbg)
bacdist = bacdist %>% filter(Genome %in% cbbg)

# Remove cyanobacterial genomes
bacdist = bacdist %>% filter(
  !(Genome %in% (
    filter(meta, grepl(";p__Cyanobacteria;", gtdb_taxonomy)) %>% pull(accession)
  ))
)

# Gather into long format
arcdist = gather(arcdist, Relative, Distance, -Genome)
bacdist = gather(bacdist, Relative, Distance, -Genome)

# Add Completeness
arcdist = arcdist %>%
  # ...for Genome in focus
  inner_join(
    select(meta, accession, checkm_completeness) %>%
      rename(Genome = accession, CompG = checkm_completeness)
  ) %>%
  # ...for Relative of Genome in focus
  inner_join(
    select(meta, accession, checkm_completeness) %>%
      rename(Relative = accession, CompR = checkm_completeness)
  )

bacdist = bacdist %>%
  # ...for Genome in focus
  inner_join(
    select(meta, accession, checkm_completeness) %>%
      rename(Genome = accession, CompG = checkm_completeness)
  ) %>%
  # ...for Relative of Genome in focus
  inner_join(
    select(meta, accession, checkm_completeness) %>%
      rename(Relative = accession, CompR = checkm_completeness)
  )

# Add CBB status
arcdist = arcdist %>%
  mutate(
    CBBG = ifelse(Genome %in% cbbg, 1, 0),
    CBBR = ifelse(Relative %in% cbbg, 1, 0)
  )

bacdist = bacdist %>%
  mutate(
    CBBG = ifelse(Genome %in% cbbg, 1, 0),
    CBBR = ifelse(Relative %in% cbbg, 1, 0)
  )

# Compare only CBB Genomes to non-CBB Relatives
arcdist = arcdist %>% filter(CBBG == 1 & CBBR == 0)
bacdist = bacdist %>% filter(CBBG == 1 & CBBR == 0)

# Require completeness of at least 90%
arcdist = arcdist %>% filter(CompG >= 90 & CompR >= 90)
bacdist = bacdist %>% filter(CompG >= 90 & CompR >= 90)

# Remove completeness and CBB status as this is now implied
arcdist = arcdist %>% select(Genome, Relative, Distance)
bacdist = bacdist %>% select(Genome, Relative, Distance)

# Arrange the Genomes by Distance to their Relatives
arcdist = arcdist %>% arrange(Distance)
bacdist = bacdist %>% arrange(Distance)

# Create iterable copies of the distance matrices
arc = arcdist
bac = bacdist

# Iterate over the distance matrices removing Genomes and Relatives until empty
arc_examples = vector("list", length(unique(arc$Genome)))
i = 1
while (nrow(arc) > 0) {
  pair = arc[1,]
  arc_examples[[i]] = pair
  arc = filter(arc, Genome != pair$Genome & Relative != pair$Relative)
  i = i + 1
}
arc_examples = bind_rows(arc_examples)

bac_examples = vector("list", length(unique(bac$Genome)))
i = 1
while (nrow(bac) > 0) {
  pair = bac[1,]
  bac_examples[[i]] = pair
  bac = filter(bac, Genome != pair$Genome & Relative != pair$Relative)
  i = i + 1
}
bac_examples = bind_rows(bac_examples)

# Add Domain label
arc_examples = arc_examples %>% mutate(Domain = "Archaea")
bac_examples = bac_examples %>% mutate(Domain = "Bacteria")

# Combine tables and save as tab-delimited file
write_tsv(
  bind_rows(arc_examples, bac_examples),
  "intermediate/rpp_examples.no_cyano.tab"
)
