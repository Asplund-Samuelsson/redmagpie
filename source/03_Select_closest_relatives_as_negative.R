options(width=150)
library(tidyverse)

# Define infiles
bacdist_file = "intermediate/bacteria.dist"
arcdist_file = "intermediate/archaea.dist"
posg_file = "data/positive_genomes.txt"

# Load data
bacdist = read_delim(bacdist_file, " ", skip=1, col_names=F)
arcdist = read_delim(arcdist_file, " ", skip=1, col_names=F)
posg = scan(posg_file, character())

# Alter column names for distance matrices
colnames(arcdist) = c("Genome", arcdist$X1)
colnames(bacdist) = c("Genome", bacdist$X1)

# Keep only Positive genomes
arcdist = arcdist %>% filter(Genome %in% posg)
bacdist = bacdist %>% filter(Genome %in% posg)

# Gather into long format
arcdist = gather(arcdist, Relative, Distance, -Genome)
bacdist = gather(bacdist, Relative, Distance, -Genome)

# Add genome status
arcdist = arcdist %>%
  mutate(
    Positive = ifelse(Genome %in% posg, 1, 0),
    Negative = ifelse(Relative %in% posg, 1, 0)
  )

bacdist = bacdist %>%
  mutate(
    Positive = ifelse(Genome %in% posg, 1, 0),
    Negative = ifelse(Relative %in% posg, 1, 0)
  )

# Compare only Positive genomes to Negative Relatives
arcdist = arcdist %>% filter(Positive == 1 & Negative == 0)
bacdist = bacdist %>% filter(Positive == 1 & Negative == 0)

# Remove genome status (Positive or Negative) as this is now implied
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
  "intermediate/example_genomes.tab"
)
