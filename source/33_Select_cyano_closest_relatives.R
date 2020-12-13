options(width=100)
library(tidyverse)

# Define infiles
bacdist_file = "intermediate/bacteria.dist"
posg_file = "intermediate/cyano_with_CBB.txt"

# Load data
bacdist = read_delim(bacdist_file, " ", skip=1, col_names=F)
posg = scan(posg_file, character())

# Alter column names for distance matrices
colnames(bacdist) = c("Genome", bacdist$X1)

# Keep only Positive genomes
bacdist = bacdist %>% filter(Genome %in% posg)

# Gather into long format
bacdist = gather(bacdist, Relative, Distance, -Genome)

# Add genome status
bacdist = bacdist %>%
  mutate(
    Positive = ifelse(Genome %in% posg, 1, 0),
    Negative = ifelse(Relative %in% posg, 1, 0)
  )

# Compare only Positive genomes to Negative Relatives
bacdist = bacdist %>% filter(Positive == 1 & Negative == 0)

# Remove genome status (Positive or Negative) as this is now implied
bacdist = bacdist %>% select(Genome, Relative, Distance)

# Arrange the Genomes by Distance to their Relatives
bacdist = bacdist %>% arrange(Distance)

# Create iterable copies of the distance matrices
bac = bacdist

# Iterate over the distance matrices removing Genomes and Relatives until empty
bac_examples = vector("list", length(unique(bac$Genome)))
i = 1
while (nrow(bac) > 0) {
  pair = bac[1,]
  bac_examples[[i]] = pair
  bac = filter(bac, Genome != pair$Genome & Relative != pair$Relative)
  i = i + 1
}
bac_examples = bind_rows(bac_examples)

# Save as tab-delimited file
write_tsv(bac_examples, "intermediate/cyano_relatives.tab")
