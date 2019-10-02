options(width=150)
library(tidyverse)

# Define infiles
ecal_file = "intermediate/map01200.no_cyano.ec_orf_ali.tab.gz"
cbbx_file = "intermediate/rpp_examples.no_cyano.tab"
meta_file = "data/gtdb_metadata.tab.gz"

# Load data
ecal = read_tsv(ecal_file, col_names = c("EC", "Sequence", "Alignment"))
cbbx = read_tsv(cbbx_file)
meta = read_tsv(meta_file)

# Prepare taxonomic information
taxo = meta %>%
  # Select relevant data
  select(accession, gtdb_taxonomy) %>%
  # Extract Domain, Phylum, and Class
  mutate(
    Domain = str_replace(
      unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 1)), "d__", ""
    ),
    Phylum = str_replace(
      unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 2)), "p__", ""
    ),
    Class = str_replace(
      unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 3)), "c__", ""
    )
  ) %>%
  # Drop the nasty taxonomy string and rename accession
  select(-gtdb_taxonomy) %>%
  rename(Accession = accession) %>%
  # Create Organism classification
  mutate(Organism = ifelse(Phylum == "Proteobacteria", Class, Phylum)) %>%
  select(Accession, Organism)

ecal = ecal %>%
  # Split alignment into amino acids
  mutate(Alignment = str_split(Alignment, "")) %>%
  unnest(Alignment)

ecal = ecal %>%
  # Rename Alignment to Residue
  rename(Residue = Alignment) %>%
  # Add position
  group_by(EC, Sequence) %>%
  mutate(Position = 1:length(Residue))

resi = ecal %>%
  # Count number of each residue per position and EC
  ungroup() %>%
  group_by(EC, Position, Residue) %>%
  summarise(Count = length(Residue)) %>%
  # Calculate fraction per position
  inner_join(
    ungroup(ecal) %>%
    select(EC, Sequence) %>%
    distinct() %>%
    group_by(EC) %>%
    summarise(Sequences = length(Sequence))
  ) %>%
  mutate(Fraction = Count / Sequences)

resf = resi %>%
  # Keep only positions with less than 5% gaps
  anti_join(
    filter(resi, Residue == "-" & Fraction >= 0.05) %>% select(EC, Position)
  )

preX = ecal %>%
  # Reduce input data to positions of interest
  inner_join(resf) %>%
  # Extract Accession
  mutate(
    Accession = unlist(lapply(str_split(Sequence, fixed("|")), "[[", 1))
  )

# Obtain "consensus" sequence
cons = resf %>%
  # Pick the most common residues per position
  group_by(EC, Position) %>%
  filter(Count == max(Count)) %>%
  # Count the number of top residues and calculate a value for the X matrix
  mutate(Value = 1/length(Residue))

preX = preX %>%
  # Normalize value for the X matrix to number of sequences per accession
  ungroup() %>%
  group_by(EC, Accession, Position) %>%
  mutate(Value = 1/length(Residue))

preX = preX %>%
  # Add consensus sequence for Accession without EC
  bind_rows(
    inner_join(
      mutate(cons, Sequence = "MOCK"),
      anti_join(
        crossing(EC=unique(preX$EC), Accession=c(cbbx$Genome, cbbx$Relative)),
        preX %>% ungroup() %>% select(EC, Accession) %>% distinct()
      )
    )
  )

preX = preX %>%
  # Create matrix column label for each Feature (gap called 0)
  mutate(
    Feature = paste("f", EC, Position, str_replace(Residue, "-", "0"), sep="_")
  ) %>%
  # Summarise Value per Accession and Feature
  ungroup() %>%
  group_by(Accession, Feature) %>%
  summarise(Value = sum(Value)) %>%
  # Spread Feature into columns
  spread(Feature, Value)

# Create X matrix
X = preX %>% ungroup() %>% select(-Accession) %>% as.matrix()

# Set row names
rownames(X) = preX$Accession

# Set NA to 0 occurrences of the Feature in question
X[is.na(X)] = 0

# Remove columns with zero variance
X = X[,apply(X, 2, var) != 0]

# Perform PCA
pcax = prcomp(X %>% scale)
pcpl = as.data.frame(pcax$x) %>% select(PC1, PC2, PC3)
pcpl$Accession = rownames(pcpl)
pcpl = as_tibble(pcpl) %>% inner_join(taxo)

library(scales)
var_pc = percent(pcax$sdev^2 / sum(pcax$sdev^2))[1:3]

gp = ggplot(pcpl, aes(x=PC1, y=PC2, size=PC3, colour=Organism))
gp = gp + geom_point(alpha=0.3)
gp = gp + scale_size(range=c(1,3))
gp = gp + labs(
            x=paste("PC1 (", var_pc[1], ")", sep=""),
            y=paste("PC2 (", var_pc[2], ")", sep=""),
            size=paste("PC3 (", var_pc[3], ")", sep="")
          )
gp = gp + theme_bw()

ggsave("results/carbon_metabolism_pca.normalized.pdf", gp, height=15/2.54, width=24/2.54)

# Write X to file
write_csv(
  as_tibble(X),
  gzfile("intermediate/carbon_metabolism_X.no_cyano.csv.gz"),
  col_names=F
)

# Create y class vector
y = ifelse(preX$Accession %in% cbbx$Genome, 1, 0)

# Write y to file
write_lines(y, "intermediate/carbon_metabolism_y.no_cyano.txt")

# Write feature name vector
write_lines(
  colnames(X), "intermediate/carbon_metabolism_feature_names.no_cyano.txt"
)

# Write data points accession ID vector
write_lines(
  rownames(X), "intermediate/carbon_metabolism_accession_ids.no_cyano.txt"
)
