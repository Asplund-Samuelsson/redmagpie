options(width=150)
library(tidyverse)

# Define infiles
dpec_file = "intermediate/deepec.tab.gz"
pfam_file = "intermediate/pfam.tab.gz"
ecal_file = "intermediate/enzyme_alignment.tab.gz"
exgn_file = "intermediate/example_genomes.tab"
acsq_file = "intermediate/enzyme_alignment.accession_sequence.tab.gz"

# Load data
dpec = read_tsv(dpec_file, col_names = c("Accession", "ORF", "EC"))
pfam = read_tsv(pfam_file, col_names = c("Accession", "ORF", "Pfam"))
ecal = read_tsv(ecal_file, col_names = c("EC", "Sequence", "Alignment"))
exgn = read_tsv(exgn_file)
acsq = read_tsv(acsq_file, col_names = c("Accession", "Sequence"))


# 1. Features from enzyme counts

# Prepare examples and features matrix
preX = dpec %>%
  # Tidy up EC
  mutate(EC = str_replace(EC, ":", "")) %>%
  # Exclude PRK and rubisco
  filter(!(EC %in% c("EC2.7.1.19", "EC4.1.1.39"))) %>%
  # Count occurrences of each EC
  group_by(Accession, EC) %>%
  summarise(Count = length(EC)) %>%
  ungroup() %>%
  # Spread ECs into columns
  spread(EC, Count)

# Create X matrix
X = preX %>% select(-Accession) %>% as.matrix()

# Set row names
rownames(X) = preX$Accession

# Set NA to 0 occurrences of the EC in question
X[is.na(X)] = 0

# Create y class vector
y = ifelse(preX$Accession %in% exgn$Genome, 1, 0)

# Write X to file
write_csv(
  as_tibble(X),
  gzfile("intermediate/EC_count_features.X.csv.gz"),
  col_names=F
)

# Write y to file
write_lines(y, "intermediate/EC_count_features.y.txt")

# Write feature name vector
write_lines(colnames(X), "intermediate/EC_count_features.feature_names.txt")

# Write data points accession ID vector
write_lines(rownames(X), "intermediate/EC_count_features.accession_ids.txt")


# 2. Features from Pfams

# Prepare examples and features matrix
preX = pfam %>%
  # Exclude PRK and rubisco
  filter(!(Pfam %in% c("PF02788.16", "PF00016.20", "PF00485.18"))) %>%
  # Count occurrences of each Pfam
  group_by(Accession, Pfam) %>%
  summarise(Count = length(Pfam)) %>%
  ungroup() %>%
  # Spread Pfams into columns
  spread(Pfam, Count)

# Create X matrix
X = preX %>% select(-Accession) %>% as.matrix()

# Set row names
rownames(X) = preX$Accession

# Set NA to 0 occurrences of the EC in question
X[is.na(X)] = 0

# Create y class vector
y = ifelse(preX$Accession %in% exgn$Genome, 1, 0)

# Write X to file
write_csv(
  as_tibble(X),
  gzfile("intermediate/pfam_features.X.csv.gz"),
  col_names=F
)

# Write y to file
write_lines(y, "intermediate/pfam_features.y.txt")

# Write feature name vector
write_lines(colnames(X), "intermediate/pfam_features.feature_names.txt")

# Write data points accession ID vector
write_lines(rownames(X), "intermediate/pfam_features.accession_ids.txt")


# 3. Features from alignments

ecal = ecal %>%
  # Split alignment into amino acids
  mutate(Alignment = str_split(Alignment, "")) %>%
  unnest(Alignment)

ecal = ecal %>%
  # Rename Alignment to Residue
  rename(Residue = Alignment) %>%
  # Add position
  group_by(EC, Sequence) %>%
  mutate(Position = 1:length(Residue)) %>%
  # Add Accession
  inner_join(acsq) %>%
  # Add positive or negative genome Class
  mutate(Class = ifelse(Accession %in% exgn$Genome, 1, 0))

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

resc = ecal %>%
  # Select relevant data
  ungroup() %>%
  select(EC, Residue, Position, Class) %>%
  # Count residues per group
  group_by(EC, Position, Class, Residue) %>%
  summarise(Count = length(Residue))

# Load parallel computing libraries
library(foreach)
library(doMC)
registerDoMC(16)

# Iterate over ECs
pcsp = lapply(unique(resc$EC), function(ec){
    # Subset to current EC
    resc_sub = filter(resc, EC == ec)
    # Iterate over all Positions in the enzyme alignment
    foreach(n=min(resc_sub$Position):max(resc_sub$Position)) %dopar% {
      # Subset to current Position
      resc_sub_sub = filter(resc_sub, Position == n)
      # Create a table of Residue vs Class
      pretbl = resc_sub_sub %>%
        ungroup() %>%
        select(-EC, -Position) %>%
        # Change name of class
        mutate(Class = ifelse(Class == 1, "Positive", "Negative")) %>%
        # Spread Count into columns for each Class
        spread(Class, Count) %>%
        # Replace NA with 0
        mutate(
          Negative = replace_na(Negative, 0),
          Positive = replace_na(Positive, 0)
        )
      # Determine the minimal sum of counts per class
      min_count = min(c(sum(pretbl$Negative), sum(pretbl$Positive)))
      pretbl = pretbl %>%
        # Normalize counts to the minimal sum of counts per class
        mutate(
          Negative = round(min_count * Negative / sum(Negative)),
          Positive = round(min_count * Positive / sum(Positive))
        ) %>%
        # Discard residues with fewer than 5 counts in each position
        filter(Negative >= 5 & Positive >= 5)
      # Transform to "table" format
      tbl = as.matrix(pretbl[,2:3])
      rownames(tbl) = pretbl$Residue
      tbl = as.table(tbl)
      # Perform Chi-square test on the table to evaluate Class separation
      p = chisq.test(tbl)$p.value
      # Return a tibble with EC, Position, and p-value for class separation
      tibble(EC = ec, Position = n, pval = p)
    } %>% bind_rows()
  }) %>% bind_rows() %>%
  # Adjust p values
  mutate(padj = p.adjust(pval, method="BH"))

# Save the p-value table
write_tsv(pcsp, gzfile("intermediate/enzyme_alignment_chi_square_tests.tab.gz"))

# Determine wanted positions (lowest p-value so that all ECs are represented)
wpos = pcsp %>%
  # Remove insignificant positions
  filter(padj < 0.001) %>%
  # Sort by increasing adjusted p-value
  arrange(padj) %>%
  filter(
    !unlist(lapply(
      1:nrow(.), function(x){sum(unique(EC) %in% EC[1:(x-1)]) == 34}
    ))
  )

preX = ecal %>%
  # Keep only wanted positions
  inner_join(select(wpos, EC, Position))

# Obtain "consensus" sequence
cons = resi %>%
  # Filter to significant positions
  inner_join(select(wpos, EC, Position)) %>%
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
        crossing(EC=unique(preX$EC), Accession=c(exgn$Genome, exgn$Relative)),
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

# Create y class vector
y = ifelse(preX$Accession %in% exgn$Genome, 1, 0)

# Write X to file
write_csv(
  as_tibble(X),
  gzfile("intermediate/alignment_features.X.csv.gz"),
  col_names=F
)

# Write y to file
write_lines(y, "intermediate/alignment_features.y.txt")

# Write feature name vector
write_lines(colnames(X), "intermediate/alignment_features.feature_names.txt")

# Write data points accession ID vector
write_lines(rownames(X), "intermediate/alignment_features.accession_ids.txt")
