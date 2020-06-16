options(width=150)
library(tidyverse)

# Define infiles
posg_file = "data/positive_genomes.txt"
rubs_file = "data/rubisco.txt"
prks_file = "data/prk.txt"
dpec_file = "intermediate/deepec.tab.gz"
pfam_file = "intermediate/pfam.tab.gz"
psmd_file = "data/assembly_table.tab"

# Load data
posg = scan(posg_file, character())
rubs = scan(rubs_file, character())
prks = scan(prks_file, character())
dpec = read_tsv(dpec_file, col_names = c("Accession", "ORF", "Feature"))
pfam = read_tsv(pfam_file, col_names = c("Accession", "ORF", "Feature"))
psmd = read_tsv(psmd_file)

# Function for getting element of split string
get_elem = function(string, pattern, position, n=Inf){
  unlist(lapply(str_split(string, pattern, n), "[[", position))
}

# Format data for analysis
rupr = tibble(
  Accession = get_elem(c(rubs, prks), "\\|", 1),
  ORF = get_elem(c(rubs, prks), "\\|", 2, 2),
  Feature = c(rep("Rubisco", length(rubs)), rep("Prk", length(prks)))
)

feat = bind_rows(rupr, dpec, pfam) %>% filter(Accession %in% posg)

feat = feat %>%
  mutate(
    Contig = get_elem(ORF,"\\|",2) %>% get_elem(":",1) %>% get_elem("_",2,2),
    xStart = as.numeric(get_elem(ORF, ":", 2)),
    xEnd = as.numeric(get_elem(ORF, ":", 3)),
    Strand = ifelse(xStart > xEnd, "-", "+"),
    Start = ifelse(Strand == "+", xStart, xEnd),
    End = ifelse(Strand == "+", xEnd, xStart)
  ) %>%
  select(-xStart, -xEnd)

# Exclude Prk and Rubisco identified by DeepEC or Pfam
feat = feat %>%
  filter(
    !(Feature %in% c(
      "EC:2.7.1.19", "EC:4.1.1.39",
      "PF02788.16", "PF00016.20", "PF00485.18", "PF00101.20"
    ))
  )

# Calculate gene midpoint
feat = feat %>%
  mutate(Midpoint = (Start + End)/2)

# Calculate position on Contig
feat = feat %>%
  select(Accession, Contig, ORF, Midpoint) %>%
  distinct() %>%
  group_by(Accession, Contig) %>%
  mutate(Position = rank(Midpoint)) %>%
  inner_join(feat) %>%
  ungroup()

# Calculate distance to Rubisco and Prk of every feature
cbbf = feat %>%
  filter(Feature %in% c("Rubisco", "Prk")) %>%
  rename(
    cStart = Start, cEnd = End,
    cFeature = Feature, cStrand = Strand,
    cPosition = Position, cMidpoint = Midpoint
  ) %>%
  select(-ORF)

cbbd = inner_join(feat, cbbf) %>% filter(Feature != cFeature)

# Only consider features that are at most 50 genes away from Rubisco or Prk
# cbbd = cbbd %>% filter(abs(Position - cPosition) <= 50)

# Determine position relative to the CBB feature and calculate distance
cbbd = cbbd %>% mutate(Distance = abs(Position - cPosition))

# Add information about plasmid status
cbbd = cbbd %>%
  left_join(
    psmd %>%
      rename(
        Molecule = `Assigned-Molecule-Location/Type`,
        Contig = `GenBank-Accn`
      ) %>%
      select(Molecule, Contig)
  ) %>%
  mutate(
    Molecule = ifelse(
      Molecule %in% c("Chromosome", "Plasmid"),
      Molecule, "Unknown"
    )
  )

# Calculate median distance and count occurrences
cbop = cbbd %>%
  mutate(Strand = ifelse(Strand == cStrand, "Same", "Opposite")) %>%
  group_by(Feature, cFeature, Strand) %>%
  summarise(
    minD = min(Distance), medD = median(Distance),
    maxD = max(Distance), meanD = round(mean(Distance)),
    Count = length(Feature),
    locChr = sum(Molecule == "Chromosome"),
    locPsm = sum(Molecule == "Plasmid"),
    locUnk = sum(Molecule == "Unknown"),
    fracPsm = locPsm / (locPsm + locChr)
  ) %>%
  arrange(medD, -Count)

# Load feature annotations from supplementary and add to table
supp = read_tsv("results/Supplementary_Consensus_rank.tab")
supp = supp %>%
  select("Rank", "Feature_Type", "Feature", "Name", "Description", "KEGG_EC")

cbop = cbop %>%
  mutate(Feature = str_remove(Feature, ":")) %>%
  left_join(supp) %>%
  select(Rank, Feature_Type, Feature, Name, everything())

# Save table
write_tsv(cbop, "results/Supplementary_Gene_proximity.tab")
