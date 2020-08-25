options(width=150)
library(tidyverse)

# Define infiles
desc_file = "data/Pfam-A_32.0.NAME_ACC_DESC.tab"
pfam_file = "intermediate/pfam.tab.gz"
exgn_file = "intermediate/example_genomes.tab"

# Load data
desc = read_tsv(desc_file, col_names=c("Name","Pfam","Description"))
pfam = read_tsv(pfam_file, col_names=c("Genome", "ORF", "Pfam"))
exgn = read_tsv(exgn_file)

# Find photosynthetic Pfams
phot = filter(desc, grepl("photosy", Description, ignore.case=T))

# Count photosynthetic and other ORFs per genome
func = pfam %>%
  left_join(select(phot, Pfam, Name)) %>%
  mutate(Function = ifelse(is.na(Name), "Other", "Photosynthesis")) %>%
  group_by(Genome, Function) %>% summarise(Count = length(ORF))

func = func %>%
  spread(Function, Count) %>%
  mutate(
    Photosynthesis = replace_na(Photosynthesis, 0),
    Other = replace_na(Other, 0)
  )

# Determine if genomes are CBB-positive or CBB-negative
func = func %>%
  mutate(CBB = ifelse(Genome %in% exgn$Genome, "CBB-positive", "CBB-negative"))

# Check on Bradyrhizobium sp. ORS 278, which is a known photosynthesizer
# filter(func, grepl("GCF_000026145", Genome))
# It has three of the phototsynthesis Pfams

# Count the number of photosynthetic genomes
# filter(func, Photosynthesis >= 3) %>% pull(CBB) %>% table
# CBB-negative CBB-positive
#          115          113

# Save the tables
phot = phot %>%
  rename(Feature = Pfam) %>%
  select(Feature, Name, Description)
write_tsv(phot, "results/Supplementary_Photosynthesis.tab")
write_tsv(func, "results/genome_photosynthesis_status.tab")
