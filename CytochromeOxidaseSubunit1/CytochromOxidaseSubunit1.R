library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)
library(ape)
library(DECIPHER)

Accession_list <- c(
  "MZ560149.1","KT753980.1","KT753957.1","KT753904.1","KT753914.1",
  "HM431881.1","KT754003.1","KT754011.1","KT753988.1","JQ950211.1",
  "KT753954.1","KT754000.1","KT753975.1","KT753940.1","KT753974.1",
  "KT753994.1","KT753978.1","JQ950201.1","JQ950231.1","KY451220.1",
  "KY451221.1","AB183328.1","AB498778.1","AB498775.1","MW077017.1",
  "MW077023.1","KY451269.1","OZ240484.1","FM999178.1","KT754022.1",
  "HM180684.1","MH581380.1","MH581378.1","KT753938.1","KT753936.1",
  "KC756048.1","KC756032.1","MN322468.1","AY885128.1","KY451284.1",
  "KY451360.1","KY451406.1","MW077025.1","MW077012.1","MW077033.1",
  "MW077031.1","KT753949.1","KP694157.1","HM180636.1"
)

fasta_data <- entrez_fetch(
  db = "nucleotide",
  id = Accession_list,
  rettype = "fasta_cds_na"
)

writeLines(fasta_data, "COI_raw.fasta")

# read ONLY this file
combined_sequences <- readDNAStringSet("COI_raw.fasta")

# clean names safely
shortnames <- gsub("^[^ ]* (([^ ]* )[^ ]*).*", "\\1", names(combined_sequences))
finishednames <- gsub(" ", "_", trimws(shortnames))
names(combined_sequences) <- make.unique(finishednames)

print(names(combined_sequences))

# codon alignment
codon_aligned <- AlignTranslation(combined_sequences)

# optional check
msa(codon_aligned)

# trim region
trimmed_dna <- subseq(codon_aligned, start = 43, end = 501)

# export FASTA
writeXStringSet(trimmed_dna, "COI.fasta")