library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)
library(ape)
library(RADami)
library(DECIPHER)

Accession_list <- c("MG194426", "KT753708", "KT753747", "KT753724", "KT753700", "KT753670", "KT753771", "KT753779", "KT753721", "KT753755", "KT753767", "KT753742", "KT753707", "KT753741", "KT753788", "KT753745", "FJ710096", "FJ710101", "FJ710107", "FM999153", "FM999152", "KT753791", "KT753705", "KT753703", "FM999148", "AB856268", "KT753716", "FJ710099")
Accession_list
fasta_data <- entrez_fetch(db = "nucleotide", 
                            id = Accession_list, 
                            rettype = "fasta")
write(fasta_data, file = "fastas.fasta")

fasta_files <- list.files(pattern = "\\.fasta$")
fasta_files
Strings<-lapply(fasta_files, readDNAStringSet)
Strings
combined_sequences <- do.call(c, Strings)
combined_sequences
width(combined_sequences)

#standardize names
names(combined_sequences)
shortnames <- gsub("^[^ ]* (([^ ]* )[^ ]*).*", "\\1", names(combined_sequences))
finishednames <- gsub(" ", "_", shortnames)
finishednames
names(combined_sequences) <- finishednames
names(combined_sequences)
combined_sequences
MSA_Alignment<-msaMuscle(combined_sequences, verbose = TRUE)
MSA_Alignment
print(MSA_Alignment, show="complete")
alignment_set <- as(MSA_Alignment, "DNAStringSet")
MSA_Subset<-subseq(alignment_set, start = 87, end = 1527)
MSA_Subset
write.DNAStringSet(x = MSA_Subset, format = "phylip", filename = "Peristernia_28S_longer.phy")
