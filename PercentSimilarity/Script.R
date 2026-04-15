#Part A
library(Biostrings)
fasta_files <- list.files(pattern = "\\.fasta$")
fasta_files
sequences <- readDNAStringSet(fasta_files)
sequences #these 3 lines of code are optional but help visualize the data
class(sequences)
names(sequences)
#Part B
library(pwalign)
aln <- pwalign::pairwiseAlignment(sequences[1], sequences[2])
aln
pwalign::pid(aln)
#90.68592 for Penion chathamensis and Penion sulcatus
aln2 <- pwalign::pairwiseAlignment(sequences[2], sequences[3])
aln2
pwalign::pid(aln2)
#99.43521 for Penion sulcatus and Peristernia sulcata
