library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)
library(ape)
library(RADami)
library(DECIPHER)

stringH3<-readDNAStringSet("H3.fasta")
stringCOI<-readDNAStringSet("COI.fasta")
write.DNAStringSet(stringH3, format="phylip", filename="Peristernia_H3.phy")
write.DNAStringSet(stringCOI, format="phylip", filename="Peristernia_COI.phy")

###these lines of code turn uncertainty into polytomies
my_tree <- read.tree("partitions.txt.treefile")
my_tree
#Define threshold
threshold <- 60
threshold
#Convert node labels to numeric (sometimes they load as characters)
bootstraps <- as.numeric(my_tree$node.label)
bootstraps
#Find the edges (branches) leading to nodes with low support
#Note: Node indices in ape start after the number of tips
nodes_to_collapse <- which(bootstraps < threshold) + length(my_tree$tip.label)
nodes_to_collapse
#Set the edge lengths of those weak nodes to 0
my_tree$edge.length[my_tree$edge[,2] %in% nodes_to_collapse] <- 0
#Collapse all zero-length branches into polytomies
final_tree <- di2multi(my_tree, tol = 1e-08)
final_tree
write.tree(my_tree, "ConcatenatedCaenogastropoda.tre")

