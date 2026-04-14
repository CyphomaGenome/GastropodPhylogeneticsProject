library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)
library(ape)
library(RADami)
library(DECIPHER)

accession_list1 <- c(
  "CM009772", "OQ682616", "KM100140", "KP716638", "MH700493", 
  "KX688549", "MZ568919", "KY196442", "KY679829", "KY679834", 
  "MZ435264", "MN557850", "MK256614", "MN082637", "MH682098", 
  "NC_037771", "KU878411", "MZ321058", "OQ695489", "KC757644", 
  "OR209184", "OM063153", "MK251987", "OQ511529", "EU827193", 
  "OP021744", "OP738004", "LC469295", "OP714184", "OP723877", 
  "OP714183", "MK783263", "OR522697", "MW376482", "MK507894", 
  "MK507893", "MK500871", "MK478016", "MK478017", "MK478019", 
  "KP716634", "KY697388", "KY697386", "OR426588", 
  "JF284698", "MZ196218", "OP503947", "EU440735", "MT410857", 
  "MW244819", "MW244822", "MW244821", "MW244818", "MW244817", 
  "MZ359282", "MW244823", "MK388726", "HM174255", "HM174254"
)
length(accession_list1)
fasta_data1 <- entrez_fetch(db = "nucleotide", 
                           id = accession_list1, 
                           rettype = "fasta")
write(fasta_data1, file = "fastas1")
accession_list2 <- c(
  "HM174252", "OR251120", "KP716635", "EU827200", "MK035573", 
  "MH198162", "MT075610", "MN462591", "MW548267", "MZ618619", 
  "KP716637", "MK335921", "MK583341", "KT826694", "MT762151", 
  "MG744570", "MT755649", "MK583343", "MH198160", "MH198164", 
  "MH931230", "MH198159", "MH198163", "MK558054", "MK558053", 
  "MN462590", "EU827199", "MN871953", "MN125492", "MT232845", 
  "MZ934412", "EU827195", "MT415926", "MN462589", "MW550295", 
  "MW550289", "MW550292", "MW550287", "MK395390", "MG786490", 
  "EU827194", "GU196685", "MN400053", "MW044625", "HQ416443", 
  "KT211493", "MK577482", "MT240815", "MT240810", "KX263257", 
  "KX263255", "MH308389", "KX263259", "MH308401", "MH308401", 
  "MH308400", "MN583349", "MH308407", "MH308406", "MH308404", 
  "MH308403", "MH308392", "MH308393", "EU827196", "MK251986", 
  "MH308395", "MH308399", "MH308398", "MH308391", "MH308390", 
  "KX263260", "EU827197", "MH308394", "OM048759", "OM048764"
)
length(accession_list2)
fasta_data2 <- entrez_fetch(db = "nucleotide", 
                            id = accession_list2, 
                            rettype = "fasta")
write(fasta_data2, file = "fastas2")

Strings<-lapply(c("fastas1", "fastas2"), readDNAStringSet)
Strings
combined_sequences <- do.call(c, Strings)
combined_sequences
#standardize names
names(combined_sequences)
shortnames <- gsub("^[^ ]* (([^ ]* )[^ ]*).*", "\\1", names(combined_sequences))
finishednames <- gsub(" ", "_", shortnames)
finishednames
names(combined_sequences) <- finishednames
names(combined_sequences)
combined_sequences
writeXStringSet(combined_sequences, "sequences.fasta", format="fasta")

#in bash Bioinformatics folder, install mafft
#wget https://mafft.cbrc.jp/alignment/software/mafft_7.526-1_amd64.deb
#sudo dpkg -i mafft_7.526-1_amd64.deb
#Navigate to proper folder in bash
#cd /mnt/c/Bioinformatics/Bioinformatics/GastropodPhylogeneticsProject/Li_et_al_full_replicate
#mafft --auto --leavegappyregion  sequences.fasta > aligned_sequences.fasta
#probably leave off --adjustdirection

#install trimal
#conda install -c bioconda trimal -y
#run trimal
#trimal -in aligned_sequences.fasta -out trimmed_alignment.fasta -automated1

#Run IQ-Tree with fast bootstraps
#iqtree2 -s trimmed_alignment.fasta -bb 10000 -nt AUTO

###the next lines of code turn uncertainty into polytomies
my_tree <- read.tree("aligned_sequences.fasta.treefile")
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
pruned_tree <- drop.tip(final_tree, "_R_Godlewskia_godlewskia")
pruned_tree
write.tree(pruned_tree, "Li_et_al_fullreplicate60threshold.tre")


###################################################################################
STORAGE

#print(MSA_Alignment, show="complete")
alignment_set <- as(MSA_Alignment, "DNAStringSet")
MSA_Subset<-subseq(alignment_set, start = 79, end = 1496)
MSA_Subset
MSA_stringset <- as(MSA_Subset, "DNAStringSet")
MSA_stringset



MSA_Subset |> 
  as.matrix() |> 
  as.DNAbin() |> 
  write.dna(file = "16S_rRNA.phy", 
            format = "sequential", 
            nbcol = -1, 
            colsep = "", 
            indent = 35)

#install.packages("ips")
library(ips)
alignment_bin <- as.DNAbin(as.matrix(MSA_Subset))
#On next step adjust multiplier for level of trimming (1=untrimmed, 0.5 means 90%=trim where 90% of sequences have gaps)
trimmed_alignment <- deleteGaps(alignment_bin, gap.max = nrow(alignment_bin) * 1)
trimmed_alignment
write.dna(trimmed_alignment, 
          file = "gastropod_16S_NoPomacea.phy", 
          format = "sequential", 
          nbcol = -1, 
          colsep = "")


library(DECIPHER)
codon_aligned <- AlignTranslation(combined_sequences)
codon_aligned
trimmed_dna <- subseq(codon_aligned, start = 1, end = 1545)
trimmed_dna
#BiocManager::install(c("IRanges", "XVector"), force=TRUE)
#remotes::install_github("andrew-hipp/RADami", upgrade = "never")
library(RADami)
write.DNAStringSet(x = trimmed_dna, format = "phylip", filename = "CytrochomeMSA.phy")


#names(combined_sequences) <- c()
#combined_sequences
#writeXStringSet(combined_sequences, "botulinumcompletesequences.fasta")


#these are valid lines of code but not neccesary
phylip_format <- msaConvert(MSA_Alignment, type="phangorn::phyDat")
phylip_format
#next single line of code writes phylip file
write.phylip(MSA_Alignment, file="botulinum.phy")

#Then go to Bash
#cd /mnt/c/Bioinformatics/botulinum
#iqtree2 -s botulinum.phy
#botulinum.phy.iqtree has the information


#Commands used for iqtree installation
#wget https://github.com/iqtree/iqtree2/releases/download/v2.2.2.7/iqtree-2.2.2.7-Linux.tar.gz
#/mnt/c/Bioinformatics/iqtree-2.2.2.7-Linux/bin
#export PATH="$PATH:/mnt/c/Bioinformatics/iqtree-2.2.2.7-Linux/bin"

#NOTE: command for DNA with codons
#iqtree -s coding_gene.phy -st CODON 
#NOTE: command for DNA with codons using invertebrate mitochondrial DNA
#iqtree -s coding_gene.phy -st CODON5

#BiocManager::install("DECIPHER")
library(DECIPHER)
codon_aligned <- AlignTranslation(combined_sequences)
codon_aligned
trimmed_dna <- subseq(codon_aligned, start = 1, end = 4119)
trimmed_dna
#BiocManager::install(c("IRanges", "XVector"), force=TRUE)
#remotes::install_github("andrew-hipp/RADami", upgrade = "never")
library(RADami)
write.DNAStringSet(x = trimmed_dna, format = "phylip", filename = "botulinumcodons.phy")

#Then go to Bash
#cd /mnt/c/Bioinformatics/botulinum
#iqtree2 -s botulinumcodons.phy -st CODON
#botulinumcodons.phy.iqtree has the information

#next lines of code try to use better codon model
write.DNAStringSet(x = trimmed_dna, format = "phylip", filename = "botulinumcodonmodel.phy")
#iqtree2 -s botulinumcodonmodel.phy -m MFP+CODON -st CODON

Strings2<-lapply(c("botulinum.fasta","boNTG.txt","boNTE.txt","boNTF.txt","boNTB.txt", "boNTEn.fasta"), readDNAStringSet)
Strings2
combined_sequences2 <- do.call(c, Strings2)
combined_sequences2
names(combined_sequences2) <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "C1", "A1", "A2", "B1", "B2", "F1", "B3", "B4", "G1", "G2", "E8", "F2", "B5", "En1")
combined_sequences2
writeXStringSet(combined_sequences2, "botulinumcompletesequenceswithancestral.fasta")
MSA_Alignment2<-msaMuscle(combined_sequences2)
MSA_Alignment2
print(MSA_Alignment2, show="complete")
codon_aligned2 <- AlignTranslation(combined_sequences2)
codon_aligned2
trimmed_dna2 <- subseq(codon_aligned2, start = 1, end = 4230)
trimmed_dna2
write.DNAStringSet(x = trimmed_dna2, format = "phylip", filename = "botulinumwithancestralcodons.phy")

#Then go to Bash
#cd /mnt/c/Bioinformatics/botulinum
#iqtree2 -s botulinumwithancestralcodons.phy -st CODON
#botulinumwithancestralcodons.phy.iqtree has the information


#additional commands
#/mnt/c/Bioinformatics/iqtree-1.6.12-Linux/iqtree-1.6.12-Linux/bin
#export PATH="$PATH:/mnt/c/Bioinformatics/iqtree-1.6.12-Linux/iqtree-1.6.12-Linux/bin"
#iqtree -s botulinumcodons.phy -st CODON
#still does not work

#Once tree is finished
# 1. Load your tree
my_tree <- read.tree("gastropod_16S_NoPomacea.phy.raxml.support")

# 2. Create a translation vector
# The names must match exactly what is currently in the tree
old_names <- my_tree$tip.label
old_names
new_names <- c("Neptunea_purpurea", "Neptunea_lamellosa", "Buccinum_aniwanum", "Buccinum_morchianum", "Japeuthria_ferrea", 
               "Siphonalia_cassidariaeformis", "Kelletia_kelletii", "Kelletia_lischkei", "Penion_chathamensis", "Peristernia_sulcata", 
               "Neobuccinum_eatoni", "Pollia_tincta", "Engina_mendicaria", "Pisania_striata", "Pisania_pusio", 
               "Cantharus_multangulus", "Granulifusus_niponicus", "Fusinus_akitai", "Colubraria_nitidula", "Colubraria_obscura", 
               "Metula_amosi", "Nassarius_semisulcatus", "Zeuxis_siquijorensis", "Nassarius_festivus", "Microvoluta_sp", 
               "Vexillum_plicarium", "Latiromitra_sp", "Ficus_subintermedia", "Oliva_mustelina", "Oliva_spicata", 
               "Olivella_volutella", "Mitra_lens", "Polystira_picta", "Comitas_kaderlyi", "Aforia_magnifica", 
               "Conus_praecellens", "Tonna_luteostoma", "Biplex_perca", "Fusitriton_oregonensis", "Plesiotriton_sp", 
               "Littorina_brevicula", "Cancellaria_cooperi", "Cancellaria_cancellata", "Cancellaria_sinensis", "Babylonia_japonica", 
               "Babylonia_lutosa", "Mitrella_bicincta", "Thalessa_savignyi", "Thais_clavigera", "Thais_haemastoma", 
               "Muricanthus_radix", "Siratus_beauii", "Cronia_sp", "Drupella_cornus", "Nucella_lapillus", "Hemifusus_tuba", "Volema_myristica", "Melongena_patula", "Cominella_adspersa", "Antillophos_laevis", 
               "Paraeuthria_plumbea", "Burnupena_cincta")
new_names

# 3. Swap the names
my_tree$tip.label <- new_names

# 4. Save the renamed tree
write.tree(my_tree, "gastropod_final_named.tre")

#iqtree2 -s 16S_rRNA.phy
#or with fast bootstraps
#iqtree2 -s 16S_rRNA.phy -B 1000 -alrt 1000 -T AUTO
