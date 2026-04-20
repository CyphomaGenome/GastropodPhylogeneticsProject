library(Biostrings)
library(seqinr)
library(pwalign)
library(msa)
library(rentrez)
library(ape)
library(DECIPHER)

Accession_list <- c(
  "KT754072.1",   # Peristernia_reincarnata
  "KT754110.1",   # Peristernia_gemmata
  "KT754088.1",   # Peristernia_nassatula
  "KT754045.1",   # Peristernia_marquesana
  "KT754034.1",   # Peristernia_forskalii
  "MW057489.1",   # Serratifusus_lineatus
  "KT754132.1",   # Triplofusus_giganteus
  "KT754140.1",   # Cinctura_hunteria
  "KT754085.1",   # Fasciolaria_tulipa
  "KT754118.1",   # Fasciolaria_bullisi
  "KT754129.1",   # Aristofusus_excavatus 
  "KT754105.1",   # Fusinus_salisburyi
  "KT754071.1",   # Fusinus_forceps
  "KT754104.1",   # Hemipolygona_armata
  "KT754148.1",   # Leucozonia_nassa
  "KT754108.1",   # Leucozonia_ocellata
  "KY489294.1",   # Buccinastrum_deforme
  "KY489295.1",   # Buccinanops_cochlidium
  "LC415350.1",   # Japeuthria_ferrea
  "MW057529.1",   # Busycon_carica
  "MW057537.1",   # Fulguropsis_pyruloides
  "KT754151.1",   # Mitrella_scripta 
  "HQ834136.1",   # Mitrella_bicincta
  "KT754069.1",   # Dolicholatirus_lancea
  "KT754067.1",   # Dolicholatirus_spiceri 
  "MW057494.1",   # Manaria_kuroharai
  "MN322542.1",   # Melongena_melongena
  "HQ834165.1",   # Nassarius_conoidalis
  "KY489355.1",   # Nassarius_semisulcatus
  "KY489371.1",   # Antillophos_beauii
  "MW057540.1",   # Pisania_pusio
  "MW057523.1",   # Engina_mendicaria
  "MW057548.1",   # Clivipollia_pulchra
  "MW057546.1",   # Caducifer_decapitatus
  "KT754080.1",   # Euthria_cumulata
  "MW057509.1"    # Euthria_walleri
)

fasta_data <- entrez_fetch(
  db = "nucleotide",
  id = Accession_list,
  rettype = "fasta_cds_na"
)

writeLines(fasta_data, "fastas.fasta")

combined_sequences <- readDNAStringSet("fastas.fasta")
combined_sequences <- DNAStringSet(gsub("-", "", as.character(combined_sequences)))

acc <- sub("^lcl\\|([^_]+)_cds_.*$", "\\1", names(combined_sequences))

name_map <- c(
  "KT754072.1" = "Peristernia_reincarnata",
  "KT754110.1" = "Peristernia_gemmata",
  "KT754088.1" = "Peristernia_nassatula",
  "KT754045.1" = "Peristernia_marquesana",
  "KT754034.1" = "Peristernia_forskalii",
  "MW057489.1" = "Serratifusus_lineatus",
  "KT754132.1" = "Triplofusus_giganteus",
  "KT754140.1" = "Cinctura_hunteria",
  "KT754085.1" = "Fasciolaria_tulipa",
  "KT754118.1" = "Fasciolaria_bullisi",
  "KT754129.1" = "Aristofusus_excavatus",
  "KT754105.1" = "Fusinus_salisburyi",
  "KT754071.1" = "Fusinus_forceps",
  "KT754104.1" = "Hemipolygona_armata",
  "KT754148.1" = "Leucozonia_nassa",
  "KT754108.1" = "Leucozonia_ocellata",
  "KY489294.1" = "Buccinastrum_deforme",
  "KY489295.1" = "Buccinanops_cochlidium",
  "LC415350.1" = "Japeuthria_ferrea",
  "MW057529.1" = "Busycon_carica",
  "MW057537.1" = "Fulguropsis_pyruloides",
  "KT754151.1" = "Mitrella_scripta",
  "HQ834136.1" = "Mitrella_bicincta",
  "KT754069.1" = "Dolicholatirus_lancea",
  "KT754067.1" = "Dolicholatirus_spiceri",
  "MW057494.1" = "Manaria_kuroharai",
  "MN322542.1" = "Melongena_melongena",
  "HQ834165.1" = "Nassarius_conoidalis",
  "KY489355.1" = "Nassarius_semisulcatus",
  "KY489371.1" = "Antillophos_beauii",
  "MW057540.1" = "Pisania_pusio",
  "MW057523.1" = "Engina_mendicaria",
  "MW057548.1" = "Clivipollia_pulchra",
  "MW057546.1" = "Caducifer_decapitatus",
  "KT754080.1" = "Euthria_cumulata",
  "MW057509.1" = "Euthria_walleri"
)

names(combined_sequences) <- unname(name_map[acc])

print(data.frame(accession = acc, assigned_name = names(combined_sequences)))

if (any(is.na(names(combined_sequences)))) {
  stop("Some accessions did not get mapped to names.")
}

codon_aligned <- AlignTranslation(combined_sequences)
codon_aligned
my_msa_object <- msa(codon_aligned)
print(my_msa_object, show = "complete")
trimmed_dna <- subseq(codon_aligned, start = 13, end = 339)
trimmed_dna
writeXStringSet(trimmed_dna, filepath = "H3.fasta")
