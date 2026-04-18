# gastropod phylogeny script
# COI + H3
# make fasta, name map, msa, and phylip for IQ-TREE

install.packages(c("rentrez", "seqinr", "ape", "BiocManager", "stringr", "dplyr", "readr"))
BiocManager::install(c("DECIPHER", "Biostrings"))

library(rentrez)
library(seqinr)
library(ape)
library(DECIPHER)
library(Biostrings)
library(stringr)
library(dplyr)
library(readr)

project_dir <- "gastropod_phylogeny"
dir.create(project_dir, showWarnings = FALSE)

markers <- c("COI", "H3")
for (m in markers) {
  dir.create(file.path(project_dir, m), recursive = TRUE, showWarnings = FALSE)
}

coi_accessions <- c(
  "MZ560149.1",
  "KT753980.1",
  "KT753957.1",
  "KT753904.1",
  "KT753914.1",
  "HM431881.1",
  "KT754003.1",
  "KT754011.1",
  "KT753988.1",
  "JQ950211.1",
  "KT753954.1",
  "KT754000.1",
  "KT753975.1",
  "KT753940.1",
  "KT753974.1",
  "KT753994.1",
  "KT753978.1",
  "JQ950201.1",
  "JQ950231.1",
  "KY451220.1",
  "KY451221.1",
  "AB183328.1",
  "AB498778.1",
  "AB498775.1",
  "MW077017.1",
  "MW077023.1",
  "KY451269.1",
  "OZ240484.1",
  "FM999178.1",
  "KT754022.1",
  "HM180684.1",
  "MH581380.1",
  "MH581378.1",
  "KT753938.1",
  "KT753936.1",
  "KC756048.1",
  "KC756032.1",
  "MN322468.1",
  "AY885128.1",
  "KY451284.1",
  "KY451360.1",
  "KY451406.1",
  "MW077025.1",
  "MW077012.1",
  "MW077033.1",
  "MW077031.1",
  "KT753949.1",
  "KP694157.1",
  "HM180636.1"
)

h3_accessions <- c(
  "KT754072.1",
  "KT754110.1",
  "KT754088.1",
  "KT754045.1",
  "KT754034.1",
  "MW057489.1",
  "KT754132.1",
  "KT754140.1",
  "KT754085.1",
  "KT754118.1",
  "KT754129.1",
  "KT754105.1",
  "KT754071.1",
  "KT754104.1",
  "KT754148.1",
  "KT754108.1",
  "KY489294.1",
  "KY489295.1",
  "LC415350.1",
  "MW057529.1",
  "MW057537.1",
  "KT754151.1",
  "HQ834136.1",
  "KT754069.1",
  "KT754067.1",
  "MW057494.1",
  "MN322542.1",
  "HQ834165.1",
  "KY489355.1",
  "KY489371.1",
  "MW057540.1",
  "MW057523.1",
  "MW057548.1",
  "MW057546.1",
  "KT754080.1",
  "MW057509.1"
)

check_duplicate_accessions <- function(accessions, label) {
  dups <- unique(accessions[duplicated(accessions)])
  if (length(dups) > 0) {
    stop(paste0("Duplicate accessions in ", label, ": ", paste(dups, collapse = ", ")))
  }
}

fetch_genbank_fasta <- function(accessions) {
  seqs <- vector("list", length(accessions))
  
  for (i in seq_along(accessions)) {
    acc <- accessions[i]
    message("Downloading FASTA: ", acc)
    
    fasta_txt <- entrez_fetch(
      db = "nucleotide",
      id = acc,
      rettype = "fasta",
      retmode = "text"
    )
    
    seq_obj <- read.fasta(
      textConnection(fasta_txt),
      as.string = FALSE,
      forceDNAtolower = FALSE,
      whole.header = TRUE
    )
    
    seqs[[i]] <- seq_obj[[1]]
  }
  
  seqs
}

fetch_species_names <- function(accessions) {
  out <- character(length(accessions))
  
  for (i in seq_along(accessions)) {
    acc <- accessions[i]
    message("Getting species name: ", acc)
    
    gb_txt <- entrez_fetch(
      db = "nucleotide",
      id = acc,
      rettype = "gb",
      retmode = "text"
    )
    
    org_line <- str_match(gb_txt, "  ORGANISM\\s+([^\\n]+)")[, 2]
    
    if (is.na(org_line) || org_line == "") {
      out[i] <- acc
    } else {
      sp <- gsub("[^A-Za-z0-9 ]", "", org_line)
      sp <- gsub("\\s+", " ", trimws(sp))
      sp <- gsub(" ", "_", sp)
      out[i] <- sp
    }
  }
  
  make.unique(out)
}

sanitize_names <- function(x, max_len = 50) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- substr(x, 1, max_len)
  make.unique(x)
}

write_fasta_list <- function(seqs, file) {
  write.fasta(
    sequences = unname(seqs),
    names = names(seqs),
    file.out = file
  )
}

align_fasta_decipher <- function(infile, outfile) {
  dna <- readDNAStringSet(infile)
  aln <- AlignSeqs(dna, verbose = FALSE)
  writeXStringSet(aln, filepath = outfile, format = "fasta")
  aln
}

save_name_map <- function(accessions, short_names, file) {
  if (length(accessions) != length(short_names)) {
    stop(paste0(
      "Length mismatch: accessions = ", length(accessions),
      ", short_names = ", length(short_names)
    ))
  }
  
  tibble(
    accession = accessions,
    short_name = short_names
  ) %>%
    write_csv(file)
}

process_marker <- function(marker_name, accessions, project_dir = ".") {
  marker_dir <- file.path(project_dir, marker_name)
  dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)
  
  check_duplicate_accessions(accessions, marker_name)
  
  raw_seqs <- fetch_genbank_fasta(accessions)
  final_names <- fetch_species_names(accessions)
  final_names <- sanitize_names(final_names, max_len = 50)
  
  renamed_seqs <- raw_seqs
  names(renamed_seqs) <- final_names
  
  raw_fasta <- file.path(marker_dir, paste0(marker_name, "_raw.fasta"))
  write_fasta_list(renamed_seqs, raw_fasta)
  
  name_map_file <- file.path(marker_dir, paste0(marker_name, "_name_map.csv"))
  save_name_map(accessions, final_names, name_map_file)
  
  cat("\nNames for", marker_name, ":\n")
  print(final_names)
  
  aligned_fasta <- file.path(marker_dir, paste0(marker_name, "_aligned.fasta"))
  aln <- align_fasta_decipher(raw_fasta, aligned_fasta)
  
  aligned_phylip <- file.path(marker_dir, paste0(marker_name, "_aligned.phy"))
  write.dna(as.DNAbin(aln), file = aligned_phylip, format = "sequential")
  
  invisible(list(
    accessions = accessions,
    short_names = final_names,
    raw_fasta = raw_fasta,
    aligned_fasta = aligned_fasta,
    aligned_phylip = aligned_phylip
  ))
}

coi_result <- process_marker(
  marker_name = "COI",
  accessions = coi_accessions,
  project_dir = project_dir
)

if (length(h3_accessions) > 0) {
  h3_result <- process_marker(
    marker_name = "H3",
    accessions = h3_accessions,
    project_dir = project_dir
  )
}
