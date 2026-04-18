# ================================
# Gastropod phylogeny starter script
# For: COI and Histone H3
# ================================

# Install if needed:
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

# ----------------
# USER SETTINGS
# ----------------

project_dir <- "gastropod_phylogeny"
dir.create(project_dir, showWarnings = FALSE)

# Set paths for each marker
markers <- c("COI", "H3")
for (m in markers) {
  dir.create(file.path(project_dir, m), recursive = TRUE, showWarnings = FALSE)
}

# Example accession lists
# Replace these with your real lists
coi_accessions <- c(
  "MZ560149.1",   # Peristernia_reincarnata
  "KT753980.1",   # Peristernia_gemmata
  "KT753957.1",   # Peristernia_nassatula
  "KT753904.1",   # Peristernia_forskalii
  "KT753914.1",   # Peristernia_marquesana
  "HM431881.1",   # Penion_chathamensis
  "KT754003.1",   # Triplofusus_giganteus
  "KT754011.1",   # Cinctura_hunteria
  "KT753988.1",   # Fasciolaria_bullisi
  "JQ950211.1",   # Serratifusus_lineatus
  "KT753954.1",   # Fasciolaria_tulipa
  "KT754000.1",   # Aristofusus_excavatus 
  "KT753975.1",   # Fusinus_salisburyi
  "KT753940.1",   # Fusinus_forceps
  "KT753974.1",   # Hemipolygona_armata
  "KT753994.1",   # Leucozonia_nassa
  "KT753978.1",   # Leucozonia_ocellata
  "JQ950201.1",   # Belomitra_caudata
  "JQ950231.1",   # Belomitra_paschalis
  "KY451220.1",   # Buccinastrum_deforme
  "KY451221.1",   # Buccinanops_cochlidium
  "AB183328.1",   # Buccinum_striatissimum
  "AB498778.1",   # Neptunea_arthritica
  "AB498775.1",   # Neptunea_polycostata 
  "MW077017.1",   # Busycon_carica
  "MW077023.1",   # Fulguropsis_pyruloides
  "KY451269.1",   # Chauvetia_candidissima
  "OZ240484.1",   # Colus_islandicus
  "FM999178.1",   # Colubraria_reticulata
  "KT754022.1",   # Mitrella_scripta 
  "HM180684.1",   # Mitrella_bicincta
  "MH581380.1",   # Tonna_sulcosa
  "MH581378.1",   # Tonna_galea
  "KT753938.1",   # Dolicholatirus_lancea
  "KT753936.1",   # Dolicholatirus_spiceri
  "KC756048.1",   # Eosipho_smithi
  "KC756032.1",   # Manaria_kuroharai
  "MN322468.1",   # Melongena_melongena
  "AY885128.1",   # Melongena_patula
  "KY451284.1",   # Nassarius_conoidalis
  "KY451360.1",   # Nassarius_semisulcatus
  "KY451406.1",   # Antillophos_beauii
  "MW077025.1",   # Pisania_pusio
  "MW077012.1",   # Engina_mendicaria
  "MW077033.1",   # Clivipollia_pulchra
  "MW077031.1",   # Caducifer_decapitatus
  "KT753949.1",   # Euthria_cumulata
  "KP694157.1",   # Buccinulum_linea
  "HM180636.1"    # Kelletia_lischkei
)

h3_accessions <- c(
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
  "MW057509.1"   # Euthria_walleri
)

# Optional manual short-name map
# If an accession is missing here, the script will generate one from GenBank metadata
manual_name_map <- c(
  "MZ560149.1" = "Peristernia_reincarnata",
  "KT753980.1" = "Peristernia_gemmata",
  "MZ560437.1" = "Peristernia_nassatula",
  "KT753904.1" = "Peristernia_forskalii",
  "KT753914.1" = "Peristernia_marquesana",
  "HM431881.1" = "Penion_chathamensis",
  "KT754003.1" = "Triplofusus_giganteus",
  "KT754011.1" = "Cinctura_hunteria",
  "KT753988.1" = "Fasciolaria_bullisi",
  "JQ950211.1" = "Serratifusus_lineatus",
  "AF373884.1" = "Fasciolaria_tulipa",
  "KT754000.1" = "Aristofusus_excavatus",
  "KT753975.1" = "Fusinus_salisburyi",
  "KT753940.1" = "Fusinus_forceps",
  "KT753974.1" = "Hemipolygona_armata",
  "KT753994.1" = "Leucozonia_nassa",
  "KT753978.1" = "Leucozonia_ocellata",
  "JQ950201.1" = "Belomitra_caudata",
  "KY451220.1" = "Buccinastrum_deforme",
  "KY451221.1" = "Buccinanops_cochlidium",
  "PX626180.1" = "Buccinum_tenuissimum",
  "AB183328.1" = "Buccinum_striatissimum",
  "AB432884.1" = "Neptunea_arthritica",
  "FJ710092.1" = "Neptunea_polycostata",
  "ON629600.1" = "Japeuthria_ferrea",
  "U86323.1" = "Busycon_carica",
  "MW077023.1" = "Fulguropsis_pyruloides",
  "KY451269.1" = "Chauvetia_candidissima",
  "OZ240484.1" = "Colus_islandicus",
  "FM999178.1" = "Colubraria_reticulata",
  "KT754022.1" = "Mitrella_scripta",
  "HQ834055.1" = "Mitrella_bicincta",
  "OQ828243.1" = "Tonna_sulcosa",
  "MW148581.1" = "Tonna_galea",
  "MW181899.1" = "Dolicholatirus_lancea",
  "JQ950198.1" = "Eosipho_smithi",
  "MN752196.1" = "Manaria_kuroharai",
  "MN322468.1" = "Melongena_melongena",
  "AY885128.1" = "Melongena_patula",
  "MN389172.1" = "Nassarius_conoidalis",
  "KY451360.1" = "Nassarius_semisulcatus",
  "KY451406.1" = "Antillophos_beauii",
  "MW077025.1" = "Pisania_pusio",
  "MN389158.1" = "Engina_mendicaria",
  "MZ559419.1" = "Clivipollia_pulchra",
  "MW278449.1" = "Caducifer_decapitatus",
  "KT753949.1" = "Euthria_cumulata",
  "KP694157.1" = "Buccinulum_linea",
  "HM180636.1" = "Kelletia_lischkei",
  
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
  "MW057494.1" = "Manaria_kuroharai",
  "MN322542.1" = "Melongena_melongena",
  "HQ834165.1" = "Nassarius_conoidalis",
  "KY489355.1" = "Nassarius_semisulcatus",
  "KY489371.1" = "Antillophos_beauii",
  "MW057540.1" = "Pisania_pusio",
  "FM999175.1" = "Pisania_striata",
  "MW057523.1" = "Engina_mendicaria",
  "MW057548.1" = "Clivipollia_pulchra",
  "MW057546.1" = "Caducifer_decapitatus",
  "KT754080.1" = "Euthria_cumulata",
  "MW057509.1" = "Euthria_walleri"
)

# Trim settings
# gap_threshold = proportion of sequences allowed to have gaps at a site
# 0.40 means remove columns where >60% are gaps
gap_threshold <- 0.40

# Bootstrap display threshold from your notes
support_threshold <- 60

# IQ-TREE executable name/path
iqtree_bin <- "iqtree2"   # change if needed, e.g. "/usr/local/bin/iqtree2"

# ----------------
# HELPER FUNCTIONS
# ----------------

fetch_genbank_fasta <- function(accessions) {
  seqs <- vector("list", length(accessions))
  names(seqs) <- accessions
  
  for (acc in accessions) {
    message("Downloading: ", acc)
    fasta_txt <- entrez_fetch(
      db = "nucleotide",
      id = acc,
      rettype = "fasta",
      retmode = "text"
    )
    seq_obj <- read.fasta(textConnection(fasta_txt), as.string = FALSE, forceDNAtolower = FALSE)
    seqs[[acc]] <- seq_obj[[1]]
  }
  
  seqs
}

fetch_genbank_metadata <- function(accessions) {
  out <- vector("list", length(accessions))
  names(out) <- accessions
  
  for (acc in accessions) {
    message("Getting metadata: ", acc)
    gb_txt <- entrez_fetch(
      db = "nucleotide",
      id = acc,
      rettype = "gb",
      retmode = "text"
    )
    
    org_line <- str_match(gb_txt, "  ORGANISM\\s+([^\\n]+)")[,2]
    if (is.na(org_line)) org_line <- acc
    
    org_line <- str_trim(org_line)
    org_line <- gsub("[^A-Za-z0-9 ]", "", org_line)
    parts <- unlist(str_split(org_line, "\\s+"))
    
    short_name <- if (length(parts) >= 2) {
      paste(parts[1], parts[2], sep = "_")
    } else {
      parts[1]
    }
    
    short_name <- make.unique(short_name)
    out[[acc]] <- short_name
  }
  
  unlist(out)
}

sanitize_names <- function(x, max_len = 18) {
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

trim_gap_columns <- function(aln, gap_threshold = 0.40) {
  # Convert aligned DNAStringSet to character matrix
  aln_char <- as.matrix(aln)
  
  # proportion gaps per column
  gap_prop <- apply(aln_char, 2, function(col) mean(col %in% c("-", ".", "N", "n")))
  
  keep <- gap_prop <= (1 - gap_threshold)
  trimmed <- aln[, keep]
  trimmed
}

save_name_map <- function(accessions, short_names, file) {
  tibble(
    accession = accessions,
    short_name = short_names
  ) %>%
    write_csv(file)
}

run_iqtree <- function(alignment_file, outdir, prefix, iqtree_bin = "iqtree2") {
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(outdir)
  
  cmd <- sprintf(
    '%s -s "%s" -m MFP -B 1000 -alrt 1000 -nt AUTO --prefix "%s"',
    iqtree_bin, basename(alignment_file), prefix
  )
  
  message("Running IQ-TREE:\n", cmd)
  system(cmd)
}

collapse_low_support <- function(treefile, threshold = 60, outfile = NULL) {
  tr <- read.tree(treefile)
  
  if (!is.null(tr$node.label)) {
    suppressWarnings({
      support_vals <- as.numeric(tr$node.label)
    })
    
    low <- which(!is.na(support_vals) & support_vals < threshold)
    if (length(low) > 0) {
      tr$node.label[low] <- NA
    }
  }
  
  if (!is.null(outfile)) {
    write.tree(tr, file = outfile)
  }
  
  tr
}

process_marker <- function(marker_name, accessions, manual_name_map = NULL,
                           project_dir = ".", gap_threshold = 0.40,
                           iqtree_bin = "iqtree2", run_tree = FALSE) {
  
  marker_dir <- file.path(project_dir, marker_name)
  dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Download sequences
  raw_seqs <- fetch_genbank_fasta(accessions)
  
  # 2. Build names
  auto_names <- fetch_genbank_metadata(accessions)
  
  final_names <- auto_names
  if (!is.null(manual_name_map)) {
    overlap <- intersect(names(manual_name_map), accessions)
    final_names[overlap] <- manual_name_map[overlap]
  }
  final_names <- sanitize_names(final_names, max_len = 18)
  
  # 3. Rename sequences
  renamed_seqs <- raw_seqs
  names(renamed_seqs) <- final_names
  
  # 4. Save raw fasta + name map
  raw_fasta <- file.path(marker_dir, paste0(marker_name, "_raw.fasta"))
  write_fasta_list(renamed_seqs, raw_fasta)
  
  name_map_file <- file.path(marker_dir, paste0(marker_name, "_name_map.csv"))
  save_name_map(accessions, final_names, name_map_file)
  
  cat("\nShort names for", marker_name, ":\n")
  print(final_names)
  
  # 5. Align
  aligned_fasta <- file.path(marker_dir, paste0(marker_name, "_aligned.fasta"))
  aln <- align_fasta_decipher(raw_fasta, aligned_fasta)
  
  # 6. Trim
  trimmed <- trim_gap_columns(aln, gap_threshold = gap_threshold)
  trimmed_fasta <- file.path(marker_dir, paste0(marker_name, "_trimmed.fasta"))
  writeXStringSet(trimmed, filepath = trimmed_fasta, format = "fasta")
  
  # 7. Also save phylip for IQ-TREE if needed
  trimmed_phylip <- file.path(marker_dir, paste0(marker_name, "_trimmed.phy"))
  write.dna(as.DNAbin(trimmed), file = trimmed_phylip, format = "sequential")
  
  # 8. Optional tree run
  if (run_tree) {
    run_iqtree(trimmed_phylip, marker_dir, prefix = marker_name, iqtree_bin = iqtree_bin)
    
    treefile <- file.path(marker_dir, paste0(marker_name, ".treefile"))
    if (file.exists(treefile)) {
      collapse_low_support(
        treefile,
        threshold = support_threshold,
        outfile = file.path(marker_dir, paste0(marker_name, "_collapsed_", support_threshold, ".tre"))
      )
    }
  }
  
  invisible(list(
    accessions = accessions,
    short_names = final_names,
    raw_fasta = raw_fasta,
    aligned_fasta = aligned_fasta,
    trimmed_fasta = trimmed_fasta,
    trimmed_phylip = trimmed_phylip
  ))
}

# ----------------
# RUN MARKERS
# ----------------

# COI
coi_result <- process_marker(
  marker_name = "COI",
  accessions = coi_accessions,
  manual_name_map = manual_name_map,
  project_dir = project_dir,
  gap_threshold = gap_threshold,
  iqtree_bin = iqtree_bin,
  run_tree = FALSE
)

# Histone H3
if (length(h3_accessions) > 0) {
  h3_result <- process_marker(
    marker_name = "H3",
    accessions = h3_accessions,
    manual_name_map = manual_name_map,
    project_dir = project_dir,
    gap_threshold = gap_threshold,
    iqtree_bin = iqtree_bin,
    run_tree = FALSE
  )
}
