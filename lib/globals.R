# Configuration ----
add.config(
  iupred2_path = file.path("lib", "iupred2a", "iupred2a.py"),
  charge_scale = "Bjellqvist" # "Bjellqvist", "Dawson", "EMBOSS", "Lehninger", "Murray", "Rodwell", "Sillero", "Solomon" or "Stryer"
)

# Amino acids ----
# list_residues_ordered <- residues <- c("S", "T", "G", "A", "P", "I", "L", "V", "M", "F", "H", "Y", "W", "C", "K", "R", "D", "E", "N", "Q")

list_residues_ordered_aa_comp <- residues <- "NQSTILVAMYFWHRKDEGPC" %>%
  str_split("") %>%
  unlist()

one_2_three_aas <- c(
  "A" = "Alanine",
  "R" = "Arginine",
  "N" = "Asparagine",
  "D" = "Aspartic acid",
  "C" = "Cysteine",
  "E" = "Glutamic acid",
  "Q" = "Glutamine",
  "G" = "Glycine",
  "H" = "Histidine",
  "I" = "Isoleucine",
  "L" = "Leucine",
  "K" = "Lysine",
  "M" = "Methionine",
  "F" = "Phenylalanine",
  "P" = "Proline",
  "S" = "Serine",
  "T" = "Threonine",
  "W" = "Tryptophan",
  "Y" = "Tyrosine",
  "V" = "Valine"
)

# peptidases ----
plasmin_regex <- "[KR]"
pepsin_enzym <- "pepsin"

# Fam20C cannonical motif ----
fam20c_motif <- regex("
(?<=E.(S.){0,10})S
| #OR
S(?=(.S){0,10}.E)
", comments = TRUE)

# Oglyc predictions paths
list_oglyc_pred_paths <- list(
  "netoglyc" = list.files("data",  pattern = "^netoglyc_[0-9A-Za-z]+.html$") %>%
    file.path("data", .),
  "glycomine" = list.files("data", pattern = "^GlycoMine.pl.html$") %>%
    file.path("data", .),
  "OglcPred" = list.files("data", pattern = "^OGlcPred.html$") %>%
    file.path("data", .),
  "glycopred" = list.files("data/GlycoPred_output/",  pattern = "^sequence[0-9]+.csv$") %>%
    file.path("data/GlycoPred_output", .)
)

assert_that(length(list_oglyc_pred_paths$netoglyc) > 0)
assert_that(length(list_oglyc_pred_paths$glycomine) > 0)
assert_that(length(list_oglyc_pred_paths$OglcPred) > 0)
assert_that(length(list_oglyc_pred_paths$glycopred) > 0)
