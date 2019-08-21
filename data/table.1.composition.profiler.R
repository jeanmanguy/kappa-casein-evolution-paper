table.1.composition.profiler <- read_tsv(
  "data/table_1_composition_profiler.csv",
  col_types = cols(
    Residue = col_character(),
    SwissProt = col_character(),
    PDB.S25 = col_character(),
    Surface.Residues = col_character(),
    DisProt = col_character()
  )
) %>%
  separate(Residue, into = c("residue_3", "residue"), sep = " ") %>%
  mutate_at(vars(residue), str_remove_all, "[()]") %>%
  gather("dataset", "average_sd", -starts_with("residue")) %>%
  separate(average_sd, into = c("average", "sd"), sep = " ± ") %>%
  mutate_at(vars(average, sd), as.numeric)
