parse_glycopred_result <- . %>% 
  map(~ {
    read_tsv(file = ., col_names = c("position", "residue", "Oglyc_pred"), col_types = "icc", n_max = R.utils::countLines(.) - 4L)
  }) %>% 
  bind_rows(.id = "id") %>% 
  group_by(id) %>% 
  mutate(protein_sequence = paste(residue, collapse = "")) %>% 
  left_join(kappa_caseins_mature, by = c("protein_sequence" = "subseq")) %>% 
  ungroup() %>% 
  select(-id, -protein_sequence) %>% 
  filter(residue %in% c("S", "T")) %>% 
  mutate(Oglyc_pred = Oglyc_pred == "G") %>% 
  select(-residue)




