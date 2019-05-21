# simple count of amino acids
# absent residues have n = 0
compute_aa_composition <- . %>%
  rowwise() %>%
  mutate(
    n = list(str_count(subseq, residues)),
    residues = list(residues),
    length = nchar(subseq)
  ) %>%
  ungroup() %>%
  select(-subseq) %>%
  unnest(residues, n) %>%
  mutate(
    freq = n / length
  ) %>%
  select(-length) %>%
  order_species_tree()


cache("kappa_caseins_splitted_cleavage_aa_composition", {
  kappa_caseins_splitted_cleavage %>%
    compute_aa_composition()
}, depends = c("kappa_caseins_splitted_cleavage"))


cache("kappa_caseins_parts_aa_composition", {
  kappa_caseins_splitted_cleavage_aa_composition %>%
    select(-n) %>%
    spread(part, freq) %>%
    mutate(residues = as.factor(residues) %>% fct_reorder(-PKC)) # plot ordered by descending amino acid frequency in PKC
}, depends = c("kappa_caseins_splitted_cleavage_aa_composition"))


cache("kappa_caseins_mature_aa_composition", {
  kappa_caseins_mature %>%
    compute_aa_composition()
}, depends = c("kappa_caseins_mature"))


kappa_caseins_aa_comp_df_plot <- kappa_caseins_splitted_cleavage_aa_composition %>%
  select(-n) %>%
  group_by(part, residues) %>%
  summarise_at(vars(freq), .funs = list(min = min, max = max)) %>%
  ungroup() %>%
  mutate(residues = fct_relevel(residues, list_residues_ordered_aa_comp))
