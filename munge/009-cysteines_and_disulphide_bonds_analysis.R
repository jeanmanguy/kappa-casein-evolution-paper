# find all positions in the alignment with at least one cysteine ----
cysteines_msa_positions_interests <- kappa_caseins_alignment_positions %>%
  filter(residue == "C") %>%
  group_by(msa_index) %>%
  pull(msa_index) %>%
  unique()

cysteines_msa_positions_interests_pseudo_seqs <- kappa_caseins_alignment_positions %>%
  filter(msa_index %in% cysteines_msa_positions_interests) %>%
  group_by(species) %>%
  summarise(pseudo_seq = glue_collapse(residue))


# we use positions in the mature kappa-casein sequence as references ----
cysteines_cattle_positions_pivot <- kappa_caseins_alignment_positions %>%
  filter(msa_index %in% cysteines_msa_positions_interests, species == "Bos_taurus") %>%
  select(msa_index, cattle_position = seq_index)

# kappa_caseins_alignment_positions_clean <-


cache("kappa_caseins_cysteines_positions_rel_cattle", {
  kappa_caseins_alignment_positions %>%
    filter(msa_index %in% cysteines_msa_positions_interests, residue != "-") %>%
    select(-seq_index) %>%
    left_join(cysteines_cattle_positions_pivot, by = "msa_index") %>%
    # mutate_at(c("msa_index", "cattle_position"), as.factor) %>%
    left_join(amino.acid.classification.cysteines.plot, by = "residue") %>%
    arrange(cattle_position)
}, depends = c("amino.acid.classification.cysteines.plot", "kappa_caseins_alignment_positions", "cysteines_msa_positions_interests", "cysteines_cattle_positions_pivot"))


# classification dimer / oligomer / intrachain ----

cache("kappa_caseins_SS_bonds_classes", {
  kappa_caseins_alignment_positions %>%
    left_join(cysteines_cattle_positions_pivot, by = "msa_index") %>%
    filter(!is.na(cattle_position)) %>%
    mutate(
      position = as.factor(cattle_position) %>%
        fct_collapse(`10` = "10", `11` = "11", `other` = c("19", "37", "55", "57", "60", "76", "88", "97", "100", "102")) %>%
        as.character()) %>%
    select(species, residue, position) %>%
    group_by(species) %>%
    nest() %>%
    mutate(data = map(data, filter, residue == "C")) %>%
    mutate(C_positions = map(data, pull, position)) %>%
    select(-data) %>%
    mutate(
      intrachain = map_lgl(C_positions, ~ "other" %in% .x),
      dimer = map_int(C_positions, length) > 0,
      oligomer = dimer & map_int(C_positions, length) > 1
    ) %>%
    select(-C_positions) %>%
    gather("bond", "possibility", -species) %>%
    mutate(bond = factor(bond, levels = c("dimer", "oligomer", "intrachain")))
}, depends = c("cysteines_cattle_positions_pivot", "kappa_caseins_alignment_positions"))


