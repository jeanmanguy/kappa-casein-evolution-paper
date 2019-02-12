# load Kyte Doolittle hydrophobicity scale from teh Peptides package
# see `?Peptides::AAdata`

data(AAdata)

cache("KyteDoolittle_scale_df", {
  AAdata$Hydrophobicity$KyteDoolittle %>%
    as.data.frame() %>%
    rownames_to_column("residue") %>%
    as_tibble() %>%
    rename("hydrophobicity" = ".") %>%
    arrange(desc(hydrophobicity))
}, depends = c("AAdata"))


# join data frame with residue positions with hydrophobicity scale ----
kappa_caseins_hydrophobicity <- kappa_caseins_alignment_positions %>%
  left_join(KyteDoolittle_scale_df, by = c("residue"))

kappa_caseins_hydrophobicity_parts <- kappa_caseins_alignment_positions_parts %>%
  left_join(KyteDoolittle_scale_df, by = c("residue"))

# means ----

# mean hydrophobicity mature
mean_hydrophobicity <- kappa_caseins_hydrophobicity %>%
  filter(!is.na(hydrophobicity)) %>%
  group_by(species) %>%
  summarise(
    mean_hydrophobicity = mean(hydrophobicity)
  ) %>%
  ungroup()

# mean hydrophobicity parts
mean_hydrophobicity_parts <- kappa_caseins_hydrophobicity_parts %>%
  group_by(species, part) %>%
  summarise(
    mean_hydrophobicity = mean(hydrophobicity)
  ) %>%
  ungroup()

