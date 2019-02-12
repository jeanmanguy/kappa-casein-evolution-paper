# IUPRED
# download IUPRED2A and compile it in lib/iupred2a

# compute disorder for each mature sequence
cache("kappa_casein_disorder", {
  kappa_caseins_mature %>%
    mutate(
      seq_index = map(subseq, ~seq(1, nchar(.x))),
      iupred_short = map(subseq, compute_iupred_short)
    ) %>%
    select(-subseq) %>%
    unnest()
}, depends = c("kappa_caseins_mature"))

# join with the positions in the alignment and the species weights
kappa_casein_disorder_align <- kappa_casein_disorder %>%
  right_join(kappa_caseins_alignment_positions, by = c("species", "seq_index")) %>% # join in this order to create NAs to create gaps in the line drawing to ggplot
  order_species_tree()
