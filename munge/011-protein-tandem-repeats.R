# use negative values for PKC (at the left/N-side of the cleavage site)
parts_positions_neg <- kappa_casein_parts_startsends %>%
  mutate(position = ifelse(part == "PKC", -end, end - start)) %>%
  order_species_tree()

df_plot_indels <- kappa.casein.manual.repeats %>%
  as_tibble() %>%
  select(-sequence) %>%
  mutate(id = row_number()) %>%
  gather("side", "position", ends_with("position")) %>%
  left_join(kappa_caseins_alignment_positions, by = c("species", "position" = "seq_index")) %>%
  select(-position, -residue) %>%
  spread(side, msa_index) %>%
  mutate(type = "tandem repeat") %>%
  full_join(trimmed_species_tree_order_pivot, by = "species") %>%
  select(-y) %>%
  order_species_tree() %>%
  filter(!is.na(type))

