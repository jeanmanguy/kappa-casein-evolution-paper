# use negative values for PKC (at the left/N-side of the cleavage site)
parts_positions_neg <- kappa_casein_parts_startsends %>%
  mutate(position = ifelse(part == "PKC", -end, end - start)) %>%
  order_species_tree()

# adjust position of PTRs for the plot
indels_repeats_adjusted_cleavage <- kappa.casein.manual.repeats %>%
  select(-sequence) %>%
  mutate(type = "repeat") %>%
  mutate(id = row_number()) %>%
  gather("side", "position", ends_with("position")) %>%
  left_join(kappa_casein_parts_startsends %>% filter(part == "PKC") %>% mutate_at("species", as.character), by = c("species")) %>%
  mutate(adjusted_position = ifelse(position < end, -(end - position), position - end))


df_plot_indels <- indels_repeats_adjusted_cleavage %>%
  select(id, species, type, side, adjusted_position) %>%
  spread(side, adjusted_position) %>%
  full_join(trimmed_species_tree_order_pivot, by = "species") %>%
  select(-y) %>%
  order_species_tree() %>%
  mutate(type = "tandem repeat") %>%
  filter(!is.na(id))

