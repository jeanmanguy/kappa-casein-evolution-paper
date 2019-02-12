# Pepsin cleavage regex ----
# see lib/globals.R
pepsin_cleavage_site <- partial(cleaver::cleavageSites, enzym = pepsin_enzym)
plasmin_cleavage_site <- partial(cleaver::cleavageSites, custom = plasmin_regex)


get_cleavage_positions <- . %>%
  mutate(
    pepsin = pepsin_cleavage_site(subseq),
    plasmin = plasmin_cleavage_site(subseq)
  ) %>%
  gather("peptidase", "cleavage_position", pepsin, plasmin) %>%
  select(-subseq)

# Find cleavage sites in each parts of kappa-casein ----
cache("kappa_casein_parts_cleavage_positions", {
  kappa_caseins_splitted_cleavage %>%
    get_cleavage_positions()
}, depends = c("kappa_caseins_splitted_cleavage", "pepsin_cleavage_site", "plasmin_cleavage_site"))

cache("mature_cleavage_sites", {
  kappa_caseins_mature %>%
    get_cleavage_positions()
}, depends = c("kappa_caseins_mature", "pepsin_cleavage_site", "plasmin_cleavage_site"))

# Cleavage sites count in each parts of kappa-casein for each peptidase ----
# count and frequency
cache("kappa_caseins_splitted_cleavage_peptidases_counts", {
  kappa_casein_parts_cleavage_positions %>%
    left_join(kappa_casein_parts_lengths, by = c("species", "part")) %>%
    mutate(
      n = map_int(cleavage_position, base::length),
      freq = n / (length - 1) # length - 1: number of peptidic bonds
    ) %>%
    select(-cleavage_position, -length) %>%
    order_species_tree()
}, depends = c("kappa_casein_parts_cleavage_positions", "kappa_casein_parts_lengths", "trimmed_species_tree_order_pivot"))

# Cleavage sites count in each parts of kappa-casein for all peptidases ----
cache("kappa_caseins_splitted_cleavage_all_peptidases_counts", {
  kappa_casein_parts_cleavage_positions %>%
    unnest() %>%
    select(-peptidase) %>%
    distinct() %>% # remove cleavage sites shared by 2+ peptidases
    group_by(species, part) %>%
    tally() %>%
    ungroup() %>%
    left_join(kappa_casein_parts_lengths, by = c("species", "part")) %>%
    transmute(species, part, frequency = n / (length - 1))
}, depends = c("kappa_casein_parts_cleavage_positions"))

# Tidy cleavage site positions ----
cache("kappa_caseins_cleavage_positions", {
  mature_cleavage_sites %>%
    unnest() %>%
    order_species_tree() %>%
    arrange(species, cleavage_position)
}, depends = c("mature_cleavage_sites", "trimmed_species_tree_order_pivot"))
