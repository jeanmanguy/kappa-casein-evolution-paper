# phosphorylation by fam20c

# find all matches with the regular expression ----
cache("kappa_casein_fam20c_matches", {
  kappa_caseins_splitted_cleavage %>%
    mutate(
      match_position = subseq %>% str_locate_all(fam20c_motif) %>% map(~ pull(as_tibble(.), "start"))
    ) %>%
    select(-subseq)
}, depends = c("kappa_caseins_splitted_cleavage", "fam20c_motif"))

# count matches ----
cache("kappa_casein_fam20c_counts", {
  kappa_casein_fam20c_matches %>%
    mutate(
      n = map_int(match_position, ~ length(.x))
    ) %>%
    select(-match_position)
}, depends = c("kappa_casein_fam20c_matches"))

# get position of matches (relative to the beginning of the mature proteins) -----
cache("kappa_casein_fam20c_matches_tidy", {
  kappa_casein_fam20c_matches %>%
    mutate_at(c("species", "part"), as.character) %>%
    left_join(kappa_casein_parts_startsends %>% mutate_at(c("species", "part"), as.character), by = c("species", "part")) %>%
    unnest() %>%
    mutate(match_position = if_else(part == "PKC", match_position, match_position + (start - 1L))) %>%
    full_join(trimmed_species_tree_order_pivot, by = "species") %>%
    select(-part, -start, -end, -y) %>%
    order_species_tree()
}, depends = c("kappa_casein_fam20c_matches", "kappa_casein_parts_lengths"))



