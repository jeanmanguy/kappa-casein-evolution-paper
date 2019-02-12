# The result of the predictions ofr each tool was saved as html or text file and need to be parsed

# netoglyc ----
cache("kappa_casein_netoglyc_predictions", {
  list_oglyc_pred_paths$netoglyc %>%
    parse_netoglyc_result() %>%
    order_species_tree()
}, depends = c("list_oglyc_pred_paths", "trimmed_species_tree_order_pivot"))

# GlycoMine ----
cache("kappa_casein_glycomine_predictions", {
  list_oglyc_pred_paths$glycomine %>%
    parse_glycomine_result() %>%
    mutate(Oglyc_pred = score > 0.5) %>%
    order_species_tree()
}, depends = c("list_oglyc_pred_paths", "trimmed_species_tree_order_pivot"))

# OGlcPred ----
cache("kappa_casein_OGlcPred_predictions", {
  list_oglyc_pred_paths$OglcPred %>%
    parse_OGlcPred_result() %>%
    order_species_tree()
}, depends = c("list_oglyc_pred_paths", "trimmed_species_tree_order_pivot"))


# GlycoPred ----
cache("kappa_casein_glycopred_predictions", {
  list_oglyc_pred_paths$glycopred %>%
    parse_glycopred_result() %>%
    order_species_tree()
}, depends = c("list_oglyc_pred_paths", "trimmed_species_tree_order_pivot"))


# cow and human known glycosylation data/known_oglyc_kappa-casein.csv
# source UniCarbKB
# http://www.unicarbkb.org/proteinsummary/P07498/annotated
# http://www.unicarbkb.org/proteinsummary/P02668/annotated
cache("known_oglyc_position_df", {
  known.oglyc.kappa.casein %>%
    mutate(
      species = str_replace(species, " ", "_"),
      position = str_split(position, ",") %>% map(str_remove, pattern = " ")
    ) %>%
    unnest() %>%
    separate(position, into = c("residue3", "position")) %>%
    mutate( # substract the length of the signal peptides to the position (Cow: 21aa, Human: 20 aa)
      position = if_else(species == "Bos_taurus", as.integer(position) - 21L, as.integer(position) - 20L),
      residue = if_else(`residue3` == "THR", "T", "S")
    ) %>%
    select(-residue3) %>%
    left_join(kappa_caseins_alignment_positions, by = c("species", "position" = "seq_index", "residue"))
}, depends = c("known.oglyc.kappa.casein", "kappa_caseins_alignment_positions"))


known_oglyc_position_grouped_df <- known_oglyc_position_df %>%
  group_by(msa_index) %>%
  summarise(label = glue_collapse(common_name, " and "))


# merge predictions ----
cache("kappa_casein_all_predictions", {
  bind_rows(
    netoglyc = kappa_casein_netoglyc_predictions %>% mutate_at("species", as.character),
    `O-GlcNAcPRED-II` = kappa_casein_OGlcPred_predictions %>% mutate_at("species", as.character),
    glycomine = kappa_casein_glycomine_predictions %>% mutate_at("species", as.character),
    glycopred = kappa_casein_glycopred_predictions %>% mutate_at("species", as.character),
    .id = "method"
  ) %>%
    order_species_tree() %>%
    left_join(kappa_caseins_alignment_positions, by = c("species", "position" = "seq_index"))
}, depends = c("kappa_casein_netoglyc_predictions", "kappa_casein_OGlcPred_predictions", "kappa_casein_glycomine_predictions", "kappa_casein_glycopred_predictions", "kappa_caseins_alignment_positions"))


# get all serine/threonine in cow and human ----
all_ST_known_oglyc <- kappa_caseins_alignment_positions %>%
  filter( # negative: S and T of human and cow
    residue %in% c("S", "T"),
    species %in% c("Bos_taurus", "Homo_sapiens")
  ) %>%
  mutate_at("species", as.character) %>%
  left_join(known_oglyc_position_df, by = c("species", "msa_index", "residue")) %>%
  mutate(known_oglyc_site = !is.na(position)) %>%
  select(-common_name, -position)


# determination of True positive, false positive, true negative, false negative ----
# for each method for each sequence

mutate_fptpfntn <- . %>%
  mutate(
      true_positive = known_oglyc_site == TRUE & Oglyc_pred == TRUE,
      true_negative = known_oglyc_site == FALSE & Oglyc_pred == FALSE,
      false_positive = known_oglyc_site == FALSE & Oglyc_pred == TRUE,
      false_negative = known_oglyc_site == TRUE & Oglyc_pred == FALSE,
    )

cache("all_predictions_fptpfntn", {
  kappa_casein_all_predictions %>%
    filter(species %in% c("Bos_taurus", "Homo_sapiens")) %>%
    mutate_at("species", as.character) %>%
    left_join(all_ST_known_oglyc, by = c("species", "msa_index", "residue")) %>%
    select(species, msa_index, method, known_oglyc_site, Oglyc_pred) %>%
    mutate_fptpfntn() %>%
    select(-species, -msa_index, -known_oglyc_site, -Oglyc_pred) %>%
    gather("state", "value", true_positive, true_negative, false_positive, false_negative) %>%
    filter(value) %>%
    group_by(method, state) %>%
    tally() %>%
    ungroup() %>%
    complete(method, state, fill = list(n = 0)) %>%
    spread(state, n)
}, depends = c("kappa_casein_all_predictions", "all_ST_known_oglyc"))


tranmute_metrics_fptpfntn <- . %>%
  transmute(
      method,
      sensivity = true_positive / (true_positive + false_negative),
      specificity = true_negative / (true_negative + false_positive),
      precision = true_positive / (true_positive + false_positive),
      accuracy = (true_positive + true_negative) / (true_positive + true_negative + false_positive + false_negative),
      MCC = (true_positive * true_negative - false_positive * false_negative) / sqrt((true_positive + false_positive) * (true_positive + false_negative) * (true_negative + false_positive) * (true_negative + false_negative))
    )

cache("all_predictions_metrics", {
  all_predictions_fptpfntn %>%
    tranmute_metrics_fptpfntn() %>%
    arrange(desc(MCC))
}, depends = c("all_predictions_fptpfntn"))


# selected method ----
# Glycomine with a lower threshold
cache("kappa_casein_oglyc_predictions", {
  list_oglyc_pred_paths$glycomine %>%
    parse_glycomine_result() %>%
    mutate(Oglyc_pred = score > 0.4) %>%
    order_species_tree()
}, depends = c("list_oglyc_pred_paths", "trimmed_species_tree_order_pivot"))

# metrics for this method and threshold

# get position in the alignment of the predictions
cache("kappa_casein_oglyc_predictions_msa", {
  kappa_casein_oglyc_predictions %>%
  right_join(kappa_caseins_alignment_positions, by = c("species", "position" = "seq_index")) %>%
  filter(Oglyc_pred) %>%
  left_join(parts_kappa_casein_alignment, by = c("msa_index")) %>%
  order_cleavage_parts() %>%
  order_species_tree()
}, depends = c("kappa_casein_oglyc_predictions", "trimmed_species_tree_order_pivot", "parts_kappa_casein_alignment", "kappa_caseins_alignment_positions"))

# count the number of O-glycosylation in each part in each species
kappa_casein_oglyc_glycomine_counts <- kappa_casein_oglyc_predictions_msa %>%
  filter(Oglyc_pred) %>%
  group_by(species, part) %>%
  tally() %>%
  ungroup() %>%
  complete(species, part, fill = list(n = 0)) # add zero counts
