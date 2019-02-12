# Find best model ----
cache("model_test", {
  phangorn::modelTest(alignment_phyDat, tree = trimmed_species_tree, model = "all", multicore = TRUE, mc.cores = 4)
}, depends = c("trimmed_species_tree", "alignment_phyDat"))

# select best model selection ----
best_model <- model_test$Model[which.min(model_test$AIC)]
env <- attr(model_test, "env") # contains trees and model parameters

# start fit ----
fit_start <- eval(get(best_model, env), env)

cache("fit", {
  phangorn::optim.pml(
    fit_start,
    model = str_split(best_model, "\\+", simplify = TRUE)[[1]],
    optBf = TRUE,
    optEdge = TRUE,
    optGamma = FALSE,
    optRate = FALSE,
    optQ = FALSE,
    optRooted = FALSE,
    optInv = TRUE,
    rearrangement = "none"
  )
}, depends = "fit_start")


# rooting of the fitted tree ----
fitted_tree <- fit$tree %>%
  ape::root(outgroup = c("Tachyglossus_aculeatus", "Ornithorhynchus_anatinus")) %>%
  phangorn::midpoint(node.labels = "label")

# save trees ----
ape::write.tree(phy = fitted_tree, file = "cache/kappa-casein_fitted_tree_alignment.nwk")
ape::write.nexus(phy = fitted_tree, file = "cache/kappa-casein_fitted_tree_alignment.nxs")


# comparison protein and phylo tree ----
cache("dist_tree_phylo", {
  trimmed_species_tree %>%
    adephylo::distTips() %>%
    tidy_dist()
}, depends = "trimmed_species_tree")


cache("dist_tree_fitted", {
  fit$tree %>%
    adephylo::distTips() %>%
    tidy_dist()
}, depends = "fit")


cache("dist_align", {
  as.AAbin(kappa_casein_alignment) %>%
    dist.aa(scaled = TRUE, pairwise.deletion = TRUE) %>%
    tidy_dist()
}, depends = "kappa_casein_alignment")


cache("compare_evodist_divtime", {
  bind_rows(
    `TRUE` = dist_tree_fitted,
    `FALSE` = dist_align,
    .id = "corrected"
  ) %>%
    inner_join(dist_tree_phylo, by = c("species.1", "species.2"), suffix = c(".fit", ".phylo")) %>%
    mutate(id = map2(species.1, species.2, c) %>% map(sort) %>% map_chr(glue_collapse, sep = "-")) %>% # remove identical measures (human - cow, and cow - human for example)
    group_by(id, corrected) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(-id) %>%
    mutate(
      divergence.time = distance.phylo / 2,
      corrected = as.logical(corrected)
    )
}, depends = c("dist_align", "dist_tree_fitted", "dist_tree_phylo"))


