# Gerstein–Sonnhammer–Chothia weight
# Gerstein, M., Sonnhammer, E. L., and Chothia, C. (1994), ‘Volume Changes in Protein Evolution’, Journal of Molecular Biology, 236/4 (4 Mar.): 1067–78. doi: 10.1016/0022-2836(94)90012-4.
# using the aphid package
# Wilkinson SP (2017) The 'aphid' package for analysis with profile hidden Markov models. R package version 1.0.0.
# https://cran.r-project.org/package=aphid


weight_phylo_ <- function(tree) {
  tree %>% ape::as.hclust.phylo() %>% as.dendrogram() %>% aphid::weight(method = "Gerstein")
}

weight_phylo <- function(trees) {
  if (class(trees) == "phylo") {
    trees %>%
      weight_phylo_()
  } else {
    trees %>%
      map(weight_phylo_)
  }
}
