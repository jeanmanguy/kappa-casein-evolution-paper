# Load supertree and correct some species names to match the NCBI taxonomy
# Fritz, S. A., Bininda-Emonds, O. R. P., and Purvis, A. (2009), ‘Geographical Variation in Predictors of Mammalian Extinction Risk: Big Is Bad, but Only in the Tropics’, Ecology Letters, 12/6 (1 June): 538–49. doi: 10.1111/j.1461-0248.2009.01307.x

cache("mammals_tree", {
  mammals_tree_ <- read.nexus(file.path("data", "ELE_1307_sm_SA1.tre"))$mammalST_MSW05_bestDates
  mammals_tree_$tip.label <- mammals_tree_$tip.label %>%
    map_chr(~ switch(.x,
      "Spermophilus_tridecemlineatus" = "Ictidomys_tridecemlineatus",
      "Cryptomys_damarensis" = "Fukomys_damarensis",
      "Vicugna_vicugna" = "Vicugna_pacos",
      "Spalax_ehrenbergi" = "Nannospalax_galili",
      "Bos_grunniens" = "Bos_mutus",
      "Monachus_schauinslandi" = "Neomonachus_schauinslandi",
      "Tarsius_syrichta" = "Carlito_syrichta",
      .x
    ))
  mammals_tree_
})
