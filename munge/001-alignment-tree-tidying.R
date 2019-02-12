# trim the mammalian tree to only use the species with a sequence ----

cache("trimmed_species_tree", {
  nodes_to_drop <- mammals_tree %>%
    ggplot2::fortify() %>%
    filter(isTip) %>%
    mutate(label = unlist(label)) %>%
    anti_join(kappa_casein, by = c("label" = "species")) %>%
    pull(node)

  tree_trimmed <- mammals_tree %>%
    drop.tip(nodes_to_drop) %>%
    multi2di() # fix multifurcating

  # replace branch length == 0 by very small value
  tree_trimmed$edge.length[tree_trimmed$edge.length == 0] <- 1e-08

  tree_trimmed$tip.label <- unlist(tree_trimmed$tip.label)
  tree_trimmed
}, depends = c("mammals_tree", "kappa_casein"))


write.nexus(trimmed_species_tree, file = "cache/trimmed_species_tree.nexus")


# tree as a tibble ----
cache("trimmed_species_tree_as_df", {
  trimmed_species_tree %>%
    ggplot2::fortify() %>%
    as_tibble()
}, depends = c("trimmed_species_tree"))


# order sequences according to the tree ----
## pivot table
cache("trimmed_species_tree_order_pivot", {
  trimmed_species_tree_as_df %>%
    filter(isTip) %>%
    mutate(label = unlist(label), y = as.integer(y)) %>%
    select(species = label, y) %>%
    arrange(y)
}, depends = c("trimmed_species_tree_as_df"))

## function working on the "species" field
order_species_tree <- . %>%
  mutate(species = as.character(species)) %>%
  left_join(trimmed_species_tree_order_pivot, by = "species") %>%
  arrange(y) %>%
  mutate(species = fct_reorder(species, y)) %>%
  select(-y)

# get table with amino acid position in the sequence and in the alignment ----
cache("kappa_caseins_alignment_positions", {
  kappa_casein_alignment %>%
    get_residue_msa_seq_index() %>%
    order_species_tree()
}, depends = c("kappa_casein_alignment", "trimmed_species_tree_order_pivot"))



# compute weights based on the tree ----
kappa_caseins_phylo_tree_weights <- trimmed_species_tree %>%
  weight_phylo() %>%
  {
    tibble(species = names(.), weight = .)
  } %>%
  mutate(weight_norm = weight / sum(weight)) %>%
  order_species_tree()


# slimmer sequence tibble ----
kappa_casein_df <- kappa_casein %>%
  select(species, sequence, length)

kappa_casein_gi <- kappa_casein %>%
  select(species, taxid, gi) %>%
  order_species_tree()

# get residue position in the alignment for cattle ----
msa_positions_cattle_only <- kappa_caseins_alignment_positions %>%
  filter(species == "Bos_taurus", residue != "-") %>%
  select(msa_index, seq_index)

# split kappa-casein between PKC and GMP ----
# using known positions in cattle

# fonction to order PKC and GMP
order_cleavage_parts <- . %>%
  mutate(part = as.factor(part) %>% fct_relevel(c("PKC", "GMP")))


# get the end of the signal peptide
pos_sigpept <- kappa.casein.parts %>% filter(part == "signal peptide") %>% pull(end)

cache("parts_kappa_casein_alignment", {
  kappa.casein.parts %>%
    filter(part != "signal peptide") %>%
    select(-species) %>%
    gather("position", "index", start, end) %>%
    mutate(index = index - pos_sigpept) %>%
    left_join(msa_positions_cattle_only, by = c("index" = "seq_index")) %>%
    select(-index) %>%
    spread(position, msa_index) %>%
    rowwise() %>%
    mutate(msa_index = list(start:end)) %>%
    select(part, msa_index) %>%
    unnest() %>%
    ungroup() %>%
    replace_na(list(part = "GMP")) %>% # fix bug with sequences with a few extra aa at the C terminus
    order_cleavage_parts() %>%
    arrange(part)
}, depends = c("kappa.casein.parts", "msa_positions_cattle_only"))



# residue at the N of the chymosin cleavage site
position_align_chymosin_cleavage_site <- parts_kappa_casein_alignment %>%
  filter(part == "PKC") %>%
  pull(msa_index) %>%
  max()

cache("kappa_caseins_mature", {
  kappa_caseins_alignment_positions %>%
    filter(!is.na(seq_index)) %>%
    group_by(species) %>%
    summarise(subseq = glue_collapse(residue) %>% as.character()) %>%
    order_species_tree()
}, depends = c("kappa_caseins_alignment_positions", "trimmed_species_tree_order_pivot"))

cache("kappa_caseins_mature_lengths", {
  kappa_caseins_mature %>%
    transmute(species, length = nchar(subseq))
}, depends = c("kappa_caseins_mature"))

# get dataframe with for each species the position in the sequence, in the alignment and appartenance to PKC or GMP
cache("kappa_caseins_alignment_positions_parts", {
  kappa_caseins_alignment_positions %>%
    mutate(part = if_else(msa_index < position_align_chymosin_cleavage_site + 1, "PKC", "GMP")) %>%
    filter(residue != "-")
}, depends = c("kappa_caseins_alignment_positions"))

# the sequence of the PKC and GMP for each species ---
cache("kappa_caseins_splitted_cleavage", {
  kappa_caseins_alignment_positions_parts %>%
    group_by(species, part) %>%
    summarise(subseq = glue_collapse(residue) %>% as.character()) %>%
    ungroup() %>%
    order_species_tree() %>%
    order_cleavage_parts() %>%
    arrange(part)
}, depends = c("kappa_caseins_alignment_positions_parts", "trimmed_species_tree_order_pivot"))


# merged data frame with mature sequences and both parts
kappa_caseins_parts_mature <- kappa_caseins_mature %>%
  mutate(part = "mature") %>%
  bind_rows(kappa_caseins_splitted_cleavage %>% mutate(part = as.character(part)))


# start and end positions of each part of kappa-casein for each species ----
cache("kappa_casein_parts_startsends", {
  kappa_caseins_alignment_positions_parts %>%
    group_by(species, part) %>%
    summarise(
      start = min(seq_index) %>% as.integer(),
      end = max(seq_index) %>% as.integer()) %>%
    ungroup()
}, depends = c("kappa_caseins_alignment_positions_parts"))


# kappa-casein parts lengths for each species ----
cache("kappa_casein_parts_lengths", {
  kappa_casein_parts_startsends %>%
    transmute(species, part, length = end - start + 1L) %>%
    ungroup() %>%
    order_cleavage_parts() %>%
    arrange(desc(part), species)
}, depends = c("kappa_casein_parts_startsends"))


# parts average length (weighted) ----
kappa_casein_parts_lengths_avg <- kappa_casein_parts_lengths %>%
  left_join(kappa_caseins_phylo_tree_weights, by = "species") %>%
  group_by(part) %>%
  summarise(
    w_mean = Hmisc::wtd.mean(length, weights = weight),
    w_var = Hmisc::wtd.var(length, weights = weight),
    w_sd = sqrt(w_var)
  )

# alignment in another format ----
alignment_phyDat <- phangorn::as.phyDat(kappa_casein_alignment, type = "AA")

# write fasta files ----
seqinr::write.fasta(sequences = as.list(kappa_caseins_mature$subseq), names = kappa_caseins_mature$species, file.out = "cache/kappa_caseins_mature.faa")
kappa_caseins_splitted_cleavage %>%
  filter(part == "GMP") %>%
  {
    seqinr::write.fasta(sequences = as.list(.$subseq), names = .$species, file.out = "cache/kappa_caseins_GMP.faa")
  }
kappa_caseins_splitted_cleavage %>%
  filter(part == "PKC") %>%
  {
    seqinr::write.fasta(sequences = as.list(.$subseq), names = .$species, file.out = "cache/kappa_caseins_PKC.faa")
  }
