get_protein_sequences_from_xml <- . %>%
  xmlParse(file = .) %>%
  {
    tibble(
      gi = xpathSApply(., "/TSeqs/TSeq/TSeq_gi", xmlValue),
      taxid = xpathSApply(., "/TSeqs/TSeq/TSeq_taxid", xmlValue),
      species = xpathSApply(., "/TSeqs/TSeq/TSeq_orgname", xmlValue),
      definition = xpathSApply(., "/TSeqs/TSeq/TSeq_defline", xmlValue),
      length = xpathSApply(., "/TSeqs/TSeq/TSeq_length", xmlValue) %>% as.integer(),
      sequence = xpathSApply(., "/TSeqs/TSeq/TSeq_sequence", xmlValue)
    )
  }

remove_subspecies_name <- . %>%
  mutate(species = str_replace(species, "^([A-Z][a-z]+ [a-z]+) [a-z]+", "\\1"))

filter_longest_protein_per_species <- . %>%
  group_by(species) %>%
  arrange(species, desc(length)) %>%
  dplyr::slice(1) %>%
  ungroup()


# manually reconstructed sequence of manatee
manually_reconstructed_sequence <- tribble(
  ~taxid, ~species, ~definition, ~sequence,
  "127582", "Trichechus_manatus", "kappa-casein [Trichechus manatus latirostris] with 1 exon recovered", Biostrings::readAAStringSet(filepath = "data/NW_004444044.1[4558376..4558432][4556874..4556906][4555358..4555819](-).faa") %>% as.character() %>% .[[1]]
)

# load the sequences from the XML files
# clean the names
# save a cached version
cache("kappa_casein", {
  file.path("data", "kappa-casein.xml") %>%
    get_protein_sequences_from_xml() %>%
    filter_longest_protein_per_species() %>%
    remove_subspecies_name() %>%
    mutate(species = str_replace(species, " ", "_")) %>%
    bind_rows(manually_reconstructed_sequence)
}, depends = c("manually_reconstructed_sequence"))
