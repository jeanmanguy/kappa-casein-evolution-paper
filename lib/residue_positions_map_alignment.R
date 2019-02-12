locate_non_gap <- . %>%
  map(str_locate_all, "[^-]") %>%
  map(purrr::pluck, 1)

map_seq_indexes_to_msa <- . %>%
  as.character() %>% # transform the alignment into a list of characters
  locate_non_gap() %>% # find all the non-gap positions (return a matrix with start and end of each match)
  map(as_tibble) %>% # transform the list of matrices to a list of tibbles
  map(select, "start") %>% # remove the end column, not necessary here
  bind_rows(.id = "species") %>% # merge the list of data_frames, and add a new column 'species' filed with the name of each data_frame
  rename("msa_index" = "start") %>% # we rename the start column
  group_by(species) %>%
  mutate(seq_index = 1:n()) %>% # create the sequence index for each species/protein
  ungroup()

get_residue_msa_index <- . %>%
  Biostrings::AAStringSet() %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  rowwise() %>%
  transmute(species, residue = str_split(x, "")) %>%
  unnest() %>%
  group_by(species) %>%
  mutate(msa_index = 1:n()) %>%
  ungroup()


get_residue_msa_seq_index <- . %>% {
  full_join(
    get_residue_msa_index(.),
    map_seq_indexes_to_msa(.),
    by = c("msa_index", "species")
  )
}
