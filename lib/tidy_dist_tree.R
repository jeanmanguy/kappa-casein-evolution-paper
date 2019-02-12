# take a dist object and return a tidy data frame with 2 columns for the 2 species and 1 column for the distance value
tidy_dist <- . %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "species.1") %>%
  gather("species.2", "distance", -species.1) %>%
  as_tibble() %>%
  filter(species.1 != species.2)
