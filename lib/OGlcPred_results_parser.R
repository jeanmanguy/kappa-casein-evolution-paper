parse_OGlcPred_result <- . %>%
  read_html() %>%
  {
    species <- html_nodes(., "h4") %>%
      html_text() %>%
      str_remove(":") %>%
      str_trim()

    result <- html_nodes(., "table") %>%
      html_table(header = TRUE) %>%
      map(as_tibble)

    set_names(result, species) %>%
      bind_rows(.id = "species") %>%
      mutate(Oglyc_pred = PositiveOrNegative == "positive", species = as.character(species)) %>%
      select(-Window, position = Position, -PositiveOrNegative)
  }
