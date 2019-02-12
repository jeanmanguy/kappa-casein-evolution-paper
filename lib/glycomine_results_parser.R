parse_glycomine_result <- . %>%
  read_html() %>%
  {
    species <- html_nodes(., "h2") %>%
      html_text() %>%
      str_match(pattern = "^>([A-Z][a-z]*_[a-z]*)") %>%
      .[, 2]

    result <- html_nodes(., "table") %>%
      html_table(header = TRUE) %>%
      map(as_tibble)

    set_names(result, species) %>%
      bind_rows(.id = "species")
  } %>%
    select(species, position = Position, score = Score) %>%
    filter(position != "Position") %>%
    mutate_at("position", as.integer) %>%
    mutate_at("score", as.double)


