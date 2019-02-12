netoglyc_name_list <- . %>%
  set_names(c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "comment"))

netoglyc_as_df <- . %>% {
  tibble(
    seqname = map_chr(., "seqname") %>% unlist(),
    source = map_chr(., "source") %>% unlist(),
    feature = map_chr(., "feature") %>% unlist(),
    start = map_chr(., "start") %>% as.integer() %>% unlist(),
    end = map_chr(., "end") %>% as.integer() %>% unlist(),
    score = map_chr(., "score") %>% as.double() %>% unlist(),
    strand = map_chr(., "strand") %>% unlist(),
    frame = map_chr(., "frame") %>% unlist(),
    comment = map_chr(., "comment") %>% unlist()
  )
}

parse_netoglyc_result <- . %>%
  map_df(~ .x %>%
    read_html() %>%
    html_text("pre") %>%
    str_split("\n") %>%
    map(str_split, pattern = "\t") %>%
    unlist(recursive = FALSE) %>%
    .[25:(length(.) - 5)] %>% # remove header and footer
    map(netoglyc_name_list) %>%
    netoglyc_as_df()) %>%
  transmute(
    species = seqname %>% str_to_title(locale = "en"), # lowercase except the first letter
    position = start, # only one residue
    score,
    Oglyc_pred = comment == "#POSITIVE"
  )
