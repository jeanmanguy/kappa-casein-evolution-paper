parse_iupred <- . %>%
  .[str_detect(., "^[^#]")] %>% # remove comments
  str_trim() %>%
  str_replace_all("     ", " ") %>% # replace 5 spaces with only 1
  str_extract("([0-9.]+)$") %>%
  as.numeric()

# sequences are submitted one by one to iupred
compute_iupred <- function(sequence, iupred_mode = "long") {
  with_tempfile("iupred_temp_fasta", {
    seqinr::write.fasta(sequence, names = "sequence_tmp", file.out = iupred_temp_fasta)
    glue("python3 {config$iupred2_path} {iupred_temp_fasta} {iupred_mode}") %>%
      system(intern = TRUE, ignore.stderr = TRUE) %>%
      parse_iupred()
  }, fileext = ".faa", pattern = "sequence_")
}

compute_iupred_long <- partial(compute_iupred, iupred_mode = "long")
compute_iupred_short <- partial(compute_iupred, iupred_mode = "short")
