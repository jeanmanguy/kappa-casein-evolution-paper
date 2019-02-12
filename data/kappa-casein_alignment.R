# load alignment
cache("kappa_casein_alignment", {
  Biostrings::readAAMultipleAlignment(file.path("data", "kappa-casein_alignment.fa"), format = "fasta")
})
