# compute net charges for each each species for each part and mature sequence
cache("kappa_caseins_ncpr", {
  kappa_caseins_parts_mature %>%
    mutate(pH = list(seq(2.5, 7.5, 0.25))) %>%
    unnest(pH) %>%
    rowwise() %>%
    mutate(
      net_charge = Peptides::charge(subseq, pH = pH, pKscale = config$charge_scale),
      length = nchar(subseq)
      ) %>%
    ungroup() %>%
    select(species, part, pH, net_charge, length) %>%
    mutate(net_charge_per_residue = net_charge / length) %>%
    select(-length) %>%
    order_cleavage_parts()
}, depends = c("kappa_caseins_parts_mature"))
