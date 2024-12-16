CSQ_genotype <- function(data) {

  alleles <- sort(unique(unlist(strsplit(as.character(data), ""))))


  genotypes <- unique(apply(expand.grid(alleles, alleles), 1, function(x) paste(sort(x), collapse = "")))


  as <- table(factor(data, levels = genotypes))


  total_alleles <- 2 * length(data)
  allele_counts <- table(unlist(strsplit(rep(names(as), as), "")))
  p <- allele_counts[alleles[1]] / total_alleles
  q <- allele_counts[alleles[2]] / total_alleles


  HO <- sum(as[grep(paste0("[", alleles[1], "]", "[", alleles[2], "]"), names(as))]) / length(data)
  HE <- 1 - (p^2 + q^2)


  E_Table <- c(p^2 * length(data), 2 * p * q * length(data), q^2 * length(data))

  O_Table <- as.vector(as)


  Result <- chisq.test(O_Table, p = E_Table / length(data))


  frekuensi_result <- data.frame(
    Frekuensi_Alel_Dominan = p,
    Frekuensi_Alel_Resesif = q,
    HO = HO,
    HE = HE
  )
  print(frekuensi_result)

  result_data <- data.frame(
    Genotype = names(as),
    Observed = O_Table,
    Freq_Observed = O_Table / length(data),
    Expected = E_Table,
    Freq_Expected = E_Table / length(data)
  )
  print(result_data)

  return(Result)
}
