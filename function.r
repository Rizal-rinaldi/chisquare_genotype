CSQ_genotype <- function (data) {
  as <- table(data)
  p <- (as.numeric((2*as[1])+as[2]))/(2*length(data))
  q <- (as.numeric((2*as[3])+as[2]))/(2*length(data))
  HO <- as.numeric(as[2]/length(data))
  HE <- 1-(p^2+q^2)
  E_Table <- as.vector(c(p^2*length(data), 
                        2*p*q*length(data), 
                        q^2*length(data)))
  O_Table <- as.vector(c((as.numeric(as[1]/length(data)))*(length(data)), 
                        (as.numeric(as[2]/length(data)))*(length(data)), 
                        (as.numeric(as[3]/length(data)))*(length(data))))
  
  Result <- chisq.test(O_Table, p = E_Table/length(data))
frquensi_result <- data.frame(
  Frekuensi_Alel_Dominan = p,
  Frekuensi_Alel_Resesif = q,
  HO = HO,
  HE = HE
)
print(frquensi_result)
result_data <- data.frame(
  genotype = c(row.names(as)),
  Observed = O_Table,
  freq_observed = O_Table/length(data) ,
  Expected = E_Table,
  freq_expected = E_Table/length(data))
print(result_data)
 return(Result)
}
