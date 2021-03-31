phaseMyAllelesBasedOnHaplotypes <- function(m, th = 0) {
  require(data.table)
  require(pryr)
  a <- copy(m$A)
  b <- copy(m$B)
  c <- copy(m$A)
  
  wi <- which(m$features$allele < th)
  a[wi, ] <- b[wi, ]
  b[wi, ] <- c[wi, ]
  a <- a[-which(m$features$allele == 0)]
  b <- b[-which(m$features$allele == 0)]
  features <- m$features[-which(m$features$allele == 0)]
  return(list(A = a, B = b, features = features))
}