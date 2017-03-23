library(geepack)
source("linearincrements.r")
source("DGM.R")
set.seed(1012)

EXPIT <- function(term) {
  return( exp(term)/(1+exp(term)) )
}

tmpdata = gendata(N=500)
write.csv(tmpdata, "demodat.csv", row.names=FALSE)
