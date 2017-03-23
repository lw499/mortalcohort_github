setwd("~/Dropbox/phd/Missing data comparisons/pek/simwithoutedu_tonly/mortalcohort_github/demo_data")

tmpdata = read.csv("demodat.csv", sep=",")
fixdatw = reshape(tmpdata, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")

library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(detectCores()-2)                                              
registerDoParallel(cl)                                                                          
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  library(geepack);library(norm);
  set.seed(30)
  seeds = floor(runif(1000)*10^8);
  
  EXPIT <- function(term) {
    return( exp(term)/(1+exp(term)) )}
  ###Bootstrap
  set.seed(seeds[m])
  random = fixdatw[sample(1:nrow(fixdatw), 500, replace=T),]
  longdat = reshape(random, varying=list(c(10:14)), v.names = c("pek"), new.row.names=1:100000, times=c(0,2,4,6,8), direction="long")
  dat = longdat[order(longdat$id),]; id=dat$id; dat$id=NULL;
  shiftage = dat$shiftage; s1=dat$s1; pek=dat$pek; dat$Case=id;
  dat$shiftage=NULL; dat$s1=NULL; dat$pek=NULL;
  dat$pek=pek; dat$shiftage=shiftage; dat$s1=s1;
  
w = reshape(dat, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")
midat = w[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","pek.6","pek.8","iadl1","Smoke","HlthPrv1","mmse1")]

imdata = as.matrix(midat)

nsim = 30
PE = VAR = matrix(0,nsim,11)
rngseed(112789)
for(i in 1:nsim)
{
  s = prelim.norm(imdata)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdata)
  tmp = as.data.frame(ximp); tmp$s1=w$s1;
  
  long = reshape(tmp, varying=list(c(4:8)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
  longsort = long[order(long$id),];
  midat = longsort[longsort$time<= longsort$s1,]
  
  fitmi = geeglm(pek ~  time +  I(time^2) + Female + shiftage + Educyrs + Smoke + time*Female + time*shiftage + time*Educyrs + time*Smoke, data = midat, id=id,na.action=na.omit)
  PE[i,]=(fitmi)$coef
  VAR[i,] = summary(fitmi)$coef[,2]^2
}

COEF = apply(PE, 2, mean)

  return(c(COEF))
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"listMI.csv")



