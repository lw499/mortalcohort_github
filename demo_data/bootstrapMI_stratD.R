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
  library(geepack); library(norm);
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
  
wd = reshape(dat, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")
wd1=wd[wd$s1==0,];wd2 = wd[wd$s1==2,];wd3 = wd[wd$s1==4,];
wd4 = wd[wd$s1==6,];wd5 = wd[wd$s1==8,];
midatd1 = wd1[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","pek.6","pek.8","s1","iadl1","Smoke","HlthPrv1","mmse1")]
midatd2 = cbind(wd2[, c("Female","Educyrs","shiftage","pek.0","pek.2","iadl1","Smoke","HlthPrv1","mmse1")])
midatd3 = cbind(wd3[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","iadl1","Smoke","HlthPrv1","mmse1")]) 
midatd4 = cbind(wd4[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","pek.6","iadl1","Smoke","HlthPrv1","mmse1")])
midatd5 = cbind(wd5[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","pek.6","pek.8","iadl1","Smoke","HlthPrv1","mmse1")])

imdatad2 = as.matrix(midatd2)
imdatad3 = as.matrix(midatd3); 
imdatad4 = as.matrix(midatd4); 
imdatad5 = as.matrix(midatd5);  

nsim = 30
PE = VAR = matrix(0,nsim,11)
rngseed(112789)
for(i in 1:nsim)
{
  s = prelim.norm(imdatad2)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdatad2)
  tmp = as.data.frame(ximp); tmp$pek.4=tmp$pek.6=tmp$pek.8=NA; tmp$s1=wd2$s1;
  
  mydat = rbind(midatd1, tmp)
  
  s = prelim.norm(imdatad3)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdatad3)
  tmp = as.data.frame(ximp); tmp$pek.6=tmp$pek.8=NA; tmp$s1=wd3$s1;
  
  mydat = rbind(mydat, tmp)
  
  s = prelim.norm(imdatad4)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdatad4)
  tmp = as.data.frame(ximp); tmp$pek.8=NA; tmp$s1=wd4$s1;
  
  mydat = rbind(mydat, tmp)
  
  s = prelim.norm(imdatad5)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdatad5)
  tmp = as.data.frame(ximp); tmp$s1=wd5$s1;
  
  mydat = rbind(mydat, tmp)
  
  long = reshape(mydat, varying=list(c(4:8)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
  longsort = long[order(long$id),];
  midatd = longsort[longsort$time<= longsort$s1,]
  
  fitmi = geeglm(pek ~  time +  I(time^2) + Female + shiftage + Educyrs + Smoke + time*Female + time*shiftage + time*Educyrs + time*Smoke, data = midatd, id=id,na.action=na.omit)
  PE[i,]=(fitmi)$coef
  VAR[i,] = summary(fitmi)$coef[,2]^2
  
}
  COEFd = apply(PE, 2, mean)


  return(c(COEFd))
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"listMId.csv")



