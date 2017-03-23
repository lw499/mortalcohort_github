##IN FILE, USING i AS THE RNGSEED

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
  source("simdatDAG1LMD.R")
  set.seed(60)
  seeds = floor(runif(1000)*10^8);
  
  EXPIT <- function(term) {
    return( exp(term)/(1+exp(term)) )
  }
  
  set.seed(seeds[m])
  tmpdata = gendata(N=500)
  tmpdata$ager = ifelse(tmpdata$shiftage<=4,1,0);
  ###MI
  wd = reshape(tmpdata, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")
  wd1=wd[wd$s1==0,];wd2 = wd[wd$s1==2,];wd3 = wd[wd$s1==4,];
  wd4 = wd[wd$s1==6,];wd5 = wd[wd$s1==8,];
  midatd1 = wd1[, c("sex","pek.0","pek.2","pek.4","pek.6","pek.8","s1")]
  midatd2 = wd2[, c("sex","pek.0","pek.2")]
  midatd3 = wd3[, c("sex","pek.0","pek.2","pek.4")]
  midatd4 = wd4[, c("sex","pek.0","pek.2","pek.4","pek.6")]
  midatd5 = wd5[, c("sex","pek.0","pek.2","pek.4","pek.6","pek.8")]
  
  imdatad2 = as.matrix(midatd2)
  imdatad3 = as.matrix(midatd3)
  imdatad4 = as.matrix(midatd4)
  imdatad5 = as.matrix(midatd5)
  
  nsim = 30
  rngseed(12345)
  
  PE = matrix(0,nsim,10); mydat=NULL
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
    
    long = reshape(mydat, varying=list(c(2:6)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
    longsort = long[order(long$id),];
    midatd = longsort[longsort$time<= longsort$s1,]
    
    #midatd$ager = ifelse(midatd$shiftage<=4,1,0); midatd$ager = as.factor(midatd$ager);
    fitmi = geeglm(pek ~  factor(time) + factor(time)*sex, data = midatd, id=id,na.action=na.omit)
    PE[i,]=summary(fitmi)$coefficient[,1]
  }
  
  COEFd = apply(PE, 2, mean)
  
  
  outtmp = c(COEFd)
  outtmp
  #cat(m, "\n")
}



test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"mi_stratD.csv")


stopCluster(cl)



