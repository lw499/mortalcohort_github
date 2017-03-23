library(doParallel)
library(foreach)
getDoParWorkers()  
detectCores()   
cl=makeCluster(detectCores()-2)          
registerDoParallel(cl)       
getDoParWorkers()   

myfunc = function(i)
{
  source("simdatDAG2LMD.R")
  library(geepack)
  source("linearincrements.r")
  set.seed(1127)
  seeds = floor(runif(1000)*10^8);
  
  EXPIT <- function(term) {
    return( exp(term)/(1+exp(term)) )
  }
  
  set.seed(seeds[i])
  tmpdata = gendata(N=500)
  tmpdata$ager = ifelse(tmpdata$shiftage<=4,1,0);
  ###MI
  wd = reshape(tmpdata, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")
  
  ###LI-D uMVN
  wd$alive1=1;wd$alive2=ifelse(wd$s1>=2,1,0);wd$alive3=ifelse(wd$s1>=4,1,0);
  wd$alive4=ifelse(wd$s1>=6,1,0);wd$alive5=ifelse(wd$s1>=8,1,0);
  wd1=wd[wd$s1==0,];wd2 = wd[wd$s1==2,];wd3 = wd[wd$s1==4,];
  wd4 = wd[wd$s1==6,];wd5 = wd[wd$s1==8,];
  wd1=as.matrix(wd1); wd2=as.matrix(wd2); wd3=as.matrix(wd3); wd4=as.matrix(wd4); wd5=as.matrix(wd5);
  tmpli2 = linearincrements(y=wd2[,7:8],x=wd2[,c(3)],alive=wd2[,12:13],method="uMVN",returnimps=T)
  tmpli3 = linearincrements(y=wd3[,7:9],x=wd3[,c(3)],alive=wd3[,12:14],method="uMVN",returnimps=T)
  tmpli4 = linearincrements(y=wd4[,7:10],x=wd4[,c(3)],alive=wd4[,12:15],method="uMVN",returnimps=T)
  tmpli5 = linearincrements(y=wd5[,7:11],x=wd5[,c(3)],alive=wd5[,12:16],method="uMVN",returnimps=T)
  
  tmpliout1 = cbind(wd1[,1:6],wd1[,7:11]); 
  tmpliout2 = cbind(wd2[,1:6],tmpli2$yimp, pek.4=NA, pek.6=NA,pek.8=NA); 
  tmpliout3 = cbind(wd3[,1:6],tmpli3$yimp, pek.6=NA,pek.8=NA); 
  tmpliout4 = cbind(wd4[,1:6],tmpli4$yimp,pek.8=NA); 
  tmpliout5 = cbind(wd5[,1:6],tmpli5$yimp);
  newdatlid = rbind(tmpliout1,tmpliout2,tmpliout3,tmpliout4,tmpliout5)
  
  new = as.data.frame(newdatlid)
  longlid = reshape(new, varying=list(c(7:11)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
  lidatd = longlid[order(longlid$Case),];
  lidatd = lidatd[lidatd$time<= lidatd$s1,]
  fitlid = geeglm(pek ~  factor(time)+ factor(time)*sex, id=Case, data = lidatd)
  PElid=summary(fitlid)$coef[,1];
  
  outtmp = c(PElid)
  outtmp
}

test = foreach(i=1:1000) %dopar% myfunc(i)
test2 = do.call("rbind", test)

write.csv(test2,"li_stratD.csv")


stopCluster(cl)