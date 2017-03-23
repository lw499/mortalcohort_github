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
  library(geepack)
  source("linearincrements.r")
  source("simdatDAG1LMD.R")
  set.seed(60)
  seeds = floor(runif(1000)*10^8);
  
  EXPIT <- function(term) {
    return( exp(term)/(1+exp(term)) )
  }
  
  set.seed(seeds[m])
  tmpdata = gendata(N=500)
  tmpdata$ager = ifelse(tmpdata$shiftage<=4,1,0);
  
  ###IPWIEE
  dat = tmpdata
  dat = dat[dat$s1>=dat$time,]
  dat$yl = c(NA, dat$pek[1:nrow(dat)-1])
  dat$yl = ifelse(dat$time==0, NA, dat$yl)
  dat$lyl = c(NA, dat$yl[1:nrow(dat)-1])
  dat$lyl = ifelse(dat$time<=2, NA, dat$lyl)
  dat$llyl = c(NA, dat$lyl[1:nrow(dat)-1])
  dat$llyl = ifelse(dat$time<=4, NA, dat$llyl)
  dat$lllyl = c(NA, dat$llyl[1:nrow(dat)-1])
  dat$lllyl = ifelse(dat$time<=6, NA, dat$lllyl)
  dat$Ri = ifelse(!is.na(dat$pek),1,0)
  dat$Rim1 = c(NA, dat$Ri[1:nrow(dat)-1])   ##Status of the last observation
  dat$Rim1 = ifelse(dat$time==0, NA, dat$Rim1)
  ###IPW-IEE-D
  ipwdat = dat[,c("Case", "time", "shiftage", "sex", "edu", "Ri","pek","s1","Rim1","yl", "lyl","llyl","lllyl","ager")]
  ipwdat$pred_obs=NA;ipwdat$Rim1=ifelse(ipwdat$time==0, 0, ipwdat$Rim1);
  
  ipwdat22 = ipwdat[ipwdat$s1==2 & ipwdat$Rim1==1,]; 
  ipwdat2 = ipwdat22[,c("yl", "shiftage", "sex","edu","time","s1")]
  logistfit = glm(Ri ~ yl +sex, family = binomial(link=logit),data = ipwdat22)
  ipwdat22$pred_obs = predict(logistfit, newdata = ipwdat2, type="response")
  
  ##########
  
  ipwdat44 = ipwdat[ipwdat$s1==4 & ipwdat$Rim1==1,];
  ipwdat442 = ipwdat44[ipwdat44$time==2,]; ipwdat444 = ipwdat44[ipwdat44$time==4,]; 
  ipwdat42 = ipwdat442[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  ipwdat44 = ipwdat444[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  logistfit2 = glm(Ri ~ yl +sex, control=list(maxit=50),family = binomial(link=logit),data = ipwdat442)
  logistfit4 = glm(Ri ~ yl +sex,control=list(maxit=50), family = binomial(link=logit),data = ipwdat444)
  ipwdat442$pred_obs = predict(logistfit2, newdata = ipwdat42, type="response")
  ipwdat444$pred_obs = predict(logistfit4, newdata = ipwdat44, type="response")
  
  ###############
  
  ipwdat66 = ipwdat[ipwdat$s1==6 & ipwdat$Rim1==1,];
  ipwdat662 = ipwdat66[ipwdat66$time==2,]; ipwdat664 = ipwdat66[ipwdat66$time==4,]; ipwdat666 = ipwdat66[ipwdat66$time==6,];
  ipwdat62 = ipwdat662[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  ipwdat64 = ipwdat664[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  ipwdat66 = ipwdat666[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  logistfit2 = glm(Ri ~ yl+sex,control=list(maxit=50), family = binomial(link=logit),data = ipwdat662)
  logistfit4 = glm(Ri ~ yl+sex,control=list(maxit=50), family = binomial(link=logit),data = ipwdat664)
  logistfit6 = glm(Ri ~ yl+sex,control=list(maxit=50), family = binomial(link=logit),data = ipwdat666)
  ipwdat662$pred_obs = predict(logistfit2, newdata = ipwdat62, type="response")
  ipwdat664$pred_obs = predict(logistfit4, newdata = ipwdat64, type="response")
  ipwdat666$pred_obs = predict(logistfit6, newdata = ipwdat66, type="response")
  
  #################
  
  ipwdat88 = ipwdat[ipwdat$s1==8 & ipwdat$Rim1==1,]
  ipwdat882 = ipwdat88[ipwdat88$time==2,]; ipwdat884 = ipwdat88[ipwdat88$time==4,]; ipwdat886 = ipwdat88[ipwdat88$time==6,]; ipwdat888 = ipwdat88[ipwdat88$time==8,];
  ipwdat82 = ipwdat882[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  ipwdat84 = ipwdat884[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  ipwdat86 = ipwdat886[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  ipwdat88 = ipwdat888[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  logistfit2 = glm(Ri ~ yl +sex, control=list(maxit=50),family = binomial(link=logit),data = ipwdat882)
  logistfit4 = glm(Ri ~ yl +sex,control=list(maxit=50), family = binomial(link=logit),data = ipwdat884)
  logistfit6 = glm(Ri ~ yl +sex,control=list(maxit=50), family = binomial(link=logit),data = ipwdat886)
  logistfit8 = glm(Ri ~ yl +sex,control=list(maxit=50), family = binomial(link=logit),data = ipwdat888)
  ipwdat882$pred_obs = predict(logistfit2, newdata = ipwdat82, type="response")
  ipwdat884$pred_obs = predict(logistfit4, newdata = ipwdat84, type="response")
  ipwdat886$pred_obs = predict(logistfit6, newdata = ipwdat86, type="response")
  ipwdat888$pred_obs = predict(logistfit8, newdata = ipwdat88, type="response")
  ##PREDICT for waves
  
  combine = rbind(ipwdat22,ipwdat442,ipwdat444,ipwdat662,ipwdat664,ipwdat666,ipwdat882,ipwdat884,ipwdat886,ipwdat888)
  combine = combine[order(combine$Case),]; combine$pred_obs=ifelse(is.na(combine$pek),NA,combine$pred_obs);
  combine$Rim1=combine$Ri=combine$yl=combine$lyl=combine$llyl=combine$lllyl=NULL;
  combw = reshape(combine, idvar = "Case", v.names=c("pek","pred_obs"),timevar="time", direction="wide")
  combw$pi.2=1/combw$pred_obs.2; combw$pi.4=1/(combw$pred_obs.2*combw$pred_obs.4);
  combw$pi.6=1/(combw$pred_obs.2*combw$pred_obs.4*combw$pred_obs.6); combw$pi.8=1/(combw$pred_obs.2*combw$pred_obs.4*combw$pred_obs.6*combw$pred_obs.8);
  combw$pred_obs.2 = combw$pred_obs.4 = combw$pred_obs.6 = combw$pred_obs.8 = NULL;
  
  long = reshape(combw, varying=list(c(7:10),c(11:14)), v.names = c("pek","pi"), new.row.names=1:1000000, times=c(2,4,6,8), direction="long")
  long$id=NULL;
  ipwdat0 = ipwdat[ipwdat$time==0,]; ipwdat0$pi=1; 
  ipwdat0$yl=ipwdat0$lyl=ipwdat0$llyl=ipwdat0$lllyl=ipwdat0$pred_obs=ipwdat0$Ri=ipwdat0$Rim1=NULL;
  restipwdat = ipwdat[ipwdat$Rim1==0 & ipwdat$time>0,]; restipwdat$pi=NA;
  restipwdat$yl=restipwdat$lyl=restipwdat$llyl=restipwdat$lllyl=restipwdat$pred_obs=restipwdat$Ri=restipwdat$Rim1=NULL;
  ripwdatd = rbind(ipwdat0,long,restipwdat)
  
  ripwdatd = ripwdatd[order(ripwdatd$Case),];
  ripwdatd = ripwdatd[ripwdatd$time<= ripwdatd$s1,]
  
  ipw.fitd = geeglm(pek ~ factor(time) + factor(time)*sex, data = ripwdatd, weights = pi, id=Case)
  
  outtmp = c(
    summary(ipw.fitd)$coefficient[,1])
  outtmp
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"ipw_stratD.csv")

