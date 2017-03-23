library(doParallel)
library(foreach)
getDoParWorkers()  
detectCores()   
cl=makeCluster(detectCores()-2)          
registerDoParallel(cl)       
getDoParWorkers()   

myfunc = function(i)
{
  source("simdatDAG3LMD.R")
  library(geepack)
  
  set.seed(1127)
  seeds = floor(runif(1000)*10^8);
  
  EXPIT <- function(term) {
    return( exp(term)/(1+exp(term)) )
  }
  
  set.seed(seeds[i])
  tmpdata = gendata(N=500)
  tmpdata$ager = ifelse(tmpdata$shiftage<=4,1,0);  
  
  ###IPWIEE
  ###IPWIEE data
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
  logistfit2 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat442)
  logistfit4 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat444)
  ipwdat442$pred_obs = predict(logistfit2, newdata = ipwdat42, type="response")
  ipwdat444$pred_obs = predict(logistfit4, newdata = ipwdat44, type="response")
  
  ###############
  
  ipwdat66 = ipwdat[ipwdat$s1==6 & ipwdat$Rim1==1,];
  ipwdat662 = ipwdat66[ipwdat66$time==2,]; ipwdat664 = ipwdat66[ipwdat66$time==4,]; ipwdat666 = ipwdat66[ipwdat66$time==6,];
  ipwdat62 = ipwdat662[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  ipwdat64 = ipwdat664[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  ipwdat66 = ipwdat666[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  logistfit2 = glm(Ri ~ yl+sex, family = binomial(link=logit),data = ipwdat662)
  logistfit4 = glm(Ri ~ yl+sex, family = binomial(link=logit),data = ipwdat664)
  logistfit6 = glm(Ri ~ yl+sex, family = binomial(link=logit),data = ipwdat666)
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
  logistfit2 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat882)
  logistfit4 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat884)
  logistfit6 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat886)
  logistfit8 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat888)
  ipwdat882$pred_obs = predict(logistfit2, newdata = ipwdat82, type="response")
  ipwdat884$pred_obs = predict(logistfit4, newdata = ipwdat84, type="response")
  ipwdat886$pred_obs = predict(logistfit6, newdata = ipwdat86, type="response")
  ipwdat888$pred_obs = predict(logistfit8, newdata = ipwdat88, type="response")
  ##PREDICT for waves
  
  combine = rbind(ipwdat22,ipwdat442,ipwdat444,ipwdat662,ipwdat664,ipwdat666,ipwdat882,ipwdat884,ipwdat886,ipwdat888)
  combine = combine[order(combine$Case),]; 
  combine$Rim1=combine$Ri=combine$yl=combine$lyl=combine$llyl=combine$lllyl=NULL;
  combw = reshape(combine, idvar = "Case", v.names=c("pek","pred_obs"),timevar="time", direction="wide")
  combw$pi.2=1/combw$pred_obs.2; combw$pi.4=1/(combw$pred_obs.2*combw$pred_obs.4);
  combw$pi.6=1/(combw$pred_obs.2*combw$pred_obs.4*combw$pred_obs.6); combw$pi.8=1/(combw$pred_obs.2*combw$pred_obs.4*combw$pred_obs.6*combw$pred_obs.8);
  combw$pred_obs.2 = combw$pred_obs.4 = combw$pred_obs.6 = combw$pred_obs.8 = NULL;
  
  long = reshape(combw, varying=list(c(7:10),c(11:14)), v.names = c("pek","pi"), new.row.names=1:1000000, times=c(2,4,6,8), direction="long")
  long$id=NULL;
  ipwdat0 = ipwdat[ipwdat$time==0,]; ipwdat0$pi=1; 
  ipwdat0$yl=ipwdat0$lyl=ipwdat0$llyl=ipwdat0$lllyl=ipwdat0$pred_obs=ipwdat0$Ri=ipwdat0$Rim1=NULL;
  ripwdatd = rbind(ipwdat0,long)
  
  ripwdatd = ripwdatd[order(ripwdatd$Case),];
  ripwdatd = ripwdatd[ripwdatd$time<= ripwdatd$s1,] 
  ripwdatd$Ri = ifelse(!is.na(ripwdatd$pek),1,0)
  
  ## Augmented IPW-D ###
  #### Imputation models (10)
  
  midat12 = dat[dat$time==2 & dat$Rim1==1 & dat$s1==2,]
  mifit1 = lm(pek  ~ yl + sex, data = midat12)
  mi12 = midat12[,c("yl", "shiftage", "sex","edu","time","s1")]
  midat12$pred_obs = predict(mifit1, newdata = mi12, type="response")
  midat12$pred_obs = ifelse(!is.na(midat12$pek), midat12$pek, midat12$pred_obs)
  midat12$pek = midat12$pred_obs; 
  
    midat122 = dat[dat$time==2 & dat$Rim1==1 & dat$s1==4,]
    mifit12 = lm(pek  ~ yl + sex, data = midat122)
    mi12 = midat122[,c("yl", "shiftage", "sex","edu","time","s1")]
    midat122$pred_obs = predict(mifit12, newdata = mi12, type="response")
    midat122$pred_obs = ifelse(!is.na(midat122$pek), midat122$pek, midat122$pred_obs)
    midat122$pek = midat122$pred_obs; 
  
    midat123 = dat[dat$time==2 & dat$Rim1==1 & dat$s1==6,]
    mifit13 = lm(pek  ~ yl + sex, data = midat123)
    mi12 = midat123[,c("yl", "shiftage", "sex","edu","time","s1")]
    midat123$pred_obs = predict(mifit13, newdata = mi12, type="response")
    midat123$pred_obs = ifelse(!is.na(midat123$pek), midat123$pek, midat123$pred_obs)
    midat123$pek = midat123$pred_obs; 
  
    midat124 = dat[dat$time==2 & dat$Rim1==1 & dat$s1==8,]
    mifit14 = lm(pek  ~ yl + sex, data = midat124)
    mi12 = midat124[,c("yl", "shiftage", "sex","edu","time","s1")]
    midat124$pred_obs = predict(mifit14, newdata = mi12, type="response")
    midat124$pred_obs = ifelse(!is.na(midat124$pek), midat124$pek, midat124$pred_obs)
    midat124$pek = midat124$pred_obs; 
  
  midat23 = dat[dat$time==4 & dat$Rim1==1 & dat$s1==4,]
  mifit2 = lm(pek  ~ yl + sex, data = midat23)
  mi23 = midat23[,c("yl", "lyl","shiftage", "sex","edu","time","s1")]
  midat23$pred_obs = predict(mifit2, newdata = mi23, type="response")
  midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)
  
    midat232 = dat[dat$time==4 & dat$Rim1==1 & dat$s1==6,]
    mifit22 = lm(pek  ~ yl + sex, data = midat232)
    mi23 = midat232[,c("yl", "lyl","shiftage", "sex","edu","time","s1")]
    midat232$pred_obs = predict(mifit22, newdata = mi23, type="response")
    midat232$pred_obs = ifelse(!is.na(midat232$pek), midat232$pek, midat232$pred_obs)
  
    midat233 = dat[dat$time==4 & dat$Rim1==1 & dat$s1==8,]
    mifit23 = lm(pek  ~ yl + sex, data = midat233)
    mi23 = midat233[,c("yl", "lyl","shiftage", "sex","edu","time","s1")]
    midat233$pred_obs = predict(mifit23, newdata = mi23, type="response")
    midat233$pred_obs = ifelse(!is.na(midat233$pek), midat233$pek, midat233$pred_obs)
  
  midat34 = dat[dat$time==6 & dat$Rim1==1 & dat$s1==6,]
  mifit3 = lm(pek  ~ yl + sex, data = midat34)
  mi34 = midat34[,c("yl", "lyl","llyl", "shiftage", "sex","edu","time","s1")]
  midat34$pred_obs = predict(mifit3, newdata = mi34, type="response")
  midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)
  
    midat342 = dat[dat$time==6 & dat$Rim1==1 & dat$s1==8,]
    mifit32 = lm(pek  ~ yl + sex, data = midat342)
    mi34 = midat342[,c("yl", "lyl","llyl", "shiftage", "sex","edu","time","s1")]
    midat342$pred_obs = predict(mifit32, newdata = mi34, type="response")
    midat342$pred_obs = ifelse(!is.na(midat342$pek), midat342$pek, midat342$pred_obs)
  
  midat45 = dat[dat$time==8 & dat$Rim1==1,]
  mifit4 = lm(pek  ~ yl + sex, data = midat45)
  mi45 = midat45[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  midat45$pred_obs = predict(mifit4, newdata = mi45, type="response")
  midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)
  
  ###########################Fitted values are not covariates but outcomes
  tmpmidat23 = dat[dat$time==4 & is.na(dat$yl) & dat$s1==4,]
  midat23$pek = midat23$pred_obs; midat23$pred_obs=NULL; midat23 = rbind(tmpmidat23,midat23);
  midat23$Ri = ifelse(!is.na(midat23$pek),1,0)
  mifit5 = lm(pek  ~  lyl + sex, data = midat23)
  mi23 = midat23[,c("lyl", "shiftage", "sex","edu","time","s1")]
  midat23$pred_obs = predict(mifit5, newdata = mi23, type="response")
  midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)
  midat23$pek = midat23$pred_obs; 
    
    tmpmidat23 = dat[dat$time==4 & is.na(dat$yl) & dat$s1==6,]
    midat232$pek = midat232$pred_obs; midat232$pred_obs=NULL; midat232 = rbind(tmpmidat23,midat232);
    midat232$Ri = ifelse(!is.na(midat232$pek),1,0)
    mifit52 = lm(pek  ~  lyl + sex, data = midat232)
    mi23 = midat232[,c("lyl", "shiftage", "sex","edu","time","s1")]
    midat232$pred_obs = predict(mifit52, newdata = mi23, type="response")
    midat232$pred_obs = ifelse(!is.na(midat232$pek), midat232$pek, midat232$pred_obs)
    midat232$pek = midat232$pred_obs;
  
    tmpmidat23 = dat[dat$time==4 & is.na(dat$yl) & dat$s1==8,]
    midat233$pek = midat233$pred_obs; midat233$pred_obs=NULL; midat233 = rbind(tmpmidat23,midat233);
    midat233$Ri = ifelse(!is.na(midat233$pek),1,0)
    mifit53 = lm(pek  ~  lyl + sex, data = midat233)
    mi23 = midat233[,c("lyl", "shiftage", "sex","edu","time","s1")]
    midat233$pred_obs = predict(mifit53, newdata = mi23, type="response")
    midat233$pred_obs = ifelse(!is.na(midat233$pek), midat233$pek, midat233$pred_obs)
    midat233$pek = midat233$pred_obs; 
  
  tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & !is.na(dat$lyl) & dat$s1==6,]
  midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
  midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
  mifit6 = lm(pek  ~ lyl + sex, data = midat34)
  mi34 = midat34[,c("lyl","llyl", "shiftage", "sex","edu","time","s1")]
  midat34$pred_obs = predict(mifit6, newdata = mi34, type="response")
  midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)
  
    tmpmidat342 = dat[dat$time==6 & is.na(dat$yl) & !is.na(dat$lyl) & dat$s1==8,]
    midat342$pek = midat342$pred_obs; midat342$pred_obs=NULL; midat342 = rbind(tmpmidat342,midat342);
    midat342$Ri = ifelse(!is.na(midat342$pek),1,0)
    mifit62 = lm(pek  ~ lyl + sex, data = midat342)
    mi34 = midat342[,c("lyl","llyl", "shiftage", "sex","edu","time","s1")]
    midat342$pred_obs = predict(mifit62, newdata = mi34, type="response")
    midat342$pred_obs = ifelse(!is.na(midat342$pek), midat342$pek, midat342$pred_obs)
  
  tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & !is.na(dat$lyl),]
  midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
  midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
  mifit7 = lm(pek  ~ lyl + sex, data = midat45)
  mi45 = midat45[,c("lyl","llyl","lllyl" , "shiftage", "sex","edu","time","s1")]
  midat45$pred_obs = predict(mifit7, newdata = mi45, type="response")
  midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)
  
  ######
  tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & is.na(dat$lyl)  & !is.na(dat$llyl) & dat$s1==6,]
  midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
  midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
  mifit8 = lm(pek  ~ llyl + sex, data = midat34)
  mi34 = midat34[,c("llyl", "shiftage", "sex","edu","time","s1")]
  midat34$pred_obs = predict(mifit8, newdata = mi34, type="response")
  midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)
  midat34$pek = midat34$pred_obs; 
    
    tmpmidat342 = dat[dat$time==6 & is.na(dat$yl) & is.na(dat$lyl)  & !is.na(dat$llyl) & dat$s1==8,]
    midat342$pek = midat342$pred_obs; midat342$pred_obs=NULL; midat342 = rbind(tmpmidat342,midat342);
    midat342$Ri = ifelse(!is.na(midat342$pek),1,0)
    mifit82 = lm(pek  ~ llyl + sex, data = midat342)
    mi34 = midat342[,c("llyl", "shiftage", "sex","edu","time","s1")]
    midat342$pred_obs = predict(mifit82, newdata = mi34, type="response")
    midat342$pred_obs = ifelse(!is.na(midat342$pek), midat342$pek, midat342$pred_obs)
    midat342$pek = midat342$pred_obs; 

  tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$lyl) & !is.na(dat$llyl),]
  midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
  midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
  mifit9 = lm(pek  ~  llyl + sex, data = midat45)
  mi45 = midat45[,c("llyl","lllyl", "shiftage", "sex","edu","time","s1")]
  midat45$pred_obs = predict(mifit9, newdata = mi45, type="response")
  midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)
  
  ######
  tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$llyl),]
  midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
  midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
  mifit10 = lm(pek  ~ lllyl + sex, data = midat45)
  mi45 = midat45[,c("lllyl", "shiftage", "sex","edu","time","s1")]
  midat45$pred_obs = predict(mifit10, newdata = mi45, type="response")
  midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)
  midat45$pek = midat45$pred_obs; 
  
  
  #####
  midat = rbind(midat12, midat122, midat123, midat124, midat23, midat232, midat233, midat34, midat342, midat45); 
  midat$pred_obs = NULL;
  midat = rbind(dat[dat$time==0,],midat)
  midatsortd = midat[order(midat$Case),]
  midatsortd = midatsortd[midatsortd$time<=midatsortd$s1,]
  
  ########## complete data set
  midatsortd$pi=1/ripwdatd$pi
  midatsortd$Rim1 = midatsortd$yl = midatsortd$lyl = midatsortd$llyl = midatsortd$lllyl=NULL;
  midatsortd$Ri = ripwdatd$Ri; midatsortd$digreal = ripwdatd$pek;
  ##NR##
  
  beta = c(18,0,0,0,0);
  diff=10;r=1;crit=0.00001
  while(diff>=crit)
  {
    phi = phic = 1;
    U = 0; dU=0;
    auniqueid = unique(midatsortd$Case);
    for(i in auniqueid)
    {
      tmp5 = midatsortd[midatsortd$Case==i,];       
      n = nrow(tmp5);
      tmpdiag = diag(5); tmpdiag = tmpdiag[1:n,!tmpdiag[,1]];
      x = matrix(c(rep(1, n),as.vector(tmpdiag)),nrow=n) #nix4
      last = tmp5[n,]; last$time=10; last2 = do.call("rbind",replicate(5-n, last, simplify = FALSE));tmp5 = rbind(tmp5,last2);
      
      ##j=1
      xtmp = as.data.frame(matrix(c(rep(tmp5$digreal[1],4), tmp5$sex[1], tmp5$edu[1]),nrow=1))
      colnames(xtmp) = c("yl","lyl","llyl","lllyl","sex","edu")    
         
      W1 = NULL;
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==0) {W1=c(tmp5$digreal[1],predict(mifit1, newdata = xtmp, type="response"), predict(mifit5, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
      
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==2) {W1=c(tmp5$digreal[1],predict(mifit1, newdata = xtmp, type="response"), predict(mifit5, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
      
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==4) {W1=c(tmp5$digreal[1],predict(mifit12, newdata = xtmp, type="response"), predict(mifit5, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
      
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==6) {W1=c(tmp5$digreal[1],predict(mifit13, newdata = xtmp, type="response"), predict(mifit52, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
 
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==8) {W1=c(tmp5$digreal[1],predict(mifit14, newdata = xtmp, type="response"), predict(mifit53, newdata = xtmp, type="response"),predict(mifit82, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
            
      #j=2
      xtmp2l = data.frame(cbind(tmp5$digreal[2],tmp5$digreal[1],tmp5$sex[1],tmp5$edu[1],tmp5$s1[1]));
      colnames(xtmp2l) = c("yl","lyl","sex","edu","s1")
      xtmp2 = data.frame(cbind(tmp5$digreal[2],tmp5$digreal[1],tmp5$sex[1],tmp5$edu[1],tmp5$s1[1]));
      colnames(xtmp2) = c("lyl","llyl","sex","edu","s1")
      xtmp2ll = data.frame(cbind(tmp5$digreal[2],tmp5$digreal[1],tmp5$sex[1],tmp5$edu[1]));
      colnames(xtmp2ll) = c("llyl","lllyl","sex","edu")
      
      W2=NULL
      if(!is.na(tmp5$digreal[2]) & tmp5$s1[1]==2){W2=c(tmp5$digreal[1], tmp5$digreal[2],predict(mifit2, newdata=xtmp2l,type="response"),predict(mifit6, newdata=xtmp2,type="response"),predict(mifit9, newdata=xtmp2ll,type="response"))[1:n]}else{W2=W2}
      
      if(!is.na(tmp5$digreal[2]) & tmp5$s1[1]==4){W2=c(tmp5$digreal[1], tmp5$digreal[2],predict(mifit2, newdata=xtmp2l,type="response"),predict(mifit6, newdata=xtmp2,type="response"),predict(mifit9, newdata=xtmp2ll,type="response"))[1:n]}else{W2=W2}
      
      if(!is.na(tmp5$digreal[2]) & tmp5$s1[1]==6){W2=c(tmp5$digreal[1], tmp5$digreal[2],predict(mifit22, newdata=xtmp2l,type="response"),predict(mifit6, newdata=xtmp2,type="response"),predict(mifit9, newdata=xtmp2ll,type="response"))[1:n]}else{W2=W2}
      
      if(!is.na(tmp5$digreal[2]) & tmp5$s1[1]==8){W2=c(tmp5$digreal[1], tmp5$digreal[2],predict(mifit23, newdata=xtmp2l,type="response"),predict(mifit62, newdata=xtmp2,type="response"),predict(mifit9, newdata=xtmp2ll,type="response"))[1:n]}else{W2=W2}
      
      #j=3
      xtmp3 = data.frame(cbind(tmp5$digreal[3],tmp5$digreal[2],tmp5$digreal[1],tmp5$sex[1],tmp5$edu[1],tmp5$s1[1]));
      colnames(xtmp3) = c("yl","lyl","llyl","sex","edu","s1")
      xtmp3l = data.frame(cbind(tmp5$digreal[3],tmp5$digreal[2],tmp5$digreal[1],tmp5$sex[1],tmp5$edu[1]));
      colnames(xtmp3l) = c("lyl","llyl","lllyl","sex","edu")
      
      W3=NULL
      if(!is.na(tmp5$digreal[3]) & tmp5$s1[1]==6){W3=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],predict(mifit3, newdata=xtmp3,type="response"),predict(mifit7, newdata=xtmp3l,type="response"))[1:n]}else{W3=W3}
      
      if(!is.na(tmp5$digreal[3]) & tmp5$s1[1]==8){W3=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],predict(mifit32, newdata=xtmp3,type="response"),predict(mifit7, newdata=xtmp3l,type="response"))[1:n]}else{W3=W3}
      
      #j=4
      xtmp4 = data.frame(cbind(tmp5$digreal[4],tmp5$digreal[3],tmp5$digreal[2],tmp5$digreal[1],tmp5$sex[1],tmp5$edu[1]));
      colnames(xtmp4) = c("yl","lyl","llyl","lllyl","sex","edu")
      
      W4=NULL
      if(!is.na(tmp5$digreal[3]) & tmp5$s1[1]==8){W4=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],tmp5$digreal[4],predict(mifit4, newdata=xtmp4,type="response"))[1:n]}else{W4=W4}
      
      fitted = x %*% t(t(beta)); tmp5 = tmp5[tmp5$s1>=tmp5$time,];
      
      ### For the left hand side of AIPW (complete side)
      var = diag(phi,n); varc = diag(phic,n)
      if(!is.na(tmp5$pi[n]) & !is.na(tmp5$digreal[n])) {uleft=(1/tmp5$pi[n])*(t(x) %*% var %*% (tmp5$digreal - x%*%t(t(beta)) ) )} else{uleft=matrix(c(rep(0,ncol(x))),ncol=1)}
      if(!is.na(tmp5$pi[n]) & !is.na(tmp5$digreal[n])) {duleft = (1/tmp5$pi[n])*(t(x) %*% varc %*% x)} else{duleft=diag(0,ncol(x))}
      
      ### For the right hand side
      if(!is.null(W1) & !is.na(tmp5$pi[2])) {c1=as.matrix((tmp5$Ri[1]/(tmp5$pi[1])-(tmp5$Ri[2])/(tmp5$pi[2]))*(W1-fitted)); p1=(tmp5$Ri[1]/(tmp5$pi[1])-(tmp5$Ri[2])/(tmp5$pi[2]))} else{c1=matrix(c(rep(0,n)),ncol=1);p1=0}
      if(!is.null(W2)& !is.na(tmp5$pi[3])) {c2=as.matrix((tmp5$Ri[2]/(tmp5$pi[2])-(tmp5$Ri[3])/(tmp5$pi[3]))*(W2-fitted));p2=(tmp5$Ri[2]/(tmp5$pi[2])-(tmp5$Ri[3])/(tmp5$pi[3]))} else{c2=matrix(c(rep(0,n)),ncol=1); p2=0;}
      if(!is.null(W3)& !is.na(tmp5$pi[4])) {c3=as.matrix((tmp5$Ri[3]/(tmp5$pi[3])-(tmp5$Ri[4])/(tmp5$pi[4]))*(W3-fitted)); p3=(tmp5$Ri[3]/(tmp5$pi[3])-(tmp5$Ri[4])/(tmp5$pi[4]))} else{c3=matrix(c(rep(0,n)),ncol=1); p3=0;}
      if(!is.null(W4)& !is.na(tmp5$pi[5])) {c4=as.matrix((tmp5$Ri[4]/(tmp5$pi[4])-(tmp5$Ri[5])/(tmp5$pi[5]))*(W4-fitted));p4=(tmp5$Ri[4]/(tmp5$pi[4])-(tmp5$Ri[5])/(tmp5$pi[5]))} else{c4=matrix(c(rep(0,n)),ncol=1);p4=0;}
      uright = t(x) %*% var %*% (c1+c2+c3+c4)
      duright = t(x) %*% var %*% x *(p1+p2+p3+p4)
      
      if(tmp5$Ri[n]==1) {Ubeta = uleft+uright; dUbeta = duleft+duright} else{Ubeta = uright; dUbeta = duright;}
      
      U = U + Ubeta
      dU = dU + dUbeta
    }
    diff = max(solve(-dU) %*% U)
    beta = beta - solve(-dU) %*% U
    #r=r+1
    #cat(r, "\n")
  }
  mybetad = beta
  t(mybetad)
  #outtmp = c(mybetad)
  #out = rbind(out, outtmp)
  #cat(m, "\n")
}

test = foreach(i=1:1000) %dopar% myfunc(i)
test2 = do.call("rbind", test)

write.csv(test2,"aipw_stratD.csv")


stopCluster(cl)