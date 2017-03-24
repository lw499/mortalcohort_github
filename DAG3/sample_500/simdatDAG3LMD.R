##Code to simulate data under DAG3

library(geepack);library(norm);library(gee);library(nlme); library(mvtnorm)

pmiss <- function(data) sapply(data,function(x) data.frame(nmiss=sum(is.na(x)), n=length(x), pmiss=sum(is.na(x))/length(x)))

logit <- function(term) {
  return( log(term/(1-term)) )
}

EXPIT <- function(term) {
  return( exp(term)/(1+exp(term)) )
}

gendata = function(N) {
  shiftagesd=2.613984 ; 
  edusd = 2.336318;
  shiftage = rnorm(N, 3.576063, shiftagesd); 
  sex = rbinom(N,1,0.5); 
  edu = rnorm(N,5.19222,edusd);   # Age, sex, Edu
  t = c(0,2,4,6,8);
  s = rnorm(N, 7.492475, 3.714435)
  s1=ifelse(s<2,0,8);s1=ifelse(s>=2 & s<4, 2, s1);
  s1=ifelse(s>=4 & s<6, 4, s1); s1=ifelse(s>=6 & s<8, 6, s1);
  
  y1 = rnorm(N, 17.63033,3.132368*0.75)
  
  b01=5.05421; b11=0.73159; b21=-0.56623; b41=0.05197; e1=1.85;
  b02=5.16183; b22=0.64555; b32=-0.83485; b52=0.19805; e2=1.804;
  b03=3.95880; b33=0.74143; b43=-0.77977; b63=0.03094; e3=1.793; 
  b04=4.20122; b44=0.71654; b54=-0.02033; e4=1.971;
  
  p01=-3.09988; p11=0.18008; p21=-0.35880; p31=0.24031;
  p02=-6.96815; p22=0.25176; p32=0.16995; p42=0.52376;
  p03=-6.45874; p33=0.14425; p43=0.38340; p53=0.71365;
  p04=-3.43201; p44=0.25073; p54=0.54346;
  
  y2 =  b01 + y1*b11 + sex*b21 + s1*b41 + rnorm(N,0,e1)
  y3 =  b02 + y2*b22 + sex*b32 + s1*b52 + rnorm(N,0,e2)
  y4 =  b03 + y3*b33 + sex*b43 + s1*b63 + rnorm(N,0,e3)
  y5 =  b04 + y4*b44 + sex*b54 + rnorm(N,0,e4)
  edu = ifelse(edu<0, 0, edu)
  Case = seq(1,N,1)
  dat = as.data.frame(cbind(Case, edu, sex, shiftage, y1, y2, y3, y4, y5 ,s1))
  
  dat$dop2 = EXPIT(p01+p11*dat$y1+p21*sex+p31*s1)
  dat$dop3 = EXPIT(p02+p22*dat$y2+p32*sex+p42*s1)
  dat$dop4 = EXPIT(p03+p33*dat$y3+p43*sex+p53*s1)
  dat$dop5 = EXPIT(p04+p44*dat$y4+p54*sex)  
  
  dat$tmpd0=1;
  
  tmp = ifelse(dat$s1>=2,1,0); dat$tmp=tmp; dat1 = dat[tmp==0,]; dattmp = dat[tmp==1,];
  tmpd = ifelse(dattmp$tmp==1 & runif(nrow(dattmp))<=dattmp$dop2, 1,0); dattmp$tmpd=tmpd; 
  ##split data sets
  tmp2 = ifelse(dattmp$s1>=4,1,0); dattmp$tmp2=tmp2;
    tmpd2= ifelse(dattmp$tmpd==1 & dattmp$tmp2==1 & runif(nrow(dattmp))<=dattmp$dop3,1,0);dattmp$tmpd2=tmpd2;dat2=dattmp[dattmp$tmp2==0,];dattmp = dattmp[dattmp$tmp2==1,];
  ##split data sets
  tmp3 = ifelse(dattmp$s1>=6,1,0); dattmp$tmp3=tmp3;
    tmpd3= ifelse(dattmp$tmpd2==1 & dattmp$tmp3==1 & runif(nrow(dattmp))<=dattmp$dop4,1,0);dattmp$tmpd3=tmpd3;dat3=dattmp[dattmp$tmp3==0,];dattmp=dattmp[dattmp$tmp3==1,];
  ##split data sets
  tmp4 = ifelse(dattmp$s1>=8,1,0); dattmp$tmp4=tmp4;
    tmpd4= ifelse(dattmp$tmpd3==1 & dattmp$tmp4==1 & runif(nrow(dattmp))<=dattmp$dop5,1,0); dattmp$tmpd4=tmpd4;
  dat4 = dattmp[dattmp$tmp4==0,]; dat5=dattmp[dattmp$tmp4==1,];
  dat1$dop2=dat1$dop3=dat1$dop4=dat1$dop5=dat1$D2=dat1$D3=dat1$D32=dat1$D4=dat1$D42=dat1$D43=dat1$D5=dat1$D52=dat1$D53=dat1$D54=dat1$tmp=NULL;
  dat1$tmpd=0;dat1$tmpd2=0;dat1$tmpd3=0;dat1$tmpd4=0;
  dat2$dop2=dat2$dop3=dat2$dop4=dat2$dop5=dat2$D2=dat2$D3=dat2$D32=dat2$D4=dat2$D42=dat2$D43=dat2$D5=dat2$D52=dat2$D53=dat2$D54=dat2$tmp=dat2$tmp2=NULL;
  dat2$tmpd3=0;dat2$tmpd4=0;
  dat3$dop2=dat3$dop3=dat3$dop4=dat3$dop5=dat3$D2=dat3$D3=dat3$D32=dat3$D4=dat3$D42=dat3$D43=dat3$D5=dat3$D52=dat3$D53=dat3$D54=dat3$tmp=dat3$tmp2=dat3$tmp3=NULL;
  dat3$tmpd4=0;
  dat4$dop2=dat4$dop3=dat4$dop4=dat4$dop5=dat4$D2=dat4$D3=dat4$D32=dat4$D4=dat4$D42=dat4$D43=dat4$D5=dat4$D52=dat4$D53=dat4$D54=dat4$tmp=dat4$tmp2=dat4$tmp3=dat4$tmp4=NULL;
  dat5$dop2=dat5$dop3=dat5$dop4=dat5$dop5=dat5$D2=dat5$D3=dat5$D32=dat5$D4=dat5$D42=dat5$D43=dat5$D5=dat5$D52=dat5$D53=dat5$D54=dat5$tmp=dat5$tmp2=dat5$tmp3=dat5$tmp4=NULL;
  dat = rbind(dat1,dat2,dat3,dat4,dat5); 
  
  longdat = reshape(dat, varying=list(c(5:9),c(11:15)), v.names = c("pek","R"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
  longdatsort = longdat[order(longdat$id),]
  dodat = longdatsort;
  #dodat = dodat[dodat$s1>=dodat$time,]
  dodat$observedy = ifelse(dodat$R==1,dodat$pek,NA)
  
  s1 = dodat$s1; dodat$s1=NULL; dodat$s1=s1; ##rearranging columns
  dodat$pek = dodat$observedy; dodat$observedy=NULL; dodat$id = NULL; dodat$R=NULL;  
  return(dodat)
}
