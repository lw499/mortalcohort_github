library(geepack);library(norm);library(gee);library(nlme); library(mvtnorm)

pmiss <- function(data) sapply(data,function(x) data.frame(nmiss=sum(is.na(x)), n=length(x), pmiss=sum(is.na(x))/length(x)))

logit <- function(term) {
  return( log(term/(1-term)) )
}

EXPIT <- function(term) {
  return( exp(term)/(1+exp(term)) )
}

gendata = function(N) {
  shiftagesd=2.613984 ; edusd = 2.613984;
  iadlsd = 3.563689; mmsesd = 2.619475;
  
  shiftage = rnorm(N, 3.576063, shiftagesd); 
  iadl1 = rnorm(N, 18.67963, iadlsd);
  mmse1 = rnorm(N, 27.40961, mmsesd); mmse1 = ifelse(mmse1>=30, 30, mmse1);
  iadl1 = trunc(iadl1); mmse1 = trunc(mmse1);
  
  Female = rbinom(N,1,0.5); 
  Smoke = rbinom(N, 1, 0.5);
  Educyrs = rnorm(N,5,edusd);   # Age, sex, Edu
  HlthPrv1 = sample(x=c(1,2,3), size=N, replace=TRUE, prob=c(0.5,0.35,0.15))
  Educyrs = ifelse(Educyrs<0, 0, Educyrs)
  
  t = c(0,2,4,6,8);
  #lamb0 = 0.07
  #mus <- lamb0*exp(0.099*shiftage - 0.605*Female)  #0.5
  #s <- rexp(N, rate=mus)
  s = rnorm(N, 7.492475, 3.714435); s=ifelse(s<=0,0,s);
  s1=ifelse(s<2,0,8);s1=ifelse(s>=2 & s<4, 2, s1);
  s1=ifelse(s>=4 & s<6, 4, s1); s1=ifelse(s>=6 & s<8, 6, s1);
  
  y1 = 13.47593 + Female*-3.33444 + Educyrs*0.22453 + shiftage*-0.18782 + Smoke*-0.80261 + HlthPrv1*0.18997 + iadl1*0.18810 + mmse1*0.07734 + I(s1==2)*0.20273 + I(s1==4)*0.16479 + I(s1==6)*-0.39598  + I(s1==8)*0.36652  + rnorm(N,0, 2.592)
  y2 =  5.24791 + y1*0.73569 + Female*-0.57920 + Educyrs*-0.05504 + shiftage*-0.01999 + Smoke*-0.12279 + I(s1==4)*0.72715 + I(s1==6)*0.32057  + I(s1==8)*0.55488 + rnorm(N,0, 1.857)
  y3 =  5.55183 + y1*0.23752 + y2*0.45543 + Female*-1.27132 + Educyrs*0.09970 + shiftage*-0.05039 + Smoke*-1.36687 + I(s1==6)*0.65842  + I(s1==8)*0.75565 + rnorm(N,0, 1.685)
  y4 =  2.48118 + y1*0.11221 + y2*0.20401 + y3*0.51076 + Female*-0.73783 + Educyrs*0.02411 + shiftage*-0.01976 + Smoke*-0.34973 + I(s1==8)*0.02962 + rnorm(N,0, 1.765)
  y5 =  0.28727 + y1*-0.01626 + y2*0.40800 + y3*0.10375 + y4*0.39561 + Female*-0.03997 + Educyrs*0.06337  + shiftage*-0.01733 + Smoke*-0.50521 + rnorm(N,0, 1.815)
  
  #y1 = 13.47593 + Female*-3.33444 + Educyrs*0.22453 + shiftage*-0.18782 + Smoke*-0.80261 + HlthPrv1*0.18997 + iadl1*0.18810 + mmse1*0.07734 + I(s1==2)*0.20273 + I(s1==4)*0.16479 + I(s1==6)*-0.39598  + I(s1==8)*0.36652  + rnorm(N,0, 2.592)
  #y2 =  5.47015 + y1*0.73862 + Female*-0.53857 + Educyrs*-0.06005 + shiftage*-0.02068 + Smoke*-0.12287 + HlthPrv1*-0.20243 + iadl1*-0.02603 + mmse1*0.02050 + I(s1==4)*0.69798 + I(s1==6)*0.35334  + I(s1==8)*0.50606 + rnorm(N,0, 1.857)
  #y3 =  3.83324 + y1*0.23514 + y2*0.45281 + Female*-1.29480 + Educyrs*0.08695 + shiftage*-0.04838 + Smoke*-1.39704 + HlthPrv1*-0.09044 + iadl1*-0.02369 + mmse1*0.09021 + I(s1==6)*0.67878  + I(s1==8)*0.70536 + rnorm(N,0, 1.685)
  #y4 =  4.867486 + y1*0.103525 + y2*0.212708 + y3*0.513147 + Female*-0.683705 + Educyrs*0.033072 + shiftage*-0.006759 + Smoke*-0.383677 + HlthPrv1*-0.076515 + iadl1*-0.070334 + mmse1*-0.034318 + I(s1==8)*0.07 + rnorm(N,0, 1.765)
  #y5 =  -5.00426 + y1*-0.01109 + y2*0.41938 + y3*0.07031 + y4*0.39111 + Female*-0.18020 + Educyrs*0.02766  + shiftage*-0.03773 + Smoke*-0.51568 + HlthPrv1*-0.26042 + iadl1*0.06612 + mmse1*0.17775 + rnorm(N,0, 1.815)

  Case = seq(1,N,1)
  dat = as.data.frame(cbind(Case, Educyrs, Female, shiftage, y1, y2, y3, y4, y5 ,s1, Smoke, HlthPrv1, iadl1, mmse1))
  
  dat$dop2 = EXPIT(-5.713027 + y1*0.136817 + Female*-0.373208 + Educyrs*0.037519 + shiftage*-0.087549 + HlthPrv1*-0.760102 + iadl1*0.002976 + mmse1*0.195774 + I(s1==4)*0.934491 + I(s1==6)*1.181298  + I(s1==8)*1.152166)
  dat$dop3 = EXPIT(-11.34620 + y2*0.23503 + Female*-0.22591 + Educyrs*-0.07684 + shiftage*0.06565 + HlthPrv1*-0.02523 + iadl1*0.18696 + mmse1*0.13019 + I(s1==6)*1.74180  + I(s1==8)*2.06026)
  dat$dop4 = EXPIT(1.51794 + y3*0.13960  + Female*0.45628 + Educyrs*0.14697 + shiftage*-0.10833 + HlthPrv1*-0.31985 + iadl1*0.02073 + mmse1*-0.11143 + I(s1==8)*1.39129)
  dat$dop5 = EXPIT(0.56406 + y4*0.28252 + Female*0.67533 + Educyrs*0.04140 + shiftage*-0.05918 + HlthPrv1*0.71468 + iadl1*0.04064 + mmse1*-0.22738)  

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
  
  longdat = reshape(dat, varying=list(c(5:9),c(15:19)), v.names = c("pek","R"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
  longdatsort = longdat[order(longdat$id),]
  dodat = longdatsort;
  #dodat = dodat[dodat$s1>=dodat$time,]
  dodat$observedy = ifelse(dodat$R==1,dodat$pek,NA)
  
  s1 = dodat$s1; dodat$s1=NULL; dodat$s1=s1; ##rearranging columns
  dodat$pek = dodat$observedy; dodat$observedy=NULL; dodat$id = NULL; dodat$R=NULL;  
  return(dodat)
}
