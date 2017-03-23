fdata = read.csv("true.csv")
n=1000
tmean = c(colMeans(fdata[,2:ncol(fdata)]))
true = tmean

tmpdat = read.csv("miss_methods.csv")

n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
geemeans = mean[1:10]; ipwmeans = mean[11:20]; ipwmeansd = mean[21:30]; aipwmeans = mean[31:40]; 
aipwmeansd = mean[41:50]; limeans = mean[51:60]; limeansd = mean[61:70];
mimeans = mean[71:80]; mimeansd = mean[81:90]

gee = tmpdat[,2:11]; ipw = tmpdat[,12:21]; ipwd = tmpdat[,22:31]; aipw = tmpdat[,32:41]; 
aipwd = tmpdat[,42:51]; li = tmpdat[,52:61]; lid = tmpdat[,62:71]; 
mi = tmpdat[,72:81]; mid = tmpdat[,82:91];
geeSE = apply(gee, 2, var)^0.5
ipwSE = apply(ipw, 2, var)^0.5
ipwdSE = apply(ipwd, 2, var)^0.5
aipwSE = apply(aipw, 2, var)^0.5
aipwdSE = apply(aipwd, 2, var)^0.5
liSE = apply(li, 2, var)^0.5
lidSE = apply(lid, 2, var)^0.5
miSE = apply(mi, 2, var)^0.5
midSE = apply(mid, 2, var)^0.5

bgeep = 100*(geemeans-true)/geeSE
bipwp = 100*(ipwmeans-true)/ipwSE
bipwdp = 100*(ipwmeansd-true)/ipwdSE
baipwp = 100*(aipwmeans-true)/aipwSE
baipwdp = 100*(aipwmeansd-true)/aipwdSE
blip = 100*(limeans-true)/liSE
blidp = 100*(limeansd-true)/lidSE
bmip = 100*(mimeans-true)/miSE
bmidp = 100*(mimeansd-true)/midSE

bgee = (geemeans-true)
bipw = (ipwmeans-true)
bipwd = (ipwmeansd-true)
baipw = (aipwmeans-true)
baipwd = (aipwmeansd-true)
bli = (limeans-true)
blid = (limeansd-true)
bmi = (mimeans-true)
bmid = (mimeansd-true)


tmpdat = read.csv("ipw_stratD.csv")
n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
ipwmeansd2 = mean[1:10]
ipwd2 = tmpdat[,2:11];
ipwdSE2 = apply(ipwd2, 2, var)^0.5

bipwdp2 = 100*(ipwmeansd2-true)/ipwdSE2
bipwd2 = (ipwmeansd2-true)

tmpdat = read.csv("ipwp.csv")
n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
ipwmeansp = mean[1:10]
ipwp = tmpdat[,2:11];
ipwpSE = apply(ipwp, 2, var)^0.5

bipwpp = 100*(ipwmeansp-true)/ipwpSE
bipw_p = (ipwmeansp-true)
#cbind(ipwmeansp, ipwpSE, bipwpp, bipwp)


tmpdat = read.csv("li_stratD.csv")
n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
limeansds = mean[1:10]
lids = tmpdat[,2:11];
lidSEs = apply(lids, 2, var)^0.5

blidps = 100*(limeansds-true)/lidSEs
blids = (limeansds-true)
#cbind(limeansds, lidSEs, blidps, blids)


tmpdat = read.csv("aipw_stratD.csv")
n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
augmean = mean[1:10]
aug = tmpdat[,2:11];
augSE = apply(aug, 2, var)^0.5

baugp = 100*(augmean-true)/augSE
baug = (augmean-true)


tmpdat = read.csv("mi_stratD.csv")
n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
mimeans = mean[1:10]
mi = tmpdat[,2:11];
midSE2 = apply(mi, 2, var)^0.5

bmidp2 = 100*(mimeans-true)/midSE2
bmid2 = (mimeans-true)

rbind(bgee*100,bipw*100,bipw_p*100, bipwd*100, bipwd2*100, baipw*100, baipwd*100, baug*100, bli*100, blid*100, blids*100, bmi*100, bmid*100, bmid2*100)

rbind(bgeep,bipwp,bipwpp,bipwdp,bipwdp2,baipwp,baipwdp,baugp,blip,blidp,blidps,bmip,bmidp,bmidp2)

rbind(geeSE*100,ipwSE*100,ipwpSE*100,ipwdSE*100,ipwdSE2*100,aipwSE*100,aipwdSE*100,augSE*100,liSE*100,lidSE*100,lidSEs*100,miSE*100,midSE*100,midSE2*100)
