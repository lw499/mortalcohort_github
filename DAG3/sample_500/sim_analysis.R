fdata = read.csv("true.csv")
n=1000
tmean = c(colMeans(fdata[,2:ncol(fdata)]))
true = tmean

tmpdat = read.csv("miss_methods.csv")

n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
geemeans = mean[1:5]; ipwmeans = mean[6:10]; ipwmeansd = mean[11:15]; aipwmeans = mean[16:20]; 
aipwmeansd = mean[21:25]; limeans = mean[26:30]; limeansd = mean[31:35];
mimeans = mean[36:40]; mimeansd = mean[41:45]; 

gee = tmpdat[,2:6]; ipw = tmpdat[,7:11]; ipwd = tmpdat[,12:16]; aipw = tmpdat[,17:21]; 
aipwd = tmpdat[,22:26]; li = tmpdat[,27:31]; lid = tmpdat[,32:36]; 
mi = tmpdat[,37:41]; mid = tmpdat[,42:46];
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
ipwmeansd2 = mean[1:5]
ipwd2 = tmpdat[,2:6];
ipwdSE2 = apply(ipwd2, 2, var)^0.5

bipwdp2 = 100*(ipwmeansd2-true)/ipwdSE2
bipwd2 = (ipwmeansd2-true)

tmpdat = read.csv("ipwp.csv")
n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
ipwmeansp = mean[1:5]
ipwp = tmpdat[,2:6];
ipwpSE = apply(ipwp, 2, var)^0.5

bipwpp = 100*(ipwmeansp-true)/ipwpSE
bipw_p = (ipwmeansp-true)


tmpdat = read.csv("li_stratD.csv")
n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
limeansds = mean[1:5]
lids = tmpdat[,2:6];
lidSEs = apply(lids, 2, var)^0.5

blidps = 100*(limeansds-true)/lidSEs
blids = (limeansds-true)


tmpdat = read.csv("aipw_stratD.csv")
n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
augmean = mean[1:5]
aug = tmpdat[,2:6];
augSE = apply(aug, 2, var)^0.5

baugp = 100*(augmean-true)/augSE
baug = (augmean-true)


tmpdat = read.csv("mi_stratD.csv")
n=1000
mean = c(colMeans(tmpdat[,2:ncol(tmpdat)]))
mimeans = mean[1:5]
mi = tmpdat[,2:6];
midSE2 = apply(mi, 2, var)^0.5

bmidp2 = 100*(mimeans-true)/midSE2
bmid2 = (mimeans-true)



cbind(bgee,bipw,bipw_p, bipwd, bipwd2, baipw, baipwd, baug, bli, blid, blids, bmi, bmid, bmid2)

cbind(bgeep,bipwp,bipwpp,bipwdp,bipwdp2,baipwp,baipwdp,baugp,blip,blidp,blidps,bmip,bmidp,bmidp2)

cbind(geeSE,ipwSE,ipwpSE,ipwdSE,ipwdSE2,aipwSE,aipwdSE,augSE,liSE,lidSE,lidSEs,miSE,midSE,midSE2)

