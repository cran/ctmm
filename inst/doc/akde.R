## ----  fig.show='hold'---------------------------------------------------
library(ctmm)
data(buffalo)
cilla <- buffalo[[1]]
M0 <- ctmm.fit(cilla) # no autocorrelation timescales
m2 <- ctmm(tau=c(6*24,1)*60^2) # ~ 6 day and 1 hour autocorrelation timescales
M2 <- ctmm.fit(cilla,m2) 

## ----  fig.show='hold', results = "hide"---------------------------------
UD0b <- akde(cilla,M0,debias=FALSE)
UD0 <- akde(cilla,M0)
UD2b <- akde(cilla,M2,debias=FALSE)
UD2 <- akde(cilla,M2)

## ----  fig.show='hold'---------------------------------------------------
plot(cilla,UD=UD0b,ylim=c(-14,12)*1000)
title("IID KDE")
plot(cilla,UD=UD2b,ylim=c(-14,12)*1000)
title("OUF AKDE")
plot(cilla,UD=UD0,ylim=c(-14,12)*1000)
title(expression("IID KDE"["C"]))
plot(cilla,UD=UD2,ylim=c(-14,12)*1000)
title(expression("OUF AKDE"["C"]))

## ----  fig.show='hold'---------------------------------------------------
summary(UD0)
summary(UD2)

