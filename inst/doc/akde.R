## ----  fig.show='hold'---------------------------------------------------
library(ctmm)
data(buffalo)
cilla <- buffalo[[1]]
M0 <- ctmm.fit(cilla) # no autocorrelation timescales
M2 <- ctmm.fit(cilla,ctmm(tau=c(6*24,1)*60^2)) # ~ 6 day and 1 hour autocorrelation timescales

## ----  fig.show='hold', results = "hide"---------------------------------
KD0 <- akde(cilla,M0)
KD2 <- akde(cilla,M2)

## ----  fig.show='hold'---------------------------------------------------
plot(cilla,AKDE=KD0,ylim=c(-14,12)*1000)
title("M0")
plot(cilla,AKDE=KD2,ylim=c(-14,12)*1000)
title("M2")

## ----  fig.show='hold'---------------------------------------------------
summary(KD0)
summary(KD2)

