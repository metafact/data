rep = integer(100)
j = 1
for (i in 1:(length(t$Year)-1))
{
  if (t$Year[i] == t$Year[i+1])
  {
    rep[j] = i
    j = j+1
  }
}



  

del2 = integer(500)
j = 1
for (i in 1:length(indeks))
{
  a = substring(t$DATE[indeks[i]+1], 7, 10)
  b = substring(t$DATE[indeks[i]], 7, 10)
  if (!is.na(t$PX_LAST[indeks[i]]) & (a==b))
  {
    del2[j] = indeks[i]
    j = j +1
  }
}
-------------------

for (i in 1:length(indeks))
{
  if (is.na(renn2$PX_LAST[indeks[i]+1]))
  {
    renn2$PX_LAST[indeks[i]+1] = renn2$PX_LAST[indeks[i]]
    renn2$PX_LAST[indeks[i]] = NA
  }
 
}




for (i in 1:length(indeks))
{
  vec[i] = unfactor(raw$Year[indeks[i]+1])-1
  raw$Year[indeks[i]] = vec[i]
}

d = data.frame(x =seq(1,10),
               n = c(0,0,1,2,3,4,4,5,6,6),
               logp = signif(-log10(runif(10)), 2))

par(mar = c(5,5,2,5))
with(d, plot(x, logp, type="l", col="red3", 
             ylab=expression(-log[10](italic(p))),
             ylim=c(0,3)))

par(new = T)
with(d, plot(x, n, pch=16, type="l", axes=F, xlab=NA, ylab=NA, cex=1.2)) #cex is for scaling 
axis(side = 4)
mtext(side = 4, line = 3, 'Number genes selected')
legend("topleft",
       legend=c(expression(-log[10](italic(p))), "N genes"),
       lty=c(1,0), pch=c(NA, 0), col=c("red3", "black")) #lty = line type, pch = is for sign of the line in the legends


####################
fil = read.delim("clipboard", dec = ",")
x = grep("([.]0[1-9][.])", fil$DATE, perl = FALSE, value = F)
zap = fil
fil = fil[-x,]
x = grep("(29|30[.])", fil$DATE, perl = FALSE, value = F)
indeks = x
del = integer(1321)
j = 1
for (i in 1:(length(indeks)))
{
  a = substring(fil$DATE[indeks[i]+1], 7, 10)
  b = substring(fil$DATE[indeks[i]], 7, 10)
  
  m = fil$PX_LAST[indeks[i]+1]
  n = fil$PX_LAST[indeks[i]]
  if ( a == b )
  { 
    if (!is.na(n) & is.na(m))
    {
      fil$PX_LAST[indeks[i]+1] = fil$PX_LAST[indeks[i]]
      fil$PX_LAST[indeks[i]] = NA
      del[j] = indeks[i]
      j = j +1
    } else if (is.na(n))
    {
      del[j] = indeks[i]
      j = j +1
    }
    else if (!is.na(n) & !is.na(m))
    {
      del[j] = indeks[i]
      j = j +1
    }
  }
}

fil = fil[-indeks,]
fil2 = fil[,2:45]
d = data.frame ( diff (as.matrix(fil2)))
b = fil2[-nrow(fil2),]
r = d/b 
r$Date =fil$DATE[-1]
r$Dummy = fil$Dummy[-1]
r$Year = substring(r$Date, 7, 10)
x = which(r$Year == "")
r = r[-x,]
r2 = merge(r, oil, by.x = "Year", by.y = "Date", all = T)
r2$StkPrem = r2$PX_LAST -r2$Tbills 
x = which(r2$TOTAL_PRODUCTION_MMBOE == "Inf")
r2 = r2[-x,]
x = which(r2$DEVEL_WELLS_DRILLED_NET == "Inf")
r2 = r2[-x,]
fit = gls(StkPrem ~ MktPremium + Oil.Return + CF_CASH_FROM_OPER+BOE_RESERVES_END_YR_US+TOTAL_PRODUCTION_MMBOE, data = r2,  na.action= na.exclude)            
summary(fit)

#Subperiod for E&P
subE = subset(r2, Dummy == 0)
gitE = gls(StkPrem ~ Oil.Return + MktPremium +BOE_RESERVES_END_YR_US +CF_CASH_FROM_OPER+TOTAL_PRODUCTION_MMBOE, data = subE,  na.action= na.exclude)
summary(gitE)

#Subperiod for Integr
subI = subset(r2, Dummy == 1)
x = which(subI$BOE_END_YEAR_WORLD == "Inf")
subI = subI[-x,]
gitI = gls(StkPrem ~ Oil.Return + MktPremium +BOE_END_YEAR_WORLD +CF_CASH_FROM_OPER+TOTAL_PRODUCTION_MMBOE, data = subI,  na.action= na.exclude)
summary(gitI)

#Subperiod for 2006-2015
sub615 = subset(r2, Year <=2015 & Year >= 2006)
git615 = gls(StkPrem ~ Oil.Return + MktPremium +BOE_RESERVES_END_YR_US +CF_CASH_FROM_OPER+TOTAL_PRODUCTION_MMBOE, data = sub615,  na.action= na.exclude)
summary(git615)





