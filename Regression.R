fil = read.delim("clipboard", dec = ",")
fil = read.csv("Data_1.csv", header = T, dec = ",", sep = ";")
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
oil = read.csv("Oil.csv", header = T, dec = ",", sep = ";")
r2 = merge(r, oil, by.x = "Year", by.y = "Date", all = T)
r2$StkPrem = r2$PX_LAST -r2$Tbills 
x = which(r2$TOTAL_PRODUCTION_MMBOE == "Inf")
r2 = r2[-x,]
x = which(r2$DEVEL_WELLS_DRILLED_NET == "Inf")
r2 = r2[-x,]
fit = lm(StkPrem ~ MktPremium + Oil.Return + CF_CASH_FROM_OPER+BOE_RESERVES_END_YR_US+TOTAL_PRODUCTION_MMBOE, data = r2,  na.action= na.exclude)            
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

#EIA data M:\Documents\MATLAB
data = read.delim("clipboard", header = T, dec = ",")
fit = lm(data$Weighted~data$MktPremium+data$ProdRet +data$ResRet)


#Find and replace in Data frame
df[df == 11929] = 1010101 #Has to be numeric then
##########################
#REGRESSION
#########
summary(fit)$coefficients[1,4] #extracts p-value of the first coefficient #$coef
coef(fit)

#Variance Inflation
library(car)
vif(fit)

qplot(x,y, color = "red", geom = c("point, "smooth")) #ggplots where "geom" smoothes points 
#or adds a trend, we need data frame for ggplot

g = qplot(Stk$Prem, r3$BOE_RESERVES_END_YR_US, color = "red")
g1 = g +geom_smooth(method = "lm", color = "blue") #Adds regression line to points



