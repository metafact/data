-options(scipen=999) #to disable scientific format of small numbers
-writeClipboard(as.character(xx$pars)) #To copy to clipboard

-plot(temp$Date, temp$Wells, t = "l", xaxt = "n", ylab = "Number of Wells")
-par (new = T)
-plot(temp$xDate temp$Oil, t = "l", col = "red", axes = F, ylab = NA)
-axis(1, at = temp$Date, labels = temp$Date, las = 2)
-axis(4)
-mtext(side = 4, line = 3, "Oil Prices")

-legend("topleft",
  legend=c(expression(-log[10](italic(p))), "N genes"),
  lty=c(1,0), pch=c(NA, 0), col=c("red3", "black")) #lty = line type [1 = line, 2 = dotted line], pch = is for sign of the line in the legends

-grid(10,10)
- abline(v = 2007,h = 60, col = "blue") # for making grids or additional lines

MERGING#VLOOKUP#MATCHING
-dfNew = merge(df1, df2, by.x = "Year", by.y = "Date", all = T) #by.x has a title "Year" while by.y has a title "Date", all = T means that those cells that are not matched                                                                                              #still will appear in a new merged df

DELETE COLUMN in DATA FRAME
-within ( df, rm ( column_Name ) )

REGEX#SUBSTRINGS#
https://regexone.com/

-grep("(.09.)", origin$DATE, perl=TRUE, value=FALSE)  #matches only those that have ".09." and return indexes not values b/c pf                                                                                                 value = FALSE
-substring("abcd",2,3) # returns "bc"

>>>>>>>>>>>>>>>>>>>>>>>ALGO>>>>>>>>>>>>>>>>

- Omit useless columns [,5:40]
- Select data that contains only dates with ."..12.Year" and not "..09.year" for examples
- origin$Year = as.numeric( substring(origin$DATE,7,10) )
- dif = data.frame(as.matrix( diff( df_origin)))
- lag = df_origin[-c(1),]; returns = dif/lag
- newdata <- df[order(df$Year),]
- retu = df[-which(is.na(df$Year)),] # delete emty cells

REGEX
2[1-9]PIPE3[0][.]  #chooses any nuber between 21 and 29 or 30 and then dot at the end.




###############

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
#Delete first all NA or similar texts so that it is not imported as FACTORS
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

oil = read.csv("Oil.csv", header = T, dec = ",", sep = ";") #Read oil data first
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

#EIA data M:\Documents\MATLAB
data = read.delim("clipboard", header = T, dec = ",")
fit = lm(data$Weighted~data$MktPremium+data$ProdRet +data$ResRet)




