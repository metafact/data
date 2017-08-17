vix= read.csv(file = "R_^VIX.csv", header = T, sep = ',', dec = '.')
ads = read.csv(file = "R_FRBP-ADS_VINTAGES_MOSTRECENT_DAILY.csv", header = T, sep = ',', dec = '.')
spare = read.csv(file = "R_OPEC-SA_SPARE.csv", header = T, sep = ';', dec = '.')
cushing = read.csv(file = "R_Weekly_Cushing_OK_Ending_Stocks_excluding_SPR_of_Crude_Oil.csv", header = T, sep = ',', dec = '.')
usprod = read.csv(file = "R_Weekly_U.S._Field_Production_of_Crude_Oil.csv", header = T, sep = ',', dec = '.')
tindex = read.csv(file = "R_Working T_Weekly.csv", header = T, sep = ';', dec = '.')
moments = read.csv(file = "Moments.csv", header = T, sep = ',', dec = '.')

test = vector(length = 10 )
test[seq(1, length(test), 4)] = 1 #Inputs every forth element with '1'
test[seq(2, length(test), 4)] = 4
test[seq(3, length(test), 4)] = 3
test[seq(4, length(test), 4)] = 2

#########

revDF = cushing[rev(rownames(cushing)),] #REVERSES ORDER OF ROWS FOR WHOLE DATA FRAME
cushD = sapply(revDF[1], as.character) #OUTPUTS/CONVERTS MATRIX/ARRAY of CHARACTERS


#EXTRACTS YEAR ONLY
DateTeller = vector(length = length(cushD))
for (i in 1:length(cushD))
{
  if ( nchar(cushD[i]) < 10 )
  {
    DateTeller[i] = substring(cushD[i], 6,9)
  }
  else if  ( nchar(cushD[i]) == 10 )
  {
    DateTeller[i] = substring(cushD[i], 7,10)
  }
}



lim = paste(DateTeller,'-', substring(cushD,1,2), sep = '')
tell = count(lim)

#Empty array
dfTeller = c() 


#To ORDER WEEKS IN EACH YEAR
for (i in 1:length(tell$freq))
{
  temp = 1:1:tell$freq[i]
  dfTeller = c(dfTeller,temp)
}

#MERGES TWO COLUMNS INTO MATRIX
var = cbind(cushD, dfTeller)

ukeN = paste(lim,'-',dfTeller, sep = '')
cushing$ukeN = rev(ukeN)

####!!!!!!!!!!!!!
revDF_us = usprod[rev(rownames(usprod)),] #REVERSES ORDER OF ROWS FOR WHOLE DATA FRAME
cushD_us = sapply(revDF_us[1], as.character) #OUTPUTS/CONVERTS MATRIX/ARRAY of CHARACTERS


#EXTRACTS YEAR ONLY
remove(DateTeller_us)
DateTeller_us = vector(length = length(cushD_us))
for (i in 1:length(cushD_us))
{
  if ( nchar(cushD_us[i]) < 10 )
  {
    DateTeller_us[i] = substring(cushD_us[i], 6,9)
  }
  else if  ( nchar(cushD_us[i]) == 10 )
  {
    DateTeller_us[i] = substring(cushD_us[i], 7,10)
  }
}

#MERGES YEAR AND MONTH
remove(lim_us, tell_us)
lim_us = paste(DateTeller_us,'-', substring(cushD_us,1,2), sep = '')
tell_us = count(lim_us)

#Empty array
dfTeller_us = c() 


#To ORDER WEEKS IN EACH YEAR
for (i in 1:length(tell_us$freq))
{
  temp_us = 1:1:tell_us$freq[i]
  dfTeller_us = c(dfTeller_us,temp_us)
}

#MERGES TWO COLUMNS INTO MATRIX
var_us = cbind(cushD_us, dfTeller_us)

ukeN_us = paste(lim_us,'-',dfTeller_us, sep = '')
usprod$ukeN = rev(ukeN_us)

####!!!!!!!!

revDF_tdx = tindex[rev(rownames(tindex)),] #REVERSES ORDER OF ROWS FOR WHOLE DATA FRAME
cushD_tdx = sapply(revDF_tdx[1], as.character) #OUTPUTS/CONVERTS MATRIX/ARRAY of CHARACTERS


#EXTRACTS YEAR ONLY
remove(DateTeller_tdx)
DateTeller_tdx = substring(cushD_tdx,1,4)


#MERGES YEAR AND MONTH
remove(lim_tdx, tell_tdx)
lim_tdx = paste(DateTeller_tdx,'-', substring(cushD_tdx,6,7), sep = '')
tell_tdx = count(lim_tdx)

#Empty array
dfTeller_tdx = c() 


#To ORDER WEEKS IN EACH YEAR
for (i in 1:length(tell_tdx$freq))
{
  temp_tdx = 1:1:tell_tdx$freq[i]
  dfTeller_tdx = c(dfTeller_tdx,temp_tdx)
}

#MERGES TWO COLUMNS INTO MATRIX
var_tdx = cbind(cushD_tdx, dfTeller_tdx)

ukeN_tdx = paste(lim_tdx,'-',dfTeller_tdx, sep = '')
tindex$ukeN = rev(ukeN_tdx)

##!!!!!
#Reduce(intersect, list(a,b,c))

a = usprod$ukeN
b = tindex$ukeN
i = match(a,b) #returns indicies of b array, by looking each a array element in b array
t1 = subset(i, (!is.na(i))) #returns only non NA elements

naa = which(!is.na(i)) #returns INDICIES of non NA elements
a[naa] == b[t1]  # Should be equal

#Making new array that is intersection btw a and b
c = cushing$ukeN
temp = a[naa] 

#Matching with third array now
j = match(temp,c)
t2 = subset(j, (!is.na(j)))
naCT = which(!is.na(j))

temp[naCT] == c[t2]

######
cushing$ComDate = sapply(cushing[1], as.character)
vix$ComDate = sapply(vix[1], as.character)
ads$ComDate = sapply(ads[1], as.character)
cushing_cut = cushing[t2,]

#Preparing Common Date column for "cushing"  data frame
ar = substring(cushing[t2,3],1,4)
maned = substring(cushing[t2,4],1,2)
#https://hotkeyit.github.io/v2/docs/misc/RegEx-QuickRef.htm
dag = regmatches(cushing[t2,4],regexpr('/[0-9]+', cushing[t2,4]))
dag = substring(dag,2, )

cushing_cut$ComDate = paste(ar,'-', maned,'-',dag, sep = '')

#FIND match/Intersection between ADS and VIX variables
vix_ads = match(vix$ComDate, ads$ComDate)
nonan_ADS = subset(vix_ads, (!is.na(vix_ads)))
naVX = which(!is.na(vix_ads))
TF = vix$ComDate[naVX] == ads$ComDate[nonan_ADS] # TO CHECK if they are equal 

#FIND match/Intersection between ADS and VIX and CUSHING variables now
cushADS = match(cushing_cut$ComDate,ads$ComDate[nonan_ADS])
nonanC_Ads = subset(cushADS, (!is.na(cushADS)))
naCUSH = which(!is.na(cushADS))

adsInd = ads[nonan_ADS,]

TF_CuAd = cushing_cut$ComDate[naCUSH] == adsInd$ComDate[nonanC_Ads]

ndf = data.frame(cushing_cut$ComDate[naCUSH], cushing_cut$KBBL[naCUSH], adsInd$Value[nonanC_Ads], vix$Close[naVX][nonanC_Ads],
                 usprod$KBBLD[naa][naCT][naCUSH],  usprod$ukeN[naa][naCT][naCUSH],tindex$ukeN[t1][naCT][naCUSH], tindex$Chosen.T[t1][naCT][naCUSH])


####
#Preparing dates for MOMENTS df
mom_dates = gsub("Jan", "01", moments$DateStr)
mom_dates = gsub("Feb", "02", mom_dates)
mom_dates = gsub("Mar", "03", mom_dates)
mom_dates = gsub("Apr", "04", mom_dates)
mom_dates = gsub("May", "05", mom_dates)
mom_dates = gsub("Jun", "06", mom_dates)
mom_dates = gsub("Jul", "07", mom_dates)
mom_dates = gsub("Aug", "08", mom_dates)
mom_dates = gsub("Sep", "09", mom_dates)
mom_dates = gsub("Oct", "10", mom_dates)
mom_dates = gsub("Mov", "11", mom_dates)
mom_dates = gsub("Dec", "12", mom_dates)

moments$YMD = paste(substring(mom_dates,7,10),'-', substring(mom_dates,4,5), '-', substring(mom_dates,1,2), sep = '')

###
colnames(ndf)[1] = "Date" 
                         
mom_ndf_i = match(ndf$Date, moments$YMD)    
mom_nonan = subset(mom_ndf_i, (!is.na(mom_ndf_i)))
ndf_nonan = which(!is.na(mom_ndf_i))

#TO check if they are equal
ndf$Date[ndf_nonan] == moments$YMD[mom_nonan]

d2 = as.Date(moments$YMD[mom_nonan])
d1 = as.Date(moments$YMD[mom_nonan-1])
diff_i = d2 -d1

#Taking one day lag in order to omit endogeneity 
ny_ndf = ndf[ndf_nonan,]
ny_ndf$VOL_One_Day_Future = moments$Vol[mom_nonan + 1]


#TTM
var1 = -1*(moments$TTM[mom_nonan + 1])[1:(s2-1)]

#ADS is in %, so take simple difference
s2 = length( ny_ndf$adsInd.Value.nonanC_Ads.)
var2 = ny_ndf$adsInd.Value.nonanC_Ads.[1:(s2 - 1)] - ny_ndf$adsInd.Value.nonanC_Ads.[ 2 : s2 ] 

# US OIL Output
s3 = s2
var3 = log( ny_ndf$usprod.KBBLD.naa..naCT..naCUSH.[1:(s3-1)] / ny_ndf$usprod.KBBLD.naa..naCT..naCUSH.[2:s3])

#
fit =lm(ny_ndf$VOL_One_Day_Future[1:(s2-1)] ~ moments$Vol[mom_nonan][1:(s2-1)] + var1 + var2 + var3 )
summary(fit)
