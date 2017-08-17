#Usefull link for parallel [http://blog.aicry.com/r-parallel-computing-in-5-minutes/]

registerDoParallel(cores=8)  #to register cores
closeAllConnections()	#To close these cores
getDoParWorkers() #Returns number of workers
 registerDoSEQ() #For sequentional, not parallel
 registerDoParallel() # This is for parallel
 ############
 
 detectCores() #But in so far as it is a useful guideline, function detectCores() tries to determine the number of CPU
 # cores on which R is running
 

#start time
strt<-Sys.time()
numWorkers <- 4
cl <- makeCluster(numWorkers, type = "PSOCK")

res <- parLapply(cl, values, workerFunc)

stopCluster(cl)

print(Sys.time()-strt)

strt<-Sys.time()
res <- lapply(values, workerFunc)

print(Sys.time()-strt)




s = "call"
workerFunc <- function( optV, F, K, mat) {AmericanOptionImpliedVolatility(type="call", value=optV, underlying = F, strike= K, dividendYield=0, riskFreeRate=0.05,
                                                                maturity=mat, volatility=0.4)}
                                                                            
                                                                           
valuevec  = fil[1:10000,c(1:3,6)]
x <- foreach(i=1:3) %do% {
  
  AmericanOptionImpliedVolatility(type="call", value=opt, underlying=F,
                                  strike=K, dividendYield=0, riskFreeRate=0.05,
                                  maturity=time, volatility=0.4) } }

setwd("M:/Matlab/")
fil = read.csv("tableY.csv", header = T, sep = ',', dec = ".")

fil2 = fil
i = which(fil[,3]>60)
fil3 = fil2[-i,]
c = fil3[,2]/fil3[,1]
j = which(c>5)
fil3 = fil3[-j,]

# plot(fil3[,2]/fil[,1], fil[,3])

st = Sys.time()
library(doParallel)
registerDoParallel(cores=4)
vec <- vector()
a =NULL

f = foreach(i=1:nrow(fil3), .errorhandling = 'pass', .packages='RQuantLib') %dopar% {
# a = try({AmericanOptionImpliedVolatility(type="call", value=fil3[i,3], underlying = fil3[i,1], strike= fil3[i,2], dividendYield=0, riskFreeRate=0.05, maturity=fil3[i,6], volatility=0.4)})
AmericanOptionImpliedVolatility(type="call", value=fil3[i,3], underlying = fil3[i,1], strike= fil3[i,2], dividendYield=0, riskFreeRate=0.05, maturity=fil3[i,6], volatility=0.4)
}
print(Sys.time()-st)
closeAllConnections()



library(doParallel)
registerDoParallel(cores=2)
x = foreach(i=1:3) %dopar% {sqrt(i)}





for(i in 1:20){
  cc = try({AmericanOptionImpliedVolatility(type="call", value=fil3[i,3], underlying = fil3[i,1], strike= fil3[i,2], dividendYield=0, riskFreeRate=0.05, maturity=fil3[i,6], volatility=0.4)}, silent = T)
  
  # print(Sys.time()-st)
}


#####02.04.2017

st = Sys.time()
library(doParallel)
registerDoParallel(cores=4)

#3955968
 
vec = foreach(i=1:3955968, .packages='RQuantLib') %dopar% {
  if ( class(F[[i]][1]) == "numeric")
  {
    F[[i]][1]
  }
   else
   {
      NaN
   }
 }
