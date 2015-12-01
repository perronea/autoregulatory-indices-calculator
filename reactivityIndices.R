#Autoregulatory Indices

#read in data from .txt and .csv
ICP = read.table("~/Desktop/Patient_Identifier/mobius/Patient_Identifier_ICP.txt", sep="\t", skip=1) #length=2119
ABP = read.table("~/Desktop/Patient_Identifier/mobius/Patient_Identifier_ABP.txt", sep="\t", skip=1) #length=271232
MAP = read.table("~/Desktop/Patient_Identifier/mobius/Patient_Identifier_MAP.txt", sep="\t", skip=1) #length=2119
CPP <- MAP$V2-ICP$V2
ECG = read.table("~/Desktop/Patient_Identifier/mobius/Patient_Identifier_ECG.txt", sep="\t", skip=1, nrows=4000) #length=1084928
totalH = read.table("~/Desktop/Patient_Identifier/nirs/Patient_Identifier_Total.csv", sep=",", skip=41)
meanTotalH <- rowMeans(totalH[2:25])
oxyH = read.table("~/Desktop/Patient_Identifier/nirs/Patient_Identifier_Oxy.csv", sep=",", skip=41)
meanOxyH <- rowMeans(oxyH[2:25])
percentTOI <- (meanOxyH/meanTotalH)*100
percentTOI[percentTOI>100] <- 100
percentTOI[percentTOI<0] <- 0

#converts a time in str format to seconds
timeToSec <- function(time) {
  secVec <- strsplit(as.character(time), ":")
  if(length(unlist(secVec)) == 4) {
    secs <- as.numeric(secVec[[1]][2])*3600+as.numeric(secVec[[1]][3])*60+as.numeric(secVec[[1]][4])
  }
  else if(length(unlist(secVec)) == 3) {
    secs <- (as.numeric(secVec[[1]][1])-2)*3600+as.numeric(secVec[[1]][2])*60+as.numeric(secVec[[1]][3]) #time off set by 2 hours
  }
  else {print("time not recognized")}
}
#ICP,ABP,TOI,MAP,CPP,ECG
ICPtime <- sapply(ICP$V1, timeToSec)
ABPtime <- sapply(ABP$V1, timeToSec)
TOItime <- sapply(totalH$V27, timeToSec)
MAPtime <- sapply(MAP$V1, timeToSec)
CPPtime <- MAPtime
ECGtime <- sapply(ECG$V1, timeToSec)

#calculate moving average
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
#eliminate na values
maICP <- na.omit(ma(ICP$V2, 10))
maICPtime <- na.omit(ma(ICPtime, 10))
maABP <- na.omit(ma(ABP$V2, 10))
maABPtime <- na.omit(ma(ABPtime, 10))
maTOI <- na.omit(ma(percentTOI, 10))
maTOItime <- na.omit(ma(TOItime, 10))
maMAP <- na.omit(ma(MAP$V2, 10))
maMAPtime <- na.omit(ma(MAPtime, 10))
maCPP <- na.omit(ma(CPP, 10))
maCPPtime <- na.omit(ma(CPPtime, 10))
maECG <- na.omit(ma(ECG$V2, 10))
maECGtime <- na.omit(ma(ECGtime, 10))

#align start time
startTime <- function(vecA, vecB, vecA1, vecB1) {
  i <- 1
  if(vecA[1] > vecB[1]) {
    while(vecB[i] < vecA[1]) {
      i <- i+1
    }
    vecB <- vecB[-(1:i)]
    vecB1 <- vecB1[-(1:i)]
  }
  else {
    while(vecA[i] < vecB[1]) {
      i <- i+1
    }
    vecA <- vecA[-(1:i)]
    vecA1 <- vecA1[-(1:i)]
  }
  return(list(vecA, vecB, vecA1, vecB1))
}
#must be returned as a list and unpacked into desired variables
#Prx
alignedStartPrx <- startTime(maICPtime, maABPtime, maICP, maABP)
startICPtimePrx <- unlist(alignedStartPrx[1])
startABPtimePrx <- unlist(alignedStartPrx[2])
startICPPrx <- unlist(alignedStartPrx[3])
startABPPrx <- unlist(alignedStartPrx[4])
#Toxa
alignedStartToxa <- startTime(maABPtime, maTOItime, maABP, maTOI)
startABPtimeToxa <- unlist(alignedStartToxa[1])
startTOItimeToxa <- unlist(alignedStartToxa[2])
startABPToxa <- unlist(alignedStartToxa[3])
startTOIToxa <- unlist(alignedStartToxa[4])
#Tox
alignedStartTox <- startTime(maTOItime, maCPPtime, maTOI, maCPP)
startTOItimeTox <- unlist(alignedStartTox[1])
startCPPtimeTox <- unlist(alignedStartTox[2])
startTOITox <- unlist(alignedStartTox[3])
startCPPTox <- unlist(alignedStartTox[4])

#calculates the mean of window size from one matching value to the next
#error if vecA is longer than vecB in certain cases must trouble shoot
align <- function(vecA, vecB, vecA1, vecB1) {
  if(vecA[2]-vecA[1] > vecB[2]-vecB[1]) {
    mvecB <- c(vecB[1])
    mvecB1 <- c(vecB1[1])
    b <- 1
    for(a in 2:length(vecA)) {
      window <- c()
      window1 <- c()
      while(vecB[b] < vecA[a] & b < length(vecB)) {
        window <- c(window, vecB[b])
        window1 <- c(window1, vecB1[b])
        b <- b+1
      }
      mvecB <- c(mvecB, mean(window))
      mvecB1 <- c(mvecB1, mean(window1))
    }
#    vecA <- head(vecA, -(length(vecA)-length(mvecB))) #trim excess elements on vecA
#    vecA1 <- head(vecA1, -(length(vecA1)-length(mvecB1)))
    return(list(vecA, mvecB, vecA1, mvecB1))
  }
  else {
    mvecA <- c(vecA[1])
    mvecA1 <- c(vecA1[1])
    a <- 1
    for(b in 2:length(vecB)) {
      window <- c()
      window1 <- c()
      while(vecA[a] < vecB[b] & a < length(vecA)) {
        window <- c(window, vecA[a])
        window1 <- c(window1, vecA1[a])
        a <- a+1
      }
      mvecA <- c(mvecA, mean(window))
      mvecA1 <- c(mvecA1, mean(window1))
    }
#    vecB <- head(vecB, -(length(vecB)-length(mvecA))) #trim excess elements on vecA
#    vecB1 <- head(vecB1, -(length(vecB1)-length(mvecA1)))    
    return(list(mvecA, vecB, mvecA1, vecB1))
  }
}
#Prx
alignedPrx <- align(startICPtimePrx, startABPtimePrx, startICPPrx, startABPPrx)
alignedICPtimePrx <- unlist(alignedPrx[1])
alignedABPtimePrx <- unlist(alignedPrx[2])
alignedICPPrx <- unlist(alignedPrx[3])
alignedABPPrx <- unlist(alignedPrx[4])
#Toxa
alignedToxa <- align(startABPtimeToxa, startTOItimeToxa, startABPToxa, startTOIToxa)
alignedABPtimeToxa <- unlist(alignedToxa[1])
alignedTOItimeToxa <- unlist(alignedToxa[2])
alignedABPToxa <- unlist(alignedToxa[3])
alignedTOIToxa <- unlist(alignedToxa[4])
#Tox
alignedTox <- align(startCPPtimeTox, startTOItimeTox, startCPPTox, startTOITox)
alignedCPPtimeTox <- unlist(alignedTox[1])
alignedTOItimeTox <- na.omit(unlist(alignedTox[2]))
alignedCPPTox <- unlist(alignedTox[3])
alignedTOITox <- na.omit(unlist(alignedTox[4]))
alignedCPPtimeTox <- head(alignedCPPtimeTox,-565) #Trim excess data
alignedCPPTox <- head(alignedCPPTox,-565)

#calculates a moving correlation from a specified lagging window size incrementing 1 data point each time
movingCorr <- function(vecA, vecB, window=6) {
  corr <- c(rep(NA, window))
  for(i in (window+1):length(vecB)) {
    windowA <- c(vecA[(i-window):(i-1)])
    windowB <- c(vecB[(i-window):(i-1)])
    corr <- c(corr, cor(windowA, windowB))
  }
  corr
}
#Prx is moving correlation between ICP and ABP
Prx <- movingCorr(alignedICPPrx, alignedABPPrx, window=10)
#replace all NA(arises when there is no change in values in one window) with 0
Prx[is.na(Prx)] <- 0
#Toxa is moving correlation between ABP and TOI
Toxa <- movingCorr(alignedABPToxa, alignedTOIToxa, window=100)
Toxa[is.na(Toxa)] <- 0
#Tox moving corr between CPP and TOI
Tox <- movingCorr(alignedCPPTox, alignedTOITox, window=100)
Tox[is.na(Tox)] <- 0

#plots
par(mfrow=c(2,1)) #plot layout
#Raw Data
plot(maICPtime, maICP, main="ICP", xlab="time(s)", ylab="ICP(mmHg)", type='l')
plot(maABPtime, maABP, main="ABP", xlab="time(s)", ylab="ABP(mmHg)", type='l')
plot(maTOItime, maTOI, main="TOI", xlab="time(s)", ylab="TOI(%)", type='l')
plot(maCPPtime, maCPP, main="CPP", xlab="time(s)", ylab="CPP(mmHg)", type='l')
#Prx
plot(alignedICPtimePrx, alignedICPPrx, main="Aligned ICP", xlab="time(s)", ylab="ICP(mmHg)", type='l')
plot(alignedABPtimePrx, alignedABPPrx, main="Aligned ABP", xlab="time(s)", ylab="ABP(mmHg)", type='l')
plot(alignedABPtimePrx, Prx, main="PRx (window=100)", xlab="time(s)", ylab="PRx", type='l')
#Toxa
plot(alignedTOItimeToxa, alignedTOIToxa, main="Aligned TOI", xlab="time(s)", ylab="TOI(%)", type='l')
plot(alignedABPtimeToxa, alignedABPToxa, main="Aligned ABP", xlab="time(s)", ylab="ABP(mmHg)", type='l')
plot(alignedTOItimeToxa, Toxa, main="TOxA (window=100)", xlab="time(s)", ylab="TOxA", type='l')
#Tox
plot(alignedCPPtimeTox, alignedCPPTox, main="Aligned CPP", xlab="time(s)", ylab="CPP(mmHg)", type='l')
plot(alignedTOItimeTox, alignedTOITox, main="Aligned TOI", xlab="time(s)", ylab="TOI(%)", type='l')
plot(alignedCPPtimeTox, Tox, main="TOx (window=100)", xlab="time(s)", ylab="TOx", type='l')


plot(maECGtime, maECG, main="ECG", xlab="time(s)", ylab="ECG(mV)", type='l')
plot(ECG$V1, ECG$V2, main="ECG", xlab="time(s)", ylab="ECG(mV)", type='l')


