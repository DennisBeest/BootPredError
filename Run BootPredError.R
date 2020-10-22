

#source("BootPredError.R") #insert file path
source("https://raw.githubusercontent.com/DennisBeest/BootPredError/master/BootPredError.R")

####################################################
#--Simulated some survival data, replace these by real data
####################################################


P <- 50
N <- 100
Train <- matrix(nrow=N,ncol=P)
Train[] <- runif(P*N)
Survival <- rep(c(0,1),each=N/2)
Train[Survival==1,1:5] <- Train[Survival==1,5]+0.25
Time <- runif(N)
Time[Survival==1] <- Time[Survival==1]/2
colnames(Train) <- 1:P
UnPen <- data.frame(U=runif(100))


####################################################
#----Run the bootstrap: Survival outcome for, respectively, penalized, unpenalized+penalized, and unpenalized data
####################################################


example1 <- BootPredError(PenData=Train,Y=Surv(Time,Survival),family="cox")
example2 <- BootPredError(PenData=Train,UnPenData=UnPen,Y=Surv(Time,Survival),family="cox")
example3 <- BootPredError(UnPenData=UnPen,Y=Surv(Time,Survival),family="cox")

example1$Lambda #lambda used
example1$Means 	#oob-auc for the bootstrap, bootstrap632, bootstrap632+, and apparent error
example1$CIs	  #95% CI's of auc estimates


####################################################
#----Run the bootstrap: Binary outcome for, respectively, penalized, unpenalized+penalized, and unpenalized data
####################################################


example4 <- BootPredError(PenData=Train,Y=Survival,family="binomial")
example5 <- BootPredError(PenData=Train,UnPenData=UnPen,Y=Survival,family="binomial")
example6 <- BootPredError(UnPenData=UnPen,Y=Survival,family="binomial")


####################################################
#----More output
####################################################


example1 <- BootPredError(PenData=Train,Y=Surv(Time,Survival),family="cox",outBig=1)

example1$EffectiveBoot    #Safety check. If data is unbalanced data one bootstrap may not contain enough cases, and the AUC is set to NA. EffectiveBoot = NBoot-number of NA's.
example1$InBagSamples     #Matrix containing the indexes of the cases that were in a particular bootstrap
example1$PredictionsOOB   #Average oob-prediction per sample
example1$PredInBag        #InBag predictions, in same order as InBagSamples
example1$PredOOB          #oob predictions per bootstrap, e.g. rowMeans(example1$PredOOB,na.rm=TRUE) == example1$PredictionsOOB
example1$CIBase           #Size of the confidence interval
example1$AUCInBag         #Inbag AUC per bootstrap
example1$RealModel        #This is the "real" model that is being bootstrapped
example1$LPModel          #Linear predictor of the RealModel


