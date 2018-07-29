
library(glmnet)
library(risksetROC)
library(survival)
library(boot)

#----References:
#Wahl et al. Assessment of predictive performance in incomplete data by combining internal validation and multiple imputation. BMC Medical Research Methodology, 16:144, (2016).
#Jiang et al. Estimating the condence interval for prediction errors of support vector machine classiers. The Journal of Machine Learning Research 9, 521-540 (2008).
#Jiang et al. Calculating confidence intervals for prediction error in microarray classification using resampling. Stat Appl Genet Mol Biol. 2008;7(1)
#Ref for .632+ bootstrap: Efron et al. Improvement on cross-validation: the 0.632+ bootstrap method. J Am Stat Assoc 92, 548-560 (1997)


##############################################
#--Start function

BootPredError <- function(PenData=NULL,Y,UnPenData=NULL,NBoot=1000,family="cox",outBig=0,alphaIn=0,Lambda=NULL,setseed1=100,setseed2=0,iAUCtmax=5,BootLab=FALSE)
{
	#PenData is the data matrix that needs penalisation
	#UnPenData is the data frame that is unpenalized. By default the matrix is used as it is, dummy variables need to be created before running the file.
	#NBoot = the number of bootstraps
	#family = "cox" or "binomial, as usual with glmnet. If PenData = NULL then glm/coxph are used 
	#alphaIn is the alpha of glmnet
	#Lambda is the lambda for glmnet
	#setseed1 set seeds of cv.glmnet
	#setseed2 sets the seeds for the bootstrap
	#iAUCtmax sets time over which risksetAUC integrates the AUC fof cox models
	#BootLab optionally estimate lambda for each bootstrap (debatable whether this is better or not)

	if(!(family=="cox" | family=="binomial")) print("No family")

	##############################################
	#--If there is penalized data:

	if(!is.null(PenData))
	{
		if(!class(PenData)=="matrix") print("make PenData a matrix")

		##############################################
		#--Size of the penalized data

		P <- dim(PenData)[2]
		N <- dim(PenData)[1]

		##############################################
		#--Include potential unpenalized data

		if(is.null(UnPenData)) PUnpen <- 0 else PUnpen <- dim(UnPenData)[2]
		X <- as.matrix(cbind(UnPenData,PenData))
		penfac <- c(rep(0,times=PUnpen),rep(1,times=P))

		##############################################
		#--This finds the best lambda. Here the average across 5 fits (for a slighlty better estimate), this could also be done once it 5 takes takes too long

		if(is.null(Lambda))
		{
			labs <- numeric()
			for(i in 1:5) 
			{
				set.seed(i+setseed1)
				labs[i] <- cv.glmnet(X,Y,family=family,standardize = FALSE, alpha = alphaIn,nfold=10,penalty.factor=penfac)$lambda.min
			}
			Lambda <- mean(labs)
		}

		##############################################
		#--Use the global lambda (Lambda) to fit the model

		model <- glmnet(X,Y,family=family,standardize = FALSE,lambda=Lambda, alpha = alphaIn,penalty.factor=penfac)
		if(family=="cox") ApparentError <- risksetAUC(Stime=Y[,1],status=Y[,2], marker=predict(model,newx=X), method="Cox", tmax=iAUCtmax,plot=FALSE)$Cindex
		if(family=="binomial") ApparentError <- as.numeric(pROC::roc(Y, as.numeric(plogis(predict(model,newx=X))))$auc)

		##############################################
		#---For saving the model, and its LP

		RealModel <- model
		LPModel <- predict(RealModel,newx=X)
		rm(model)
	}

	##############################################
	#--If there is no penalized data:

	if(is.null(PenData))
	{	
		X <- UnPenData
		N <- dim(X)[1]

		if(!class(X)=="data.frame") print("make UnpenData a dataframe")

		if(family=="cox")
		{
			model <- coxph(Y ~ .,data=X)
			ApparentError <- risksetAUC(Stime=Y[,1],status=Y[,2], marker=predict(model), method="Cox", tmax=iAUCtmax,plot=FALSE)$Cindex
		}
		if(family=="binomial")
		{
			model <- glm(Y ~ .,data=X,family="binomial")
			ApparentError <- as.numeric(pROC::roc(Y, as.numeric(plogis(predict(model))))$auc)
		}

		Lambda <- NULL
		RealModel <- model
		LPModel <- NULL

		rm(model)
	}


	##############################################
	#--Output saving for the bootstrap

	PredInBag=PredOOB <- matrix(nrow=N,ncol=NBoot)
	InBagSamples <- matrix(nrow=N,ncol=NBoot)
	AUCOOB=AUCInBag <- numeric(NBoot)+NA

	##############################################
	#--Here starts the bootstrap

	for(b in 1:NBoot)
	{
		if(b %% 100 == 0) print(b)

		##############################################
		#--Get a bootstrap sample

		set.seed(b+setseed2)
		bootsample <- sort(sample(1:N,N,replace=TRUE))

		XBoot <- X[bootsample,,drop=FALSE]
		XTest <- X[-unique(bootsample),,drop=FALSE]

		if(family=="cox") 
		{
			YTest <- Y[-unique(bootsample),]
			YBoot <- Y[bootsample,]
		}

		if(family=="binomial") 
		{
			YTest <- Y[-unique(bootsample)]
			YBoot <- Y[bootsample]
		}


		##############################################
		#--Fit a model to the bootstrapped sample

		if(!is.null(PenData))
		{
			if(BootLab==TRUE)
			{
				LambdaInBoot <- cv.glmnet(XBoot,YBoot,family=family,standardize = FALSE, alpha = alphaIn,nfold=10,penalty.factor=penfac)$lambda.min
			}
			if(BootLab==FALSE)
			{
				LambdaInBoot <- Lambda
			}

			model <- glmnet(XBoot,YBoot,family=family,standardize = FALSE,lambda=LambdaInBoot,alpha = alphaIn,penalty.factor=penfac)

			PredOOB[-unique(bootsample),b] <- predict(model,newx=XTest)
			PredInBag[,b] <- as.numeric(predict(model,newx=XBoot))
			rm(model)
		}

		if(is.null(PenData))
		{
			if(family=="cox") model <- coxph(YBoot ~ .,data=XBoot)
			if(family=="binomial") model <- glm(YBoot ~ .,data=XBoot,family="binomial")

			PredOOB[-unique(bootsample),b] <- predict(model,newdata=XTest)
			PredInBag[,b] <- as.numeric(predict(model,newdata=XBoot))
			rm(model)
		}

		##############################################
		#--Get AUC/iAUC for the bootstrap fit
		#--AUCOOB: AUC of outofbag samples
		#--AUCInBag: AUC of inbag samples 

		if(family=="cox") 
		{
			if(sum(Y[,2][-unique(bootsample)])>=3) AUCOOB[b] <- risksetAUC(Stime=Y[,1][-unique(bootsample)],status=Y[,2][-unique(bootsample)], marker=PredOOB[-unique(bootsample),b], method="Cox", tmax=iAUCtmax,plot=FALSE)$Cindex
			AUCInBag[b] <- risksetAUC(Stime=YBoot[,1],status=YBoot[,2], marker=PredInBag[,b], method="Cox", tmax=iAUCtmax,plot=FALSE)$Cindex
		}

		if(family=="binomial") 
		{
			if(sum(Y[-unique(bootsample)])>=5) AUCOOB[b] <- as.numeric(pROC::roc(Y[-unique(bootsample)], as.numeric(plogis(PredOOB[-unique(bootsample),b])))$auc)
			AUCInBag[b] <- as.numeric(pROC::roc(YBoot, as.numeric(plogis(PredInBag[,b])))$auc)
		}

		InBagSamples[,b] <- bootsample
	}


	#####################################################
	#---Calculate the AUC, AUC632, and AUC632plus


	PredictionsOOB <- rowMeans(PredOOB,na.rm=TRUE)
	if(family=="cox") AUC <- risksetAUC(Stime=Y[,1],status=Y[,2], marker=PredictionsOOB, method="Cox", tmax=iAUCtmax,plot=FALSE)$Cindex
	if(family=="binomial") AUC <- as.numeric(pROC::roc(Y, as.numeric(plogis(PredictionsOOB)))$auc)

	Fac632 <-  1-1/exp(1)
	AUC632 <- AUC*Fac632 + (1-Fac632)*ApparentError

	R <- (AUC-ApparentError)/(0.5-ApparentError)
	W <- Fac632/(1-(1-Fac632)*R)
	W <- pmin(W,1)
	AUC632plus <- W*AUC + (1-W)*ApparentError

	Means <- c(AUC,AUC632,AUC632plus,ApparentError)
	names(Means) <- c("BootstrapOOB","Bootstrap632","Bootstrap632plus","ApparentError")


	#####################################################
	#---Determine size of the CI's and place these around the AUC, AUC632, and AUC632plus
	#---The idea is that the inbag samples give the best indication of the size of the confidence intervals (but the are centered around the apparant error)
	#---If the apparent error's AUC is >0.90 high this (probably) does not work
	#---The size of the these confidence intervals is then placed around the OOB estimate


	mAUCInBag <- mean(AUCInBag,na.rm=TRUE)
	CIBase <- boot:::norm.inter(AUCInBag, c(0.025,0.975))[,2]-mAUCInBag

	CIs <- matrix(nrow=3,ncol=2)
	CIs[1,] <- AUC+CIBase
	CIs[2,] <- AUC632+CIBase
	CIs[3,] <- AUC632plus+CIBase
	colnames(CIs) <- c("Lower","Upper")
	rownames(CIs) <- c("BootstrapOOB","Bootstrap632","Bootstrap632plus")


	#####################################################
	#---Output

	if(outBig==0) out <- list(Lambda=Lambda,Means=Means,CIs=CIs)
	if(outBig==1) out <- list(EffectiveBoot=sum(!is.na(AUCOOB)),InBagSamples=InBagSamples,PredictionsOOB=PredictionsOOB,PredInBag=PredInBag,PredOOB=PredOOB,Lambda=Lambda,CIBase=CIBase,AUCInBag=AUCInBag,RealModel=RealModel,LPModel=LPModel,Means=Means,CIs=CIs)
	return(out)
}

#####################################################
#---end of function

