# aiming at elucidating the predictive power
# comparing one measurement and linear prediction

#######
### 1. Build the best model using a break point at 1-1
### 2. Model evalulation using mean and sd
### 3. Compare the model with the best model using quadratic
### 3. Prove that adding more data is not necessary
### 4. Quality control of the model


#### Inputs ####
# 1 When performing prediction using only fitness data
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Data2 <- read.table(file="ExpDataAllConditions100.txt", header=T, sep="\t", stringsAsFactors = F);
DataCore <- Data2[Data2$T0All>=500,]

require(MASS)
setwd("C:/Users/lichuan/Dropbox/1.3 Mutiple conditions/7Analysis/Plots/7PredictFitness")

type <- "NonZero"; # Choose from "N1", "N1N2", "N1N2NoBound", "NonZero"

nAll <- 37:40; # measurement for average fitness
nColumns <- list(21:25, 26:30, 31:33, 34:36)# measurement for individual columns for fitness

getModIndex <- function(robustModel1){
  k1 <- summary(robustModel1)
  k1CI <- confint.default(robustModel1) 
  # k1CI[1]
  # k1CI[2]
  return(c(k1$coefficients[1], k1$coefficients[2], k1CI[1], k1CI[2]))
}

orgModPerf <- function(returnList){
  k1 <- returnList[[1]]
  k2 <- returnList[[2]]
  DirFit <- returnList[[3]]
  CorPear <- cor(DirFit$Y, DirFit$predY)
  CorSpea <- cor(DirFit$Y, DirFit$predY, method = "spearman")
  return(c(k1, k2, CorPear, CorSpea))
}

yRange = 0.1
TableS2 <- c()

### S1 PredFitDir - use all available to make a good model
PredFitDir <- function(Train, Test){
  # Train = Data[Data$Num==1, c(nY, nX)]
  # Test = Data[Data$Num<3, c(nY, nX)]
  colnames(Train) <- c("Y", "X")
  colnames(Test) <- c("Y", "X")
  
  ### record the order of Train and Test
  # orderTrain <- 1:(dim(Train)[1])
  orderTest <- 1:(dim(Test)[1])
  
  ### Adjust values of X and Y
  Train <- Train - 1
  Test <- Test - 1
  
  ###### Part1 predict deleterious mutations ###
  ### Use only deleterious mutations
  Train1 <- Train[Train$X<0,] # use only deleterious mutation at 30C to predict the model
  # orderTrain1 <- orderTrain[Train$X<0]
  Test1 <- Test[Test$X<0,] # to predict the fitness at 23C
  orderTest1 <- as.numeric(orderTest[Test$X<0])
  
  ### rlm force to go through origin
  robustModel1 = rlm(Y ~ X - 1, data = Train1)
  
  TestX1 <- data.frame(X=Test1[,2])
  TestYPred1 <- predict(robustModel1, TestX1)
  
  ###### Part2 predict beneficial mutations ###
  ### Use only beneficial mutations
  Train2 <- Train[Train$X>0,] # use only deleterious mutation at 30C to predict the model
  # orderTrain2 <- orderTrain[Train$X>0]
  Test2 <- Test[Test$X>0,] # to predict the fitness at 23C
  orderTest2 <- as.numeric(orderTest[Test$X>0])
  
  ### rlm force to go through origin
  robustModel2 = rlm(Y ~ X - 1, data = Train2)
  # summary(robustModel)
  
  TestX2 <- data.frame(X=Test2[,2])
  TestYPred2 <- predict(robustModel2, TestX2)
  
  ###### Part 3 add in the case where fitness = 1 for the test set
  Test3 <- Test[Test$X==0,]
  TestYPred3 <- rep(0,sum(Test$X==0))
  orderTest3 <- as.numeric(orderTest[Test$X==0])
  
  ### combine the two sets
  Test <- rbind(Test1, Test2, Test3)
  TestYPred <- as.vector(c(TestYPred1, TestYPred2, TestYPred3, use.names=F))
  orderTest <- c(orderTest1, orderTest2, orderTest3)
  
  ### Add back the adjusted values
  Test <- Test + 1
  TestYPred <- TestYPred + 1
  
  TestYPred[TestYPred < 0.5] <- 0.5
  
  ### reorder test
  Result <- cbind(Test, TestYPred)
  Result <- Result[order(orderTest),]
  colnames(Result) <- c("Y", "X", "predY")
  # get the point when X=0 and X=2
  
  # https://www.quora.com/How-do-you-find-the-confidence-intervals-for-robust-regression-Tukey-Least-Trimmed-Square-M-estimation-in-R
  k1 <- getModIndex(robustModel1)
  
  k2 <- getModIndex(robustModel2)
  
  # return multiple values
  # http://stackoverflow.com/questions/1826519/how-to-assign-from-a-function-which-returns-more-than-one-value
  return(list(round(k1,4), round(k2,4), Result))
}

for (nXEnv in 1:4){
  for (nYEnv in 1:4){
    if (nXEnv == nYEnv){
      next;
    }    
    
    #### A set of variable to choose from ####
    nX <- nAll[nXEnv]; # 30C as predictor
    nY <- nAll[nYEnv]; # 23C to be predicted
    
    Xcolumns <- nColumns[[nXEnv]]
    Ycolumns <- nColumns[[nYEnv]]
    
    myBool <- apply(DataCore[,c(Xcolumns,Ycolumns)],1, min) > 0.5
    
    # sum(myBool) # 6319
    
    myBool[1] <- F
    Data <- DataCore[myBool,] # Data used to build the model
    
    FileNameAll <-  c("23","30","DM","37")
    XEnv <- FileNameAll[nXEnv]
    YEnv <- FileNameAll[nYEnv]
    FileAdd <- paste0(FileNameAll[nXEnv], "_",FileNameAll[nYEnv])
    
    # load packages
    
    #### Build the model ####
    if (type=="N1"){
      returnList <- PredFitDir(Train = Data[Data$Num<2, c(nY, nX)], Test = Data[Data$Num>-1, c(nY, nX)])
    }else if(type == "N1N2"){
      returnList <- PredFitDir(Train = Data[Data$Num<3, c(nY, nX)], Test = Data[Data$Num>-1, c(nY, nX)])
    }else if(type == "NonZero"){
      returnList <- PredFitDir(Train = Data[, c(nY, nX)], Test = Data[, c(nY, nX)]) # 
    }else if (type == "Mixture"){
      returnList <- PredFitDir(Train = Data[, c(nY, nX)], Test = Data2[Data2$T0All>0, c(nY, nX)])  
    }
    
    
    DirFit <- returnList[[3]]
    write.table(x = DirFit, file=paste0("FitnessObsAndPred_",FileAdd,".txt"),  
                quote=FALSE, sep="\t",row.names = T);
    
    #### 2.2 Calculate model and variance for each point ####
    ### plot out the prediction and actual values, mark different number of mutant with different colors and plot N1 the last 
    FitRange <- c(0.5,1.4)
    
    ### Calculate the deviation of measurement
    getMeasureDev <- function(nX, Xcolumns){
      X30 <- Data[,nX] # average
      Se30 <- rowMeans(abs(Data[,Xcolumns]-X30))
      return(Se30)
    }
    
    ### Calculate bias and measurement bias
    BiasDir <- DirFit$predY - DirFit$Y
    
    Se30 <- getMeasureDev(nX, Xcolumns)
    Se23 <- getMeasureDev(nY, Ycolumns)
    
    BiasDevOne <- function(ResultAll_All=DirFit, 
                           Biasfilename="FigureBias.pdf", 
                           Varfilename ="FigureError.pdf"){
      X <- ResultAll_All$X
      Bias <- ResultAll_All$predY - ResultAll_All$Y
      
      setwd("C:/Users/lichuan/Dropbox/1.3 Mutiple conditions/7Analysis/Plots/7PredictFitnessFinal/")
      pdf(Biasfilename, width = 5, height = 5)
      plot(x=1, y=0, xlim = FitRange, ylim = c(-yRange, yRange), cex=0, xlab = paste0("Fitness observed at ", XEnv ,"°C"), ylab = "Bias in model prediction")
      lines(c(-1,3),c(0,0),lty=2, lwd=2)
      lines(lowess(Bias ~ X, f = .2), col = "red", lwd=2)
      
      dev.off()
      
      ## 2.2.2 Variance ##
      # Plot error against X
      setwd("C:/Users/lichuan/Dropbox/1.3 Mutiple conditions/7Analysis/Plots/7PredictFitnessFinal/")
      Var <- abs(Bias)
      
      pdf(Varfilename, width = 5, height = 5)
      plot(x=1, y=0,xlim = FitRange, ylim = c(0, yRange), cex=0, xlab = paste0("Fitness observed at ", XEnv ,"°C"), ylab = "Deviation in model prediction")
      
      lines(lowess(Var ~ X, f = .2), col = "red", lwd=2)
      lines(lowess(Se23 ~ X, f = .2), col = "blue", lwd=2, lty=5)
      # lines(lowess(Se30 ~ X, f = .2), col = "gray", lwd=2, lty=5)
      
      legend("topright", 
             c("Prediction error", paste0("Measurement error")),
             lty = c(1,5), lwd=2, col = c("red","blue"))
      
      dev.off()
    }
    BiasDevOne(DirFit, Biasfilename=paste0("FigureBias_", FileAdd, ".pdf"), 
               Varfilename=paste0("FigureError_", FileAdd, ".pdf"))
    
    getTwoCor <- function(){
      Result <- c()
      
      # correlation
      Result[1] <- cor(DirFit$Y, DirFit$predY)
      
      # Rank correlation
      Result[2] <- cor(DirFit$Y, DirFit$predY, method = "spearman")
      
      # mean bias
      Result[3] <- mean(BiasDir)
      
      # mean dev
      Result[4] <- mean(abs(BiasDir))
      
      # median bias
      Result[5] <- median(BiasDir)
      
      # median dev
      Result[6] <- median(abs(BiasDir))
      
      return(Result)
    }
    
    Result <- c(XEnv, YEnv, sum(myBool), returnList[[1]][1:2], returnList[[2]][1:2])
    Result <- c(Result, round(c(getTwoCor(), mean(Se30), median(Se30), mean(Se23), median(Se23)),4))
    TableS2 <- rbind(TableS2, Result)
  }
}

colnames(TableS2) <- c("EnvX",	"EnvY",	"n", "k1", "k1_sd", "k2", "k2_sd", 
                   "LinCor","RankCor","MeanBias","MeanDev","MedianBias","MedianDev",
                   "meanErrX", "medianErrX", "meanErrY", "medianErrY")

setwd("C:/Users/lichuan/Dropbox/1.3 Mutiple conditions/7Analysis/Plots/7PredictFitnessFinal/")

write.table(x = TableS2, file="TableS2",  quote=FALSE, sep="\t",row.names = T);


