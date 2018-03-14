rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

Data <- read.table(file="ExpDataAllConditions100.txt", header=T, sep="\t", stringsAsFactors = F);

#### 4. Calculates epistasis for double mutants ####
### S4: calcEpi - Calculate epistasis for a column of fitness values
calcEpi <- function(Data, ColFit){ # ColFit means the column showing fitness value
  Fitness1 <- Data[Data$Num==1,]
  Fitness2 <- Data[Data$Num==2,]
  
  ExpFit <- c();
  
  for (i in 1:length(Fitness2[,1])){
    ObsFit <- Fitness2[i,ColFit];
    MutPos <- unlist(strsplit(x=Fitness2$Pos[i],split = " "));
    MutNuc <- unlist(strsplit(x=Fitness2$Nuc[i],split = " "));
    bool1 <- (Fitness1$Pos == as.numeric(MutPos[1]))&(Fitness1$Nuc == MutNuc[1]);
    bool2 <- (Fitness1$Pos == as.numeric(MutPos[2]))&(Fitness1$Nuc == MutNuc[2]);
    
    Fit1 <- Fitness1[bool1, ColFit]
    Fit2 <- Fitness1[bool2, ColFit]
    ExpFit <- c(ExpFit, max(0.5, Fit1*Fit2)) # apply cutoff of 0.5
  }
  Epistasis <- Fitness2[,ColFit] - ExpFit;
  return(Epistasis)
}

### S5: calcEpiExpfit - calculate both epistasis and expected fitness
calcEpiExpfit <- function(Data, ColFit){ # ColFit means the column showing fitness value
  Fitness1 <- Data[Data$Num==1,]
  Fitness2 <- Data[Data$Num==2,]
  
  ExpFit <- c();
  
  for (i in 1:length(Fitness2[,1])){
    ObsFit <- Fitness2[i,ColFit];
    MutPos <- unlist(strsplit(x=Fitness2$Pos[i],split = " "));
    MutNuc <- unlist(strsplit(x=Fitness2$Nuc[i],split = " "));
    bool1 <- (Fitness1$Pos == as.numeric(MutPos[1]))&(Fitness1$Nuc == MutNuc[1]);
    bool2 <- (Fitness1$Pos == as.numeric(MutPos[2]))&(Fitness1$Nuc == MutNuc[2]);
    
    Fit1 <- Fitness1[bool1, ColFit]
    Fit2 <- Fitness1[bool2, ColFit]
    ExpFit <- c(ExpFit, max(0.5, Fit1*Fit2))
  }
  Epistasis <- Fitness2[,ColFit] - ExpFit;
  return(cbind(Epistasis,ExpFit))
}

### 4.1 Calculate epistasis at all four conditions
# When we focus on all N2 mutants with over 100 reads
head(Data[,37:40]) 
Epist23 <- calcEpi(Data, 37)
Epist30 <- calcEpi(Data, 38)
EpistDM <- calcEpi(Data, 39)
Epist37 <- calcEpi(Data, 40)

EpistAll <- cbind(Epist23, Epist30, EpistDM, Epist37)
CorAll <- matrix(nrow = 4, ncol = 4)
row.names(CorAll) <- c("23","30","DM","37")
colnames(CorAll) <- c("23","30","DM","37")

for (i in 1:3){
  for (j in (i+1):4){
    if (i != j){
      k <- cor.test(EpistAll[,i], EpistAll[,j], method = "pearson")
      CorAll[i,j] <- k$estimate
      k <- cor.test(EpistAll[,i], EpistAll[,j], method = "spearman")
      CorAll[j,i] <- k$estimate
    }
  }
}

range(CorAll[is.finite(CorAll)])


### 4.2 Plot epistasis
emin <- -0.7
emax <- 0.7
plotdist <- function(diff,xlab="Epistasis",bin=30){
  par(mar=c(5,6,1,1));
  k<-hist(diff,bin, main=NULL, xlab = xlab, xlim = c(emin, emax),
          col="grey",border="white",cex.lab=2, cex.axis=1.5);
  lines(c(0,0),c(0,max(k$counts)*1.02),col="red",lwd=3)
}


pdf(file="Epistasis23.pdf", width = 8,height = 5)
plotdist(EpiData$Epist23)
dev.off()

pdf(file="Epistasis30.pdf", width = 8,height = 5)
plotdist(Epist30)
dev.off()

pdf(file="EpistasisDMSO.pdf", width = 8,height = 5)
plotdist(EpistDM)
dev.off()

pdf(file="Epistasis37.pdf", width = 8,height = 5)
plotdist(Epist37)
dev.off()


### S6: CalcMeanSd - Calculate and get printable versions of mean and sd
calcMeanSd <- function(Fitness){
  meanFit <- mean(Fitness)
  sdFit <- sd(Fitness)
  str <- paste(c(round(meanFit,3), "+/-" , round(sdFit,3)),sep = "", collapse = "") 
  return(str)
}

calcMeanSd(Epist23)
calcMeanSd(Epist30)
calcMeanSd(EpistDM)
calcMeanSd(Epist37)

### S7: compare0 - Calculate the fraction of positive and negative values
compare0 <- function(Epi){
  Pos <- sum(Epi>0)
  Neg <- sum(Epi<0)
  return(c(Pos/length(Epi), Neg/length(Epi)))
}

compare0(EpiData$Epist23)
compare0(EpiData$Epist30)
compare0(EpiData$EpistDM)
compare0(EpiData$Epist37)

#### 5. Compare epistasis in different conditions ####
require(MASS) 
### S8: drawcorEpi - draw the epistasis vs epistasis
drawcorEpi <- function(figname, column1, column2, xlabtext, ylabtext){ # 30 first, 37 second
  pdf(figname, width=4, height=4) 
  cor1 <- cor.test(column1, column2) 
  cor1 <- cor1$estimate
  cor2 <- cor.test(column1, column2, method = "spearman") 
  cor2 <- cor2$estimate
  
  plot(column1, column2, pch=19, cex=0.2, 
       col=rgb(0,0,0,0.2), 
       xlab = xlabtext, 
       ylab = ylabtext, 
       xlim = c(emin, emax), ylim = c(emin, emax))
  
  lines(x=c(-2,1.4),y=c(-2,1.4),col="red",lwd=1)
  text(-0.45, 0.55, labels = paste0("r=", round(cor1,3)))
  text(-0.45, 0.4, labels = paste0("o=", round(cor2,3)))
  z <- rlm(column2 ~ column1)
  abline(z,col="blue",lty="dashed", lwd=2)
  dev.off()
}

### 5.1 plot all pair-wise correlation of epistasis 
drawcorEpi("EpiCorrelation23_30Robust.pdf", Epist23, Epist30, "Epistasis at 23C", "Epistasis at 30C")
drawcorEpi("EpiCorrelation23_DMSORobust2.pdf", Epist23, EpistDM, "Epistasis at 23C", "Epistasis at DMSO")

### 5.4 print out epistasis result
EpiData <- cbind(Data[Data$Num==2,], Epist23 = Epist23, Epist30 = Epist30, EpistDM = EpistDM, Epist37 = Epist37)
write.table(x=EpiData,file = "ExpDataDouble.txt",quote=FALSE, sep="\t");

Fitness1 <- Data[Data$Num == 1,]
### S10: getPforEpistasis - calculate the p-value for epistasis
getPforEpistasis <- function(EpiData, Fitness1, colNum){
  pvalue <- c()
  for (i in 1:dim(EpiData)[1]){
    ObsFit <- EpiData[i,colNum];
    MutPos <- unlist(strsplit(x=EpiData$Pos[i],split = " "));
    MutNuc <- unlist(strsplit(x=EpiData$Nuc[i],split = " "));
    bool1 <- (Fitness1$Pos == as.numeric(MutPos[1])) & (Fitness1$Nuc == MutNuc[1]);
    bool2 <- (Fitness1$Pos == as.numeric(MutPos[2])) & (Fitness1$Nuc == MutNuc[2]);
    
    Fit1 <- Fitness1[bool1,colNum]
    Fit2 <- Fitness1[bool2,colNum]
    ExpAdd <- Fit1*Fit2
    ExpAdd[ExpAdd < 0.5] <- 0.5
    
    Epistasis <- ObsFit-ExpAdd;
    if (length(unique(as.numeric(Epistasis)))==1){ # if all values are the same
      if (unique(as.numeric(Epistasis))==0){
        pvalue <- c(pvalue, 1) # If all values are 0, not significantly different from 0
      }else{
        pvalue <- c(pvalue, 0) # if all values are non-0 constant, will be significant
      }
    }else{
      k <- t.test(Epistasis) 
      pvalue <- c(pvalue, k$p.value)
    }
  }
  return(pvalue)
}

### 6.1 calculate p value for epistasis in each condition
colNum <- 21:25 # one condition at a time
head(Data[colNum])

pvalEpi23 <- getPforEpistasis(EpiData, Fitness1, 21:25)
pvalEpi30 <- getPforEpistasis(EpiData, Fitness1, 26:30)
pvalEpiDM <- getPforEpistasis(EpiData, Fitness1, 31:33)
pvalEpi37 <- getPforEpistasis(EpiData, Fitness1, 34:36)

EpiData <- cbind(EpiData, pvalEpi23, pvalEpi30, pvalEpiDM, pvalEpi37)


### S11: plotEpi - Plot individually significant and nonsignificant ones
plotEpi <- function(Epistasis, pvalue, yrange){
  par(mar=c(4,4.1,1,1));
  par(oma = c(3, 4, 0, 0));
  
  Mybreak <- seq(from=emin, to=emax, by=0.05);
  Xrange <- c(emin,emax)
  Yrange <- yrange
  
  p2 <- hist(Epistasis, plot=F,breaks =Mybreak);
  plot(p2, col=rgb(0,0,0,1/8), xlim=Xrange, ylim=Yrange, xlab = "Epistasis", cex.lab=1.4,
       main=NULL,freq = T, border="gray", ylab="Number of pairs of mutations")
  
  FitSig2 <- Epistasis[pvalue < 0.05];
  p1 <- hist(FitSig2, plot=F, breaks =Mybreak);
  plot(p1, col=rgb(0,0,1,1/2), breaks =Mybreak, xlim=Xrange, add=T, lty="blank") 
  lines(x=c(0,0),y=c(0,max(p2$counts)*1.05),col="red",lwd=2)
  return(1)
}

# plot significant epistasis
plotEpi(EpiData$Epist23, EpiData$pvalEpi23, c(0,1500))

plotEpi(EpiData$Epist30, EpiData$pvalEpi30, c(0,1500))
EpiData <- EpiData[EpiData$T0All > 500,]
write.table(EpiData, file = "EpiDataAll.txt",quote=FALSE, sep="\t",row.names = T);

### Count the fraction that change signs across environment
EpiData <- read.table(file="EpiDataAll.txt", header=T, sep="\t", stringsAsFactors = F);

Env1 <- 46
Env2 <- 47

Result <- matrix(nrow=2, ncol=2)
rownames(Result) <- c("PosEnv2","NegEnv2")
colnames(Result) <- c("PosEnv1","NegEnv1")

myBool <- EpiData[Env1] > 0 & EpiData[Env1+4] < 0.05 & EpiData[Env2] > 0 & EpiData[Env2+4] < 0.05
Result[1,1] <- sum(myBool)
myBool <- EpiData[Env1] < 0 & EpiData[Env1+4] < 0.05 & EpiData[Env2] < 0 & EpiData[Env2+4] < 0.05
Result[2,2] <- sum(myBool)
myBool <- EpiData[Env1] < 0 & EpiData[Env1+4] < 0.05 & EpiData[Env2] > 0 & EpiData[Env2+4] < 0.05
Result[1,2] <- sum(myBool) # 0
myBool <- EpiData[Env1] > 0 & EpiData[Env1+4] < 0.05 & EpiData[Env2] < 0 & EpiData[Env2+4] < 0.05
Result[2,1] <- sum(myBool) # 7

# calculate the fraction for other environment pairs as well

