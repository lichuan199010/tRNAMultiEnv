### percentage of GxGxE ##### 
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Data <- read.table(file="EpiDataAll.txt", header=T, sep="\t", stringsAsFactors = F);

head(Data[46:49])
colAll <- 46:49

Result1 <- matrix(nrow = 4, ncol = 4)
row.names(Result1) <- colnames(Data[colAll])
colnames(Result1) <- colnames(Data[colAll])

Result2 <- matrix(nrow = 4, ncol = 4)
row.names(Result2) <- colnames(Data[colAll])
colnames(Result2) <- colnames(Data[colAll])

Result3 <- matrix(nrow = 4, ncol = 4)
row.names(Result3) <- colnames(Data[colAll])
colnames(Result3) <- colnames(Data[colAll])

for (i in 1:3){
  for (j in (i+1):4){
    # Method1: proportion by comparing which one is bigger
    MyBool1 <- (Data[,colAll[i]] < Data[,colAll[j]]) & (Data[,colAll[i]]*Data[,colAll[j]]>0)
    MyBool2 <- (Data[,colAll[i]] > Data[,colAll[j]]) & (Data[,colAll[i]]*Data[,colAll[j]]>0)
    Result1[i,j] <- sum(MyBool1)/(sum(MyBool1)+sum(MyBool2))
    
    # Method2: Larger fitness changes (Use 1/8 of the plot for magnitude change)
    MyBool1 <- (abs(Data[,colAll[i]]) - abs(Data[,colAll[j]])>0) & (Data[,colAll[i]]*Data[,colAll[j]]>0)
    MyBool2 <- (abs(Data[,colAll[i]]) - abs(Data[,colAll[j]])<0) & (Data[,colAll[i]]*Data[,colAll[j]]>0)
    Result2[i,j] <- sum(MyBool1)/(sum(MyBool1)+sum(MyBool2))
  }
}


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

colFitInd <- 21:36
EpiAll <- c()
for (i in colFitInd){
  EpiAll <- cbind(EpiAll, calcEpi(Data, i))
}

colAll <- list(1:5, 6:10, 11:13, 14:16)

CalcTtestP <- function(Data, columnX, columnY){
  pVal <- c()
  
  for (i in 1:dim(Data)[1]){
    Set1 <- as.numeric(Data[i, columnX])
    Set2 <- as.numeric(Data[i, columnY])
    if (length(unique(Set1)) == 1 & length(unique(Set2)) == 1){
      if (Set1[1] == Set2[1]){
        pVal <- c(pVal, 1)
      }else{
        pVal <- c(pVal, 0)
      }
      
    }else{
      k <- t.test(Data[i, columnY], Data[i, columnX])
      pVal <- c(pVal, k$p.value)
    }
  }
  return(pVal)
}

for (i in 1:3){
  for (j in (i+1):4){
    Allp <- CalcTtestP(EpiAll, colAll[[i]], colAll[[j]])
    Result3[i,j] <- sum(Allp<0.05)/length(Allp)
  }
}

setwd("C:/Users/lichuan/Dropbox/1.3 Mutiple conditions/7Analysis/Plots/2GxEandEpixE")
write.table(x=Result1, file = "EpixE_Result1.txt",quote=FALSE, row.names = F, sep="\t");
write.table(x=Result2, file = "EpixE_Result2.txt",quote=FALSE, row.names = F, sep="\t");
write.table(x=Result3, file = "EpixE_Result3.txt",quote=FALSE, row.names = F, sep="\t");


### Do a color plot ####
PlotMat <- Result3
for (i in 1:3){
  for (j in (i+1):4){
    PlotMat[5-i, 5-j] <- abs(2*Result2[i,j]-1)
  }
}

Mycol <- rev(rainbow(n=32, start=0,end=0.2))

pdf("Percentage of Epistasis.pdf",width = 7,height = 7)
image(PlotMat, col=Mycol,xlab = "",ylab = "", xaxt='n', yaxt='n',zlim = c(0,1))
xaxis = c(0,1/3,2/3,1)
xaxislab = c("23C","30C","DMSO","37C")
mtext(side=3, at=xaxis, text=xaxislab, cex=2.0, line=0.2, col=1)
mtext(side=2, at=xaxis, text=xaxislab, cex=2.0, line=0.2, col=1)
legend(grconvertX(0.5, "device"), grconvertY(1, "device"),
       c("0",".5","1"), fill=Mycol)
dev.off()

# order on the plot
t(round(PlotMat, 3)*100)
