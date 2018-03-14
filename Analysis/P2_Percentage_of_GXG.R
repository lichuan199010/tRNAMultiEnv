rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Data <- read.table(file="ExpDataAllConditions100.txt", header=T, sep="\t", stringsAsFactors = F);

### GxE by percentage ####
head(Data[37:40])
colAll <- 37:40

# three different ways to calculate percentages
Result1 <- matrix(nrow = 4, ncol = 4)
row.names(Result1) <- colnames(Data[37:40])
colnames(Result1) <- colnames(Data[37:40])

Result2 <- matrix(nrow = 4, ncol = 4)
row.names(Result2) <- colnames(Data[37:40])
colnames(Result2) <- colnames(Data[37:40])


for (i in 1:3){
  for (j in (i+1):4){
    # Method1: proportion by comparing which one is bigger
    MyBool1 <- (Data[,colAll[i]] < Data[,colAll[j]]) & ((Data[,colAll[i]]-1)*(Data[,colAll[j]]-1)>0)
    MyBool2 <- (Data[,colAll[i]] > Data[,colAll[j]]) & ((Data[,colAll[i]]-1)*(Data[,colAll[j]]-1)>0)
    Result1[i,j] <- sum(MyBool1)/(sum(MyBool1)+sum(MyBool2))
    
    # Method2: Larger fitness changes (Use 1/8 of the plot for magnitude change)
    MyBool1 <- (abs(Data[,colAll[i]]-1) - abs(Data[,colAll[j]]-1)>0) & ((Data[,colAll[i]]-1)*(Data[,colAll[j]]-1)>0)
    MyBool2 <- (abs(Data[,colAll[i]]-1) - abs(Data[,colAll[j]]-1)<0) & ((Data[,colAll[i]]-1)*(Data[,colAll[j]]-1)>0)
    Result2[i,j] <- sum(MyBool1)/(sum(MyBool1)+sum(MyBool2))
  }
}


### GxE from direct comparison ####
Result3 <- matrix(nrow = 4, ncol = 4) # t.test
row.names(Result3) <- colnames(Data[37:40])
colnames(Result3) <- colnames(Data[37:40])
head(Data[21:36])

colAll <- list(21:25, 26:30, 31:33, 34:36)

CalcTtestP <- function(Data, columnX=colAll[[1]], columnY=colAll[[3]]){
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
    Allp <- CalcTtestP(Data, colAll[[i]], colAll[[j]])
    Result3[i,j] <- sum(Allp<0.05)/length(Allp)
  }
}


### Do a color plot ####
PlotMat <- Result3

for (i in 1:3){
  for (j in (i+1):4){
    PlotMat[5-i, 5-j] <- abs(2*Result2[i,j]-1)
  }
}
t(round(PlotMat, 3)*100)

xaxislab = c("23C","30C","DMSO","37C")

Mycol <- rev(rainbow(n=32, start=0,end=0.2))

pdf("PercentSignificantByTtest.pdf",width = 7,height = 7)
image(PlotMat, col=Mycol,xlab = "",ylab = "", xaxt='n', yaxt='n',zlim = c(0,1))
xaxis = c(0,1/3,2/3,1)
mtext(side=3, at=xaxis, text=xaxislab, cex=1.0, line=0.2, col=1)
mtext(side=2, at=xaxis, text=xaxislab, cex=1.0, line=0.2, col=1)
legend(grconvertX(0.5, "device"), grconvertY(1, "device"),
       c("0",".5","1"), fill=Mycol)
dev.off()
