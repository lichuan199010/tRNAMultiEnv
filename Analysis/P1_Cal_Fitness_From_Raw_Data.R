### 1. Raw data compilation ###
### 2. Calculate # of generations ###
### 3. Remove unwanted variants ###
### 4. Calculate Fitness ###
### 5. Calculate Fitness ###
### 7. Add cutoff to fitness values ###

#=== Main Program ===#
#### 1. Raw data compilation ####
# read in raw data and extract useful information
rm(list=ls())
# 
# setwd("E:/Dropbox/1.3 Mutiple conditions/7Analysis/Outputs"); # at home
setwd("D:/Outputs");

# read in all data
wt <- read.table(file="Mut0.txt", header=F, sep="\t", stringsAsFactors = F);
mut1 <- read.table(file="Mut1.txt", header=F, sep="\t", stringsAsFactors = F);
mut2 <- read.table(file="Mut2.txt", header=F, sep="\t", stringsAsFactors = F);
mut3 <- read.table(file="Mut3.txt", header=F, sep="\t", stringsAsFactors = F);
mut4 <- read.table(file="Mut4.txt", header=F, sep="\t", stringsAsFactors = F);
mut5 <- read.table(file="Mut5.txt", header=F, sep="\t", stringsAsFactors = F);
mut6 <- read.table(file="Mut6.txt", header=F, sep="\t", stringsAsFactors = F);
mut7 <- read.table(file="Mut7.txt", header=F, sep="\t", stringsAsFactors = F);
mut8 <- read.table(file="Mut8.txt", header=F, sep="\t", stringsAsFactors = F);
mut8Plus <- read.table(file="Mut9.txt", header=F, sep="\t", stringsAsFactors = F);


groupMut <- c(0,rep(1,nrow(mut1)),rep(2,nrow(mut2)),rep(3,nrow(mut3)),
              rep(4,nrow(mut4)),rep(5,nrow(mut5)),rep(6,nrow(mut6)),
              rep(7,nrow(mut7)),rep(8,nrow(mut8)),rep(9,nrow(mut8Plus)));

length(groupMut) 

Data <- cbind(groupMut, rbind(wt, mut1, mut2, mut3, mut4, mut5, mut6, mut7, mut8, mut8Plus));
head(Data)

names(Data) <- c("Num", "Pos","Nuc", "T0All","T30C1","T30C2","T30C3",
                 "T37C1","T37C2","T37C3","T0R1", # Lane 1 above
                 "D1", "D2", "D3", 
                 "T23C1", "T23C2", "T23C3", "T0R2", "Ctrl",# Lane 2 above
                 "T0R1Add","T30C5","T30C4","T37C1Add", "T37C3Add", # Lane 3 above
                 "T0R2Add","T23C5","T23C4",
                 "T37C2Add","T30C1Add","T30C2Add","T23C3Add"); # Lane 4 above

# 2. Combine the datasets together
# 5 rep for 23 and 30C, 3 rep for DMSO and 37C
Data2 <- Data[,1:20]
names(Data2) <- c("Num", "Pos","Nuc", "T0All", # 1-4
                  "T23C1", "T23C2", "T23C3", "T23C4", "T23C5", # 5-9
                  "T30C1", "T30C2", "T30C3", "T30C4", "T30C5", # 10-14
                  "D1", "D2", "D3",  # 15-17
                  "T37C1","T37C2","T37C3") # 18-20
Data2[,5:9] <- cbind(Data$T23C1, Data$T23C2, Data$T23C3+Data$T23C3Add, Data$T23C4, Data$T23C5)
Data2[,10:14] <- cbind(Data$T30C1+Data$T30C1Add, Data$T30C2+Data$T30C2Add, Data$T30C3, Data$T30C4, Data$T30C5)
Data2[,15:17] <- cbind(Data$D1, Data$D2, Data$D3)
Data2[,18:20] <- cbind(Data$T37C1+Data$T37C1Add, Data$T37C2+Data$T37C2Add, Data$T37C3+Data$T37C3Add)

head(Data2[4,])

write.table(x=Data2,file = "CountDataAllConditions.txt",quote=FALSE, sep="\t",row.names = F);

Data <- Data2

#### 2. Calculate # of generations ####
### 2.1 previously calculated generation numbers from OC changes
G23 <- 13.666
G30 <- 11.305
GDM <- 11.894
G37 <- 11.506

### 2.2 If we look at the porportion of t0 among all the cells with the correct length of insertion
# t0
t0All <- sum(Data$T0All)
t0Wt <- max(Data$T0All)

### 2.3 add in extra generation for wt genotype
# for each condition calculate the total number of reads
head(Data[,5:9])
t1All23 <- sum(sum(Data[,5:9]))
t1Wt23 <-  sum(Data[1,5:9])
ngenadd23 <- log2(t1Wt23/t0Wt*t0All/t1All23)

head(Data[,10:14])
t1All30 <- sum(sum(Data[,10:14]))
t1Wt30 <- sum(Data[1,10:14])
ngenadd30 <- log2(t1Wt30/t0Wt*t0All/t1All30)

head(Data[,15:17])
t1AllDM <- sum(sum(Data[,15:17]))
t1WtDM <- sum(Data[1,15:17])
ngenaddDM <- log2(t1WtDM/t0Wt*t0All/t1AllDM)

head(Data[,18:20])
t1All37 <- sum(sum(Data[,18:20]))
t1Wt37 <- sum(Data[1,18:20])
ngenadd37 <- log2(t1Wt37/t0Wt*t0All/t1All37)

### 2.4 add in these generations after correcting for the fraction of WT cells
ngen23 <- G23+ngenadd23
ngen30 <- G30+ngenadd30
ngenDM <- GDM+ngenaddDM
ngen37 <- G37+ngenadd37


#### 3. Remove unwanted variants ####
# with mutation at invariable sites CG...G
unexpPos <- c();

for (i in 1:dim(Data)[1]){
  k <- as.numeric(unlist(strsplit(x=toString(Data$Pos[i]),split = " ")))
  if(sum(k>69)|sum(k<1)){
    unexpPos  <- c(unexpPos, i)   
  }
}

Data2 <- Data[-unexpPos,] # not including wt

# add back wt row to the front
Data2 <- rbind(Data[1,], Data2)

#### 4. Calculate Fitness ####
### 4.1 Calculate ratio
colT1 <- c(5:20) # row with fitness values at the 4 conditions
wt <- Data2[1,]

RatioEach <-Data2[,colT1]/Data2[,4]
RatioSave <- RatioEach

for (i in 1:dim(RatioEach)[2]){
  RatioEach[,i] <- RatioEach[,i]/as.numeric(wt[colT1][i]/wt[4]);
}


#### 5. Calculate Fitness ####
### 5.1 Calculate fitness for all
ngenAll <- c(rep(ngen23, 5), rep(ngen30,5),rep(ngenDM, 3), rep(ngen37,3))

RatioEach2 <- RatioEach;
FitnessEach <- RatioEach2;

for (i in 1:dim(RatioEach2)[2]){
  FitnessEach[,i] <- RatioEach2[,i]^(1/ngenAll[i])
}

dim(FitnessEach)
head(FitnessEach)

FitnessEach[FitnessEach < 0.5] <- 0.5

### 5.2 calculate mean and standard error for each sample
Mean23C <- apply(FitnessEach[,1:5], MARGIN = 1, FUN = mean)
Mean30C <- apply(FitnessEach[,6:10], MARGIN = 1, FUN = mean)
MeanDM <- apply(FitnessEach[,11:13], MARGIN = 1, FUN = mean)
Mean37C <- apply(FitnessEach[,14:16], MARGIN = 1, FUN = mean)

sd23C <- apply(FitnessEach[,1:5], MARGIN = 1, FUN = sd)
sd30C <- apply(FitnessEach[,6:10], MARGIN = 1, FUN = sd)
sdDM  <- apply(FitnessEach[,11:13], MARGIN = 1, FUN = sd)
sd37C <- apply(FitnessEach[,14:16], MARGIN = 1, FUN = sd)

### 5.3 Combine data for output
dataOutEach <- cbind(Data2, FitnessEach, 
                     Mean23C, Mean30C, MeanDM, Mean37C, sd23C, sd30C, sdDM, sd37C)

colnames(dataOutEach) <- c(names(Data2), 
                           "FitT23C1", "FitT23C2", "FitT23C3", "FitT23C4", "FitT23C5",
                           "FitT30C1", "FitT30C2", "FitT30C3", "FitT30C4", "FitT30C5", 
                           "FitD1", "FitD2", "FitD3", 
                           "FitT37C1", "FitT37C2", "FitT37C3", 
                           "Fit23","Fit30","FitDM","Fit37",
                           "Fit23Sd","Fit30Sd","FitDMSd","Fit37Sd")
# output
write.table(x=dataOutEach,file = "ExpDataAllConditions50.txt",quote=FALSE, sep="\t",row.names = F);

dataOutEach2 <- dataOutEach[dataOutEach$T0All>=100,]
write.table(x=dataOutEach2,file = "ExpDataAllConditions100.txt",quote=FALSE, sep="\t",row.names = F);
