## ---------------------------
##
## Script name:  Multiarm SWDs
##
## Purpose of script: calculate power for multi-arm SWDs
##
## Author: Phil Sundin
##
## Date Created: 2021-04-22
##
## Email: phillip1492@ucla.edu
##
## ---------------------------
##
## Calculate the power for detecting main treatment effects for
## multi-arm stepped wedge factorial design (SWFD) 
## We only consider complete designs.
## Model choice: RCS = repeated cross-sectional, 
## NEC = nested exchangeable, Cohort = Cohort.
## delta_1, delta_2, and delta_3 are standardized effect sizes for 
## the main effect of treatments 1 and 2 and 3, respectively
## Sequence1 is a vector representing the time points (periods) at
## which each cluster transitions from control to treatment 1.
## Sequence2 is a vector representing the time points (periods) at
## which each cluster transitions from control to treatment 2.
## For example, the 3 cluster, 4 period SWD, 
##     0 1 1   1 
##     0 0 1+2 1+2
##     0 0 2   1+2
## where 0 denotes control condition, 1 denotes the condition with only 
## treatment 1, 2 denote the condition with only treatment 2, and 1+2 denotes the condition
## where a cluster receives both treatment 1 and 2
## would have Sequence1 <- c(2,3,4) and Sequence2 <- c(NA,3,3)
## -----------


PowerFunction_MultiArm <- function(RhoW, ModelChoice,n.individuals,  delta_1,
delta_2,delta_3, RhoA, IAC, typeIerror, Sequence_Tx1,Sequence_Tx2)
{
  Sequence_Tx1 <- unlist(Sequence_Tx1)
  Sequence_Tx2 <- unlist(Sequence_Tx2)

  ##variance parameters for repeated cross-sectional design
if (ModelChoice == "RCS")
  { Diagonal <-  (1-RhoW) / n.individuals
OffDiagonal <- (RhoW) }

##variance parameters for cohort design
if (ModelChoice == "Cohort")

{   Diagonal <-  (1 - IAC - RhoW +  IAC*RhoW ) / n.individuals
OffDiagonal <- RhoW + IAC*(1 - RhoW)/ n.individuals }
  
##variance parameters for nested exchangeable model
 if (ModelChoice == "NEC")
{   Diagonal <-  RhoW + (1-RhoW)/n.individuals  - RhoA 
OffDiagonal <- RhoA}
  
StudyDesign <- data.frame(site = c(1:n.clusters))

StudyDesign$TimeA <- Sequence_Tx1
StudyDesign$TimeB <- Sequence_Tx2

##repeat each cluster for each time period
StudyDesign <- 
  StudyDesign[rep(seq_len(nrow(StudyDesign)), each = n.periods), ]
StudyDesign$Time <- rep(1:n.periods)

StudyDesign$Tx_A <- 0
StudyDesign$Tx_B <- 0
StudyDesign$Tx_A[StudyDesign$Time >= StudyDesign$TimeA] <- 1
StudyDesign$Tx_B[StudyDesign$Time >= StudyDesign$TimeB] <- 1
StudyDesign$Tx_AB <- StudyDesign$Tx_A * StudyDesign$Tx_B

##treating combined condition as separate
StudyDesign$Tx_A[StudyDesign$Tx_AB == 1] <- 0
StudyDesign$Tx_B[StudyDesign$Tx_AB == 1] <- 0


rownames(StudyDesign) <- 1:nrow(StudyDesign)

#####################
##########Calculate precision matrix
#####################

##quantities in paper
f <- n.clusters /   (Diagonal + n.periods * OffDiagonal)
g <- n.clusters*OffDiagonal / (Diagonal*(Diagonal + n.periods * OffDiagonal))

SumofXijbyperiodbycluster <- sum(StudyDesign$Tx_A)
y_1 <- SumofXijbyperiodbycluster / (Diagonal + n.periods * OffDiagonal)
h_1 <- SumofXijbyperiodbycluster / (Diagonal*(Diagonal + n.periods * OffDiagonal))
l_1 <-  SumofXijbyperiodbycluster / Diagonal 

SumofYijbyperiodbycluster <- sum(StudyDesign$Tx_B)
y_2 <- SumofYijbyperiodbycluster / (Diagonal + n.periods * OffDiagonal)
h_2 <- SumofYijbyperiodbycluster / (Diagonal*(Diagonal + n.periods * OffDiagonal))
l_2 <- SumofYijbyperiodbycluster / Diagonal

SumofXYijbyperiodbycluster <- sum(StudyDesign$Tx_AB)
y_3 <- SumofXYijbyperiodbycluster / (Diagonal + n.periods * OffDiagonal)
h_3 <- SumofXijbyperiodbycluster / (Diagonal*(Diagonal + n.periods * OffDiagonal))
l_3 <-  SumofXYijbyperiodbycluster / Diagonal 

z_1 <-  0
AggregateData_TreatmentASquared_ByCluster <- aggregate(Tx_A ~ site, 
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2))
z_1 <- sum(AggregateData_TreatmentASquared_ByCluster$Tx_A) * 
  OffDiagonal / (Diagonal*(Diagonal + n.periods*OffDiagonal))

z_2 <-  0
AggregateData_TreatmentBSquared_ByCluster <- aggregate(Tx_B ~ site, 
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2))
z_2 <- sum(AggregateData_TreatmentBSquared_ByCluster$Tx_B) * 
  OffDiagonal / (Diagonal*(Diagonal + n.periods*OffDiagonal))

z_3 <-  0
AggregateData_TreatmentABSquared_ByCluster <- aggregate(Tx_AB ~ site, 
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2))
z_3 <- sum(AggregateData_TreatmentABSquared_ByCluster$Tx_AB) * 
  OffDiagonal / (Diagonal*(Diagonal + n.periods*OffDiagonal))


AggregateData_TreatmentA_ByCluster <- aggregate(Tx_A ~ site, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))))
AggregateData_TreatmentB_ByCluster <- aggregate(Tx_B~ site, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))))
AggregateData_TreatmentAB_ByCluster <- aggregate(Tx_AB~ site, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))))


q_part1 <- 0 # sum(StudyDesign$Tx_AB) / Diagonal

q_part2 <-  as.numeric(t(AggregateData_TreatmentA_ByCluster$Tx_A) %*% AggregateData_TreatmentB_ByCluster$Tx_B *
OffDiagonal / ((Diagonal*(Diagonal + (n.periods * OffDiagonal)))))

q_1 <- q_part1 - q_part2

q_2_part2 <- as.numeric(t(AggregateData_TreatmentA_ByCluster$Tx_A) %*% AggregateData_TreatmentAB_ByCluster$Tx_AB *      
                     OffDiagonal/(Diagonal*(Diagonal+n.periods*OffDiagonal)))

q_3_part2 <- as.numeric(t(AggregateData_TreatmentB_ByCluster$Tx_B) %*% AggregateData_TreatmentAB_ByCluster$Tx_AB *      
                     OffDiagonal/(Diagonal*(Diagonal+n.periods*OffDiagonal)))

q_2 <- #l_3 
  0 - q_2_part2
q_3 <- #l_3
  0 - q_3_part2

AggregateData_TreatmentASquared_ByPeriod <- aggregate(Tx_A ~ Time, 
  data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))^2/Diagonal^2))
w_1 <- sum(AggregateData_TreatmentASquared_ByPeriod$Tx_A)

AggregateData_TreatmentBSquared_ByPeriod <- aggregate(Tx_B ~ Time, 
  data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))^2))
w_2 <- sum(AggregateData_TreatmentBSquared_ByPeriod$Tx_B) / Diagonal^2

AggregateData_TreatmentABSquared_ByPeriod <- aggregate(Tx_AB ~ Time, 
  data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))^2))
w_3 <- sum(AggregateData_TreatmentABSquared_ByPeriod$Tx_AB) / Diagonal^2

AggregateData_TxA <- aggregate(Tx_A ~ Time, 
data=StudyDesign,  FUN = function(x) c(mn = (sum(x))))

AggregateData_TxB <- aggregate(Tx_B ~ Time, 
      data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

AggregateData_TxAB <- aggregate(Tx_AB ~ Time, 
data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

w_xy <- t(AggregateData_TxA$Tx_A) %*% AggregateData_TxB$Tx_B * (1/Diagonal^2)

w_x_xy <- t(AggregateData_TxA$Tx_A) %*% AggregateData_TxAB$Tx_AB *  (1/Diagonal^2)

w_y_xy <- t(AggregateData_TxB$Tx_B) %*% AggregateData_TxAB$Tx_AB *  (1/Diagonal^2)

A_11 <-  ((l_1 - z_1)*(f*n.periods*(f+g*n.periods)) -
    y_1^2*(f+g*n.periods) - f*(n.periods*w_1 - l_1^2)) / (f*n.periods*(f+g*n.periods))

A_22 <- ((l_2 - z_2)*(f*n.periods*(f+g*n.periods)) -
    y_2^2*(f+g*n.periods) - f*(n.periods*w_2 - l_2^2)) /  (f*n.periods*(f+g*n.periods))

A_33 <- ((l_3 - z_3)*(f*n.periods*(f+g*n.periods)) -
    y_3^2*(f+g*n.periods) - f*(n.periods*w_3 - l_3^2)) /  (f*n.periods*(f+g*n.periods))


A_12 <- q_1 - y_1*y_2/(f*n.periods) - (1/(f+g*n.periods))*(
  w_xy - l_1*l_2/n.periods)
A_13 <- q_2 - y_1*y_3/(f*n.periods) - 
  (1/(f+g*n.periods))*(w_x_xy - l_1*l_3/n.periods)
A_23 <- q_3 - y_2*y_3/(f*n.periods) - 
  (1/(f+g*n.periods))*(w_y_xy - l_2*l_3/n.periods)


PrecisionMatrix <- matrix(nrow=3,ncol=3)
PrecisionMatrix[1,1] <- A_11
PrecisionMatrix[1,2] <-  PrecisionMatrix[2,1] <- A_12
PrecisionMatrix[2,2] <- A_22
PrecisionMatrix[1,3] <- PrecisionMatrix[3,1] <- A_13
PrecisionMatrix[3,3] <- A_33
PrecisionMatrix[2,3] <- PrecisionMatrix[3,2] <- A_23

InversePrecision <- solve (PrecisionMatrix)

ResultsTable <- matrix(0,nrow=3,ncol=3)
colnames(ResultsTable) <- c("Main Effect of Treatment 1",
                    "Main Effect of Treatment 2",
                    "Main Effect of Treatment 3")
rownames(ResultsTable) <- c("Standard Error (Sqrt of Variances)",
                    "ResultsTable Statistic (Estimated Value / SE)",
                     "Power with alpha=0.05 "                    )

##Reading in the Standard Errors
ResultsTable[1,1] <- round(sqrt(InversePrecision[1,1]),4)
ResultsTable[1,2] <- round(sqrt(InversePrecision[2,2]),4)
ResultsTable[1,3] <- round(sqrt(InversePrecision[3,3]),4)

ResultsTable[2,1] <- delta_1 / ResultsTable[1,1]
ResultsTable[2,2] <- delta_2 / ResultsTable[1,2]
ResultsTable[2,3] <- delta_3 / ResultsTable[1,3]
   
############calculating power for treatment A
ResultsTable[3,1] <- pnorm(ResultsTable[2,1]  - qnorm(1-typeIerror/2),)+
 pnorm(qnorm(typeIerror/2) - ResultsTable[2,1] )
 
  ############calculating power for treatment B
ResultsTable[3,2] <- pnorm(ResultsTable[2,2] - qnorm(1-typeIerror/2))+
pnorm(qnorm(typeIerror/2) - ResultsTable[2,2])

  ############calculating power for treatment AB
ResultsTable[3,3] <- pnorm(ResultsTable[2,3] - qnorm(1-typeIerror/2),)+
  pnorm(qnorm(typeIerror/2) - ResultsTable[2,3])
 
return(ResultsTable)
}

require(ggplot2)


#Example
delta_1 <- delta_2 <- delta_3 <- 0.6
n.clusters <- 8
n.periods <- 5
n.individuals <- 15
typeIerror <- 0.05/3

#Creating data frame for calculating power
PowerTable <- data.frame(IAC = c(rep(seq(0,0.3,by=0.01),times=4)))
PowerTable$RhoW = 0.1
PowerTable$RhoA = 0.15

#each separate list corresponds to different sequencing of treatments for both 
#treatments 1 and 2, as seen in Section 3.3 of the paper
##the condition with both treatments is treated as a separate intervention
  PowerTable$Sequence1 <- c(rep(list(c(2,3,2,3,4,5,NA,NA)),times=31),
                                 rep(list(c(2,2,3,4,NA,5,5,4)),times=31),
                                 rep(list(c(2,3,3,4,4,5,5,5)),times=31),
                                 rep(list(c(2,3,4,3,3,5,5,NA)),times=31
                                ))
                     
  PowerTable$Sequence2 <- c(rep(list(c(NA,NA,2,3,4,5,3,2)),times=31),
                                 rep(list(c(3,4,5,NA,4,4,3,2)),times=31),
                                 rep(list(c(5,5,5,4,4,3,3,2)),times=31),
                                 rep(list(rev(c(2,3,4,3,3,5,5,NA))),times=31
                                ))
                            
#Determine power for each value of RhoW and each design 
Results <-  mapply(PowerFunction_MultiArm,
        RhoW = PowerTable$RhoW,
        ModelChoice = "Cohort",
        RhoA= PowerTable$RhoA,
        IAC = PowerTable$IAC,
        n.individuals =  n.individuals,
        delta_1 = delta_1,
        delta_2 = delta_2,
        delta_3 = delta_3,
        typeIerror = typeIerror,
        Sequence_Tx1 = PowerTable$Sequence1,
        Sequence_Tx2 = PowerTable$Sequence2)

#read in appropriate value of the results table
PowerTable$Power_Treatment1 <- Results[3,]
PowerTable$Power_Treatment2 <- Results[6,]
PowerTable$Power_Treatment3 <- Results[9,]

PowerTable$DesignChoice <- c(rep(1,times=31),
                                 rep(2,times=31),
                                 rep(3,times=31),
                                 rep(4,times=31)
                                 )
#removing unnecessary lines
PowerTable <- PowerTable[ , -which(names(PowerTable) %in% c("Sequence1","Sequence2"))]

library(reshape2)

FinalPowerTable <- melt(PowerTable,id=c("IAC","DesignChoice"))
library(latex2exp)

FinalPowerTable$DesignChoice <-
    as.character(FinalPowerTable$DesignChoice)

##Plot for main effects
PowerTable_MainEffect1 <- 
  FinalPowerTable[FinalPowerTable$variable == "Power_Treatment1" ,]
PowerTable_MainEffect1$variable <- as.character(PowerTable_MainEffect1$variable)

PowerTable_MainEffect1$variable[  PowerTable_MainEffect1$DesignChoice == 1 ] <-  "Main Effect Design #1"
PowerTable_MainEffect1$variable[  PowerTable_MainEffect1$DesignChoice == 2 ] <- "Main Effect Design #2"
PowerTable_MainEffect1$variable[  PowerTable_MainEffect1$DesignChoice == 3 ] <- "Main Effect Design #3"
PowerTable_MainEffect1$variable[  PowerTable_MainEffect1$DesignChoice == 4] <- "Main Effect Design #4"


PowerPlot  <- ggplot(data= 
          PowerTable_MainEffect1[PowerTable_MainEffect1$IAC <= 0.30,], 
          aes(x = IAC, y=value,
  colour=variable )) +
  geom_line() +
  guides(color = guide_legend(override.aes = list(size = 2)))+
  labs(y= "Power", x = TeX("$IAC$")) + 
  scale_y_continuous(limits=c(0.7, 1)) + 
  labs(color= "Study Design") +
  scale_color_manual(values=
  c("green1","dodgerblue","red","magenta2")) +
  theme(text = element_text(size=24)) 
plot(PowerPlot)
  
ggsave("Powerplot_section34_tx1.png",PowerPlot,
       width = 40,height=15,units="cm") 




##Plot for main effects 2
PowerTable_MainEffect2 <- 
  FinalPowerTable[FinalPowerTable$variable == "Power_Treatment2" ,]
PowerTable_MainEffect2$variable <- as.character(PowerTable_MainEffect2$variable)

PowerTable_MainEffect2$variable[  PowerTable_MainEffect2$DesignChoice == 1 ] <-  "Main Effect Design #1"
PowerTable_MainEffect2$variable[  PowerTable_MainEffect2$DesignChoice == 2 ] <- "Main Effect Design #2"
PowerTable_MainEffect2$variable[  PowerTable_MainEffect2$DesignChoice == 3 ] <- "Main Effect Design #3"
PowerTable_MainEffect2$variable[  PowerTable_MainEffect2$DesignChoice == 4] <- "Main Effect Design #4"


PowerPlot  <- ggplot(data= 
          PowerTable_MainEffect2[PowerTable_MainEffect2$IAC <= 0.30,], 
          aes(x = IAC, y=value,
  colour=variable )) +
  geom_line() +
  guides(color = guide_legend(override.aes = list(size = 2)))+
  labs(y= "Power", x = TeX("$IAC$")) + 
  scale_y_continuous(limits=c(0.7, 0.95)) + 
  labs(color= "Study Design") +
  scale_color_manual(values=
  c("green1","dodgerblue","red","magenta2")) +
  theme(text = element_text(size=24)) 
  plot(PowerPlot)
  
ggsave("Powerplot_section34_tx2.png",PowerPlot,
       width = 40,height=15,units="cm") 



##Plot for main effects 3
PowerTable_MainEffect3 <- 
  FinalPowerTable[FinalPowerTable$variable == "Power_Treatment3" ,]
PowerTable_MainEffect3$variable <- as.character(PowerTable_MainEffect3$variable)

PowerTable_MainEffect3$variable[  PowerTable_MainEffect3$DesignChoice == 1 ] <-  "Main Effect Design #1"
PowerTable_MainEffect3$variable[  PowerTable_MainEffect3$DesignChoice == 2 ] <- "Main Effect Design #2"
PowerTable_MainEffect3$variable[  PowerTable_MainEffect3$DesignChoice == 3 ] <- "Main Effect Design #3"
PowerTable_MainEffect3$variable[  PowerTable_MainEffect3$DesignChoice == 4] <- "Main Effect Design #4"


PowerPlot  <- ggplot(data= 
          PowerTable_MainEffect3[PowerTable_MainEffect3$IAC <= 0.30,], 
          aes(x = IAC, y=value,
  colour=variable )) +
  geom_line() +
  guides(color = guide_legend(override.aes = list(size = 2)))+
  labs(y= "Power", x = TeX("$IAC$")) + 
  scale_y_continuous(limits=c(0.5, 1)) + 
  labs(color= "Study Design") +
  scale_color_manual(values=
  c("green1","dodgerblue","red","magenta2")) +
  theme(text = element_text(size=24)) 
  plot(PowerPlot)
  
ggsave("Powerplot_section34_tx3.png",PowerPlot,
       width = 40,height=15,units="cm") 

