## ---------------------------
##
## Script name: Linear Contrasts of Treatment Effects in SWD
##
## Purpose of script: calculate power for linear contrasts of treatment effects in a SWD
## 
## Author: Phil Sundin
##
## Date Created: 2020-09-28
##
## Email: phillip1492@ucla.edu
##
## ---------------------------
##
## Calculate the power for detecting linear contrasts of two treatment effects in a
## two-treatment stepped wedge design (SWD) with additive treatment effects. 
## We only consider complete designs.
## Model choice: RCS = repeated cross-sectional, 
## NEC = nested exchangeable, Cohort = Cohort.
## delta_1 and delta_2 are standardized effect sizes for treatments 1 and 2, respectively
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

####################


PowerFunction_AdditiveTreatments <- function(RhoW, ModelChoice,DesignChoice,
n.clusters, RhoA, IAC,Sequence_Tx1,Sequence_Tx2,typeIerror,n.individuals, delta_12)
{
  Sequence_Tx1 <- unlist(Sequence_Tx1)
  Sequence_Tx2 <- unlist(Sequence_Tx2)
##variance parameters for repeated cross-sectional design
if (ModelChoice == "RCS")
  { Diagonal <-  (1-RhoW) / n.individuals
OffDiagonal <- (RhoW) }

##variance parameters for cohort design
else if (ModelChoice == "Cohort")
{   Diagonal <-  RhoW + IAC * (1 - RhoW) / n.individuals
OffDiagonal <- (1 - RhoW - IAC + (IAC * RhoW)) / n.individuals }
  
##variance parameters for nested exchangeable model
else if (ModelChoice == "NEC")
{   Diagonal <-  RhoA
OffDiagonal <- RhoW + (1-RhoW)/n.individuals  - RhoA }
  
StudyDesign <- data.frame(site = c(1:n.clusters))

StudyDesign$TimeA <- Sequence_Tx1
StudyDesign$TimeB <- Sequence_Tx2


##repeat each cluster for each time period
StudyDesign <- 
  StudyDesign[rep(seq_len(nrow(StudyDesign)), each = n.periods), ]
StudyDesign$Time <- rep(1:n.periods)

#Creating treatment indicators
StudyDesign$Tx_A <- 0
StudyDesign$Tx_B <- 0
StudyDesign$Tx_A[StudyDesign$Time >= StudyDesign$TimeA] <- 1
StudyDesign$Tx_B[StudyDesign$Time >= StudyDesign$TimeB] <- 1
StudyDesign$Tx_AB <- StudyDesign$Tx_A * StudyDesign$Tx_B

#####################
##########Calculate precision matrix
#####################

##quantities in paper
f <- n.clusters / 
  (Diagonal + n.periods * OffDiagonal)
g <- n.clusters*OffDiagonal / 
  (Diagonal*(Diagonal + n.periods * OffDiagonal))

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
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))^2))
z_1 <- sum(AggregateData_TreatmentASquared_ByCluster$Tx_A) * 
  OffDiagonal / (Diagonal*(
  Diagonal + n.periods*OffDiagonal
))

z_2 <-  0
AggregateData_TreatmentBSquared_ByCluster <- aggregate(Tx_B ~ site, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))^2))
z_2 <- sum(AggregateData_TreatmentBSquared_ByCluster$Tx_B) * 
  OffDiagonal / (Diagonal*(
  Diagonal + n.periods*OffDiagonal
))

z_3 <-  0
AggregateData_TreatmentABSquared_ByCluster <- aggregate(Tx_AB ~ site, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))^2))
z_3 <- sum(AggregateData_TreatmentABSquared_ByCluster$Tx_AB) * 
  OffDiagonal / (Diagonal*(
  Diagonal + n.periods*OffDiagonal
))

q_part1 <- 0
q_part1 <-  sum(StudyDesign$Tx_AB) /Diagonal

q_part2 <- 0
AggregateData_TreatmentA_ByCluster <- aggregate(Tx_A ~ site, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))))
AggregateData_TreatmentB_ByCluster <- aggregate(Tx_B~ site, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))))
AggregateData_TreatmentAB_ByCluster <- aggregate(Tx_AB~ site, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))))

q_part2 <-  as.numeric(t(AggregateData_TreatmentA_ByCluster$Tx_A) %*% AggregateData_TreatmentB_ByCluster$Tx_B *
OffDiagonal / ((Diagonal*(Diagonal + (n.periods * OffDiagonal)))))

q_1 <- q_part1 - q_part2

w_1 <- 0
AggregateData_TreatmentASquared_ByPeriod <- aggregate(Tx_A ~ Time, 
  data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))^2/Diagonal^2))
w_1 <- sum(AggregateData_TreatmentASquared_ByPeriod$Tx_A)

AggregateData_TreatmentBSquared_ByPeriod <- aggregate(Tx_B ~ Time, 
  data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))^2))
w_2 <- sum(AggregateData_TreatmentBSquared_ByPeriod$Tx_B) / Diagonal^2

AggregateData_TxA <- aggregate(Tx_A ~ Time, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))))

AggregateData_TxB <- aggregate(Tx_B ~ Time, 
          data=StudyDesign, 
      FUN = function(x) c(mn = (sum(x))))


w_xy <- t(AggregateData_TxA$Tx_A) %*% AggregateData_TxB$Tx_B *
  (1/Diagonal^2)


A_11 <-  ((l_1 - z_1)*(f*n.periods*(f+g*n.periods)) -
    y_1^2*(f+g*n.periods) - f*(n.periods*w_1 - l_1^2)) / 
  (f*n.periods*(f+g*n.periods))

A_22 <- ((l_2 - z_2)*(f*n.periods*(f+g*n.periods)) -
    y_2^2*(f+g*n.periods) - f*(n.periods*w_2 - l_2^2)) /
  (f*n.periods*(f+g*n.periods))


A_12 <- q_1 - y_1*y_2/(f*n.periods) - (1/(f+g*n.periods))*(
  w_xy - l_1*l_2/n.periods)

PrecisionMatrix <- matrix(nrow=2,ncol=2)
PrecisionMatrix[1,1] <- A_11
PrecisionMatrix[1,2] <- PrecisionMatrix[2,1] <- A_12
PrecisionMatrix[2,2] <- A_22


InversePrecision2x2 <- (solve(PrecisionMatrix))
SE_TxAMinusTxB <- sqrt(InversePrecision2x2[1,1] + 
  InversePrecision2x2[2,2] - 2 * InversePrecision2x2[2,1])

TestStat_TxAminusTxB <- delta_12 / SE_TxAMinusTxB

PowerTxAMinusTxB <- pnorm(TestStat_TxAminusTxB - qnorm(1-typeIerror/2))+
 pnorm(qnorm(typeIerror/2) - TestStat_TxAminusTxB)

return(PowerTxAMinusTxB)
}

#Examples
delta_12 <- 0.4
n.periods <- 4
n.individuals <- 15

#creating data frame for power as a function of rho_w
PowerTable <- data.frame(RhoW = c(rep(seq(0,0.4,by=0.01),times=3)))
PowerTable$DesignChoice <- c(rep('"Late" Factorial Design, 12 clusters',times=41),
                             rep('"Early" Factorial Design, 10 clusters',times=41),
                             rep("Concurrent Design, 12 clusters",times=41))

PowerTable$n.clusters <- c(rep(12,times=41),
                          rep(10,times=41),
                          rep(12,times=41))

PowerTable$Sequence1 <- c(rep(list(c(2,2,3,3,4,4,4,4,4,4,4,4)),times=41),
                          rep(list(c(2,2,2,3,4,4,3,4,4,4)),times=41),
                          rep(list(c(2,2,3,3,4,4,NA,NA,NA,NA,NA,NA)),times=41)
                          )

PowerTable$Sequence2 <- c(rep(list(c(4,4,4,4,4,4,2,2,3,3,4,4)),times=41),
                          rep(list(c(4,4,4,3,4,4,3,2,2,2)),times=41),
                          rep(list(c(NA,NA,NA,NA,NA,NA,4,4,3,3,2,2)),times=41)
                          )

typeIerror <- 0.05

PowerTable$Power <- mapply(PowerFunction_AdditiveTreatments,
                        RhoW = PowerTable$RhoW,
                        ModelChoice="RCS",
                        IAC = NA,
                        RhoA = NA,
                        n.individuals =  n.individuals,
                        delta_12 = delta_12,
                        DesignChoice = PowerTable$DesignChoice,
                        n.clusters = PowerTable$n.clusters,
                        Sequence_Tx1 = PowerTable$Sequence1, 
                        Sequence_Tx2 = PowerTable$Sequence2,
                        typeIerror = typeIerror) 

##removing unnecessary columns of clusters
PowerTable <- PowerTable[ ,-which(names(PowerTable) %in%
            c("Sequence1","Sequence2","n.clusters"))]

library(reshape2)
FinalPowerTable <- melt(PowerTable,id=c("RhoW","DesignChoice"))

#reorder to levels
FinalPowerTable$DesignChoice <- factor(FinalPowerTable$DesignChoice, 
                      levels = c("Concurrent Design, 12 clusters", 
                                 '"Early" Factorial Design, 10 clusters', 
                                 '"Late" Factorial Design, 12 clusters'))

library(latex2exp)
  
PowerPlot <- ggplot(data=FinalPowerTable[FinalPowerTable$RhoW <= 0.35,],
  aes(x = RhoW, y=value,colour=DesignChoice)) +
  geom_line() +
  labs(y= "Power", x = TeX("$\\rho_w$")) + 
  labs(colour="Study Design") +
	theme(text = element_text(size=24)) +
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(text = element_text(size=36))

plot(PowerPlot)
