## ---------------------------
##
## Script name: Additive Treatments SWD
##
## Purpose of script: calculate power for a two-treatment stepped wedge design with additive treatment effects
##
## Author: Phil Sundin
##
## Date Created: 2020-09-28
##
## Email: phillip1492@ucla.edu
##
## ---------------------------
#####################

## Calculate the power for detecting treatment effects in a
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



library(ggplot2)
#Generate a separate function that calculates power for SWDs
PowerFunction_Additive <- function(RhoW, ModelChoice,n.individuals, delta_1,
              delta_2, RhoA, IAC, Sequence_Tx1,Sequence_Tx2, typeIerror)
{
##variance parameters for repeated cross-sectional design
if (ModelChoice == "RCS")
  { Diagonal <-  (1-RhoW) / n.individuals
OffDiagonal <- (RhoW) }

##variance parameters for cohort design
if (ModelChoice == "Cohort")

{   Diagonal <-  (1 - IAC - RhoW +  IAC*RhoW ) / n.individuals
OffDiagonal <- RhoW + IAC*(1 - RhoW)/ n.individuals }
  
##variance parameters for nested exchangeable model
else if (ModelChoice == "NEC")
{   Diagonal <-  RhoW + (1-RhoW)/n.individuals  - RhoA 
OffDiagonal <- RhoA}
  
#creating data frame for design matrix  
StudyDesign <- data.frame(site = c(1:n.clusters))

StudyDesign$TimeA <- NA
StudyDesign$TimeB <- NA

#assignment of sequences
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

z_1 <-  0
AggregateData_TreatmentASquared_ByCluster <- aggregate(Tx_A ~ site, 
          data=StudyDesign,   FUN = function(x) c(mn = (sum(x))^2))
z_1 <- sum(AggregateData_TreatmentASquared_ByCluster$Tx_A) * 
  OffDiagonal / (Diagonal*( Diagonal + n.periods*OffDiagonal))

z_2 <-  0
AggregateData_TreatmentBSquared_ByCluster <- aggregate(Tx_B ~ site, 
          data=StudyDesign,  FUN = function(x) c(mn = (sum(x))^2))
z_2 <- sum(AggregateData_TreatmentBSquared_ByCluster$Tx_B) * 
  OffDiagonal / (Diagonal*( Diagonal + n.periods*OffDiagonal))

z_3 <-  0
AggregateData_TreatmentABSquared_ByCluster <- aggregate(Tx_AB ~ site, 
          data=StudyDesign,   FUN = function(x) c(mn = (sum(x))^2))
z_3 <- sum(AggregateData_TreatmentABSquared_ByCluster$Tx_AB) * 
  OffDiagonal / (Diagonal*( Diagonal + n.periods*OffDiagonal))

q_part1 <-  sum(StudyDesign$Tx_AB) / Diagonal
q_part2 <- 0
AggregateData_TreatmentA_ByCluster <- aggregate(Tx_A ~ site, 
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))))
AggregateData_TreatmentB_ByCluster <- aggregate(Tx_B~ site, 
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))))
AggregateData_TreatmentAB_ByCluster <- aggregate(Tx_AB~ site, 
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

q_part2 <-  as.numeric(t(AggregateData_TreatmentA_ByCluster$Tx_A) %*% AggregateData_TreatmentB_ByCluster$Tx_B *
OffDiagonal / ((Diagonal*(Diagonal + (n.periods * OffDiagonal)))))

q_1 <- q_part1 - q_part2
AggregateData_TreatmentASquared_ByPeriod <- aggregate(Tx_A ~ Time, 
  data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2/Diagonal^2))
w_1 <- sum(AggregateData_TreatmentASquared_ByPeriod$Tx_A)

AggregateData_TreatmentBSquared_ByPeriod <- aggregate(Tx_B ~ Time, 
  data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2))
w_2 <- sum(AggregateData_TreatmentBSquared_ByPeriod$Tx_B) / Diagonal^2

AggregateData_TxA <- aggregate(Tx_A ~ Time, 
    data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

AggregateData_TxB <- aggregate(Tx_B ~ Time, 
  data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

w_xy <- t(AggregateData_TxA$Tx_A) %*% AggregateData_TxB$Tx_B * (1/Diagonal^2)

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


#get standard error of main effect of treatment 1
SE_Tx1 <- sqrt((solve(PrecisionMatrix))[1,1])
TestStatistic_Tx1 <- delta_1 / SE_Tx1
PowerTx1 <- pnorm(TestStatistic_Tx1 - qnorm(1-typeIerror/2))+
 pnorm(qnorm(typeIerror/2) - TestStatistic_Tx1)

#get standard error of main effect of treatment 1
SE_Tx2 <- sqrt((solve(PrecisionMatrix))[2,2])
TestStatistic_Tx2 <- delta_2 / SE_Tx2
PowerTx2 <- pnorm(TestStatistic_Tx2 - qnorm(1-typeIerror/2))+
 pnorm(qnorm(typeIerror/2) - TestStatistic_Tx2)

return( list(PowerTx1 = PowerTx1, PowerTx2 = PowerTx2)) 
       }


#Example
delta_1 <- delta_2 <- 0.4 
n.periods <- 4
n.individuals <- 15
typeIerror  <- 0.05/2

#create data frame for power
PowerTable <- data.frame(RhoW = c(seq(0,0.4,by=0.01)))

#12-cluster concurrent design
n.clusters <- 12
# Sequencing has 2 clusters transition to treatment 1 at time 2,
# 2 clusters transition to treatment 1 at time 3, and 2 clusters transition
# to treatment 1 at time 3. These clusters never receieve treatment 2
Sequence_Tx1 <- c(2,2,3,3,4,4,NA,NA,NA,NA,NA,NA)

#similar sequencing for treatment 2
Sequence_Tx2 <- c(NA,NA,NA,NA,NA,NA,4,4,3,3,2,2)

Results_Concurrent_12clusters <- mapply(PowerFunction_Additive,
                                 RhoW = PowerTable$RhoW,
                                 ModelChoice = "RCS",
                                 RhoA= NA,
                                 n.individuals =  n.individuals, 
                                 delta_1 = delta_1,
                                 delta_2 = delta_2,
                                 IAC = NA,
                                 Sequence_Tx1 = list(Sequence_Tx1),
                                 Sequence_Tx2 = list(Sequence_Tx2), 
                                 typeIerror = typeIerror)

#Read in power for treatment 1
PowerTable$'Concurrent Design, 12 clusters' <- 
  as.numeric(t((data.frame(Results_Concurrent_12clusters))[1,]))


#calculate power 10 cluster concurrent design
Sequence_Tx1 <- c(2,2,3,3,4, NA,NA,NA,NA,NA)
Sequence_Tx2 <- c(NA,NA,NA,NA,NA,4,3,3,2,2)
n.clusters <- 10
Results_Concurrent_10clusters <-  mapply(PowerFunction_Additive,
                                 RhoW = PowerTable$RhoW,
                                 ModelChoice = "RCS",
                                 RhoA= NA,
                                 IAC = NA,
                                 n.individuals =  n.individuals, 
                                 delta_1 = delta_1,
                                 delta_2 = delta_2,
                                 Sequence_Tx1 = list(Sequence_Tx1),
                                 Sequence_Tx2 = list(Sequence_Tx2), 
                                 typeIerror = typeIerror)

#Read in power for treatment 1
PowerTable$'Concurrent Design, 10 clusters' <- 
  as.numeric(t((data.frame(Results_Concurrent_10clusters))[1,]))


library(reshape2)
PowerTable$'Single Treatment SWDs' <- PowerTable_OneTx$Power
PowerTable <- melt(PowerTable,id="RhoW")

FinalPowerTable <- PowerTable

##Plot power as a function of rho_w
library(latex2exp)
PowerPlot <- ggplot(data=FinalPowerTable[FinalPowerTable$RhoW <= 0.35,],
       aes(x = RhoW, y=value, colour=variable)) +
       geom_line() +
    labs(y= "Power", x = TeX("$\\rho_w$")) + 
 labs(colour="Study Design") +
scale_colour_manual(values = c("#F8766D",
  "#00B6EB","#53B400")) + 
    theme(text = element_text(size=24)) +
 guides(color = guide_legend(override.aes = list(size = 2)))+
	 theme(text = element_text(size=36)) 
plot(PowerPlot)

ggsave("Powerplot_Section31.png",PowerPlot,
       width = 40,height=15,units="cm")
