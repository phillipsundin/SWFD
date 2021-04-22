## ---------------------------
##
## Script name: Single Treatment SWD
##
## Purpose of script: calculate power for a 
## single treatment stepped wedge design
##
## Author: Phil Sundin
##
## Date Created: 2020-09-28
##
## Email: phillip1492@ucla.edu
##
## ---------------------------
#####################
## Calculate the power for detecting a treatment effect in a
## single treatment stepped wedge design (SWD). 
## We only consider complete designs.
## Model choice: RCS = repeated cross-sectional, 
## NEC = nested exchangeable, Cohort = Cohort.
## delta_1 is standardized effect size for treatment 1.
## Sequence1 is a vector representing the time points (periods) at
## which each cluster transitions from control to treatment condition.
## For example, the 3 cluster, 4 period SWD, 
##     0 1 1 1
##     0 0 1 1
##     0 0 0 1
## where 0 denotes control condition and 1 denotes the treatment condition
## would have Sequence1 <- c(2,3,4)
####################

PowerFunction_OneTx <- function(RhoW, ModelChoice,
        RhoA, IAC, n.clusters, n.periods, n.individuals, delta_1,
        Sequence1, typeIerror)
{
  Sequence1 <- unlist(Sequence1)

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
  
#dataframe used to construct the SWD 
StudyDesign <- data.frame(site = c(1:n.clusters))

## The TimeA variable represents when cluster switches from the control 
## condition to treatment condition
StudyDesign$TimeA <- Sequence1 

##repeat each cluster for each time period
StudyDesign <-  StudyDesign[rep(seq_len(nrow(StudyDesign)), 
                  each = n.periods), ]
StudyDesign$Time <- rep(1:n.periods)

#creating a treatment indicator
StudyDesign$Tx_1 <- 0
StudyDesign$Tx_1[StudyDesign$Time >= StudyDesign$TimeA] <- 1

rownames(StudyDesign) <- 1:nrow(StudyDesign)

###############
#calculate precision and var-cov matrix
###############

##quantities in paper
f <- n.clusters / 
  (Diagonal + n.periods * OffDiagonal)
g <- n.clusters * OffDiagonal / 
  (Diagonal * (Diagonal + n.periods * OffDiagonal))

SumofXijbyperiodbycluster <- sum(StudyDesign$Tx_1)
y_1 <- SumofXijbyperiodbycluster / (Diagonal + n.periods * OffDiagonal)
h_1 <- SumofXijbyperiodbycluster / (Diagonal*(Diagonal + n.periods * OffDiagonal))
l_1 <-  SumofXijbyperiodbycluster / Diagonal 

z_1 <-  0
AggregateData_TreatmentASquared_ByCluster <- aggregate(Tx_1 ~ site, 
data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2))

z_1 <- sum(AggregateData_TreatmentASquared_ByCluster$Tx_1) * 
  OffDiagonal / (Diagonal * (Diagonal + n.periods*OffDiagonal))

w_1 <- 0
AggregateData_TreatmentASquared_ByPeriod <- aggregate(Tx_1 ~ Time, 
  data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2/Diagonal^2))
w_1 <- sum(AggregateData_TreatmentASquared_ByPeriod$Tx_1)

AggregateData_TxA <- aggregate(Tx_1 ~ Time, 
data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

A_11 <-  ((l_1 - z_1) * (f * n.periods * (f + g * n.periods)) -
    y_1^2 * (f + g * n.periods) - f * (n.periods * w_1 - l_1^2)) / 
  (f * n.periods * (f + g * n.periods))

InversePrecision <- solve(A_11)

#calculate power
SE_Tx1 <- sqrt(InversePrecision)
TestStat_Tx1 <- delta_1 / SE_Tx1

PowerTxA <- pnorm(TestStat_Tx1 - qnorm(1-typeIerror/2)) +
                  pnorm(qnorm(typeIerror/2) - TestStat_Tx1)
return(PowerTxA)
}


#Example
delta_1 <- 0.4

#6 clusters total and 4 time periods
n.clusters <- 6
n.periods <- 4

#2 clusters transition from control to treatment at time 2,
#2 clusters to transition at time 3, and 2 clusters transition at time 4
Sequence1 <-  c(2,2,3,3,4,4) 

n.individuals <- 15
typeIerror <- 0.05

#creating a data frame for power
PowerTable_OneTx <- data.frame(RhoW = c(seq(0,0.4,by=0.01)))
PowerTable_OneTx$ModelChoice <- "RCS"
PowerTable_OneTx$Sequence1 <- c(rep(list(c(2,2,3,3,4,4))))
                                   
#calculate power as a function of rho_w for the repeated cross-sectional design
PowerTable_OneTx$Power <- mapply(PowerFunction_OneTx,
                                 RhoW = PowerTable_OneTx$RhoW,
                                 ModelChoice = PowerTable_OneTx$ModelChoice,
                                 RhoA= NA,
                                 IAC = NA,
                                 n.clusters = n.clusters, 
                                 delta_1 = delta_1,
                                 n.individuals = n.individuals,
                                 n.periods = n.periods,
                                 Sequence1 = PowerTable_OneTx$Sequence1,
                                 typeIerror = typeIerror
                                 )

#plot power as a function of ICC
plot(PowerTable_OneTx$RhoW,PowerTable_OneTx$Power, type="l")
