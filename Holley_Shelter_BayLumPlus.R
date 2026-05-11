
library(BayLumPlus)
library(rjags)
library(ggplot2)

### HOLLEY SHELTER - BayLumPLUS ###

## NORTH COAST 14C AND OSL##

# Set the working directory
setwd("C:/Users/Amelie/Documents/Holley shelter/North Coast/14C_OSL")

# STEP 1: Bayesian Computation of unconstrained ages #

# Radiocarbon ages and luminescence De (Q + F)
Holley_NC_Meas=c(30730,30300,
                 76.5,79.2,
                 29200,30330,31380,30050,30640,29880,30980,31750,32710,30980,33300,
                 35600,31800,33300,36300,31400,37000,34300,37500,37500,37500,37000,35400,33800,34900,
                 37500,37200,37200,43500,38500,39300,43000,43500,41200,47000,47000
)

# Luminescence De err (Q + F)
Holley_NC_dosesErr=c(2.7,9.4)

# Luminescence Dr (Q + F)
Holley_NC_doseRates=c(1.92,1.96)

# Radiocarbon err
C14_Er=c(370,520,1500,330,470,330,360,330,350,370,590,350,1800,2300,1500,1800,2600,1500,
         2800,2000,2900,2800,2700,2800,1600,1400,1500,2800,2000,2700,5000,2200,3300,3500,
         3000,3800,6500,4300)

# Create a list 
dt=list(M = Holley_NC_Meas, sD = Holley_NC_dosesErr, ddot = Holley_NC_doseRates,sigmaC14Cal=C14_Er)

# Assign OSL and C14 data (0 = 14C, 1 = OSL)
OSLorC14=c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0)

# Sample Names in stratigraphical order (14C and OSL)
SampleNames=c("118SQ12","120SQ12","HRS7-Q","HRS7-F","240SQ12","241SQ12","364SQ12","386SQ12","388SQ12","390SQ12","463SQ12","492SQ12",
              "518SQ12","609SQ12","966HSD12_57DF","967HSD12_57DF","602GHM12_56EDG","601GHM12_56EDG","1059GHM12_57EDG","633GHM12_56FG",
              "1089GHM12_57FC","1111GHM12_57FC","693GDB12_56GH","723GHM12_56GH","1246NS12_57GH",
              "1276EMC12_57GH", "809DGB12_56GH","929GHM12_56HW","989GDB12_56HW","1546GHM12_57HW","1191GDB12_56HW","1189GDB12_56HW",
              "1612GDB12_57KM","1636GDB12_57LM","1647GDB12_57_1M","HS1443SQ12_56JP", "HS1484SQ12_56JP","HS1513SQ12_56JP",
              "HS1578SQ12_56KL","HS1603SQ12_56KL")
nb_sample=length(SampleNames)

# Assign systematic errors in %
sigma_s =  c(s_betaK = 0.11,
             s_betaU = 0.11,
             s_betaTh = 0.11,
             s_gammaK = 0,
             s_gammaU = 0,
             s_gammaTh = 0,
             s_gammaDR = 0.05,
             s_CAL = 0.04,
             s_intDR = 0)

# Path to the .csv file with corresponding Dr values
PathThetaInput=c("C:/Users/Amelie/Documents/Holley shelter/North Coast/14C_OSL/")

### NOTE: for the North Coast - only two rows for the original matrix => error when running create_ThetaMatrix() ###
## Solution: duplicate the original matrix to obtain 4 rows (2x2) and crop to obtain the two original rows ###

# Read the original 2-row matrix from the .csv file
ThetaInputMatrix <- read.csv("Holley_NC_matrix.csv", sep=",")

# Duplicate the rows to avoid the 2-sample bug
ThetaInputMatrix_expanded <- rbind(ThetaInputMatrix, ThetaInputMatrix)

# Run create_ThetaMatrix on the expanded matrix
Theta_expanded <- create_ThetaMatrix(input = ThetaInputMatrix_expanded, 
                                     output_file = NULL, 
                                     sigma_s = sigma_s)

# Extract only the first 2×2 submatrix (the real samples)
Theta <- Theta_expanded[1:2, 1:2]

# Fix the row and column names
rownames(Theta) <- colnames(Theta) <- ThetaInputMatrix$SAMPLE_ID

# Check the submatrix values
print(Theta)

# Model Bayesian ages without stratigraphical constraints
Holley_NC_14C_OSL=Compute_AgeS_DC14(DATA=dt, Nb_sample = nb_sample, SampleNames = SampleNames,
                                    encoding = OSLorC14, ThetaMatrix = Theta, prior = "Unconstrained", 
                                    CalibrationCurve = "C:/Users/Amelie/Documents/Holley shelter/SHCal20.csv",
                                    PriorAge =rep(c(1,60),nb_sample), Iter = 10000, burnin = 5000, t = 5)

#### Optional steps
# Save the output as an .RData file to be re-used if necessary
save(Holley_NC_14C_OSL, file = "MCMC_Holley_NC_14C-OSL_NoStratiNoCov.RData")

# Load the ouput to re-run the next steps (Isotonic Distorsion) if necessary
load("MCMC_Holley_NC_14C-OSL_NoStratiNoCov.RData")
#### End Optional steps

# Plot and save the unconstrained ages
AgeWith_NoConstraints = plot_Ages(Holley_NC_14C_OSL, plot_mode = "density")
p <- recordPlot()
pdf("Holley_NorthCoast_Ages_unconstrained.pdf", width = 10, height = 8)
replayPlot(p)
dev.off()

# STEP 2: Isotonic Regression on Constrained Ages #

# Define a matrix with stratigraphic constraints
SC = matrix(data=0,ncol=40,nrow=41)
SC[1,]=rep(1,40)
SC[2,]=c(rep(0,2),rep(1,38)) # 118SQ12
SC[3,]=c(rep(0,2),rep(1,38)) # 120SQ12
SC[4,]=c(rep(0,14),rep(1,26)) # HRS7-Q
SC[5,]=c(rep(0,14),rep(1,26)) # HRS7-F
SC[6,]=c(rep(0,14),rep(1,26)) # 240SQ12
SC[7,]=c(rep(0,14),rep(1,26)) # 241SQ12
SC[8,]=c(rep(0,14),rep(1,26)) # 364SQ12
SC[9,]=c(rep(0,14),rep(1,26)) # 386SQ12
SC[10,]=c(rep(0,14),rep(1,26)) # 388SQ12
SC[11,]=c(rep(0,14),rep(1,26)) # 390SQ12
SC[12,]=c(rep(0,14),rep(1,26)) # 463SQ12
SC[13,]=c(rep(0,14),rep(1,26)) # 492SQ12
SC[14,]=c(rep(0,14),rep(1,26)) # 518SQ12
SC[15,]=c(rep(0,14),rep(1,26)) # 609SQ12
SC[16,]=c(rep(0,16),rep(1,24)) # 966HSD12_57DF
SC[17,]=c(rep(0,16),rep(1,24)) # 967HSD12_57DF
SC[18,]=c(rep(0,19),rep(1,21)) # 602GHM12_56EDG
SC[19,]=c(rep(0,19),rep(1,21)) # 601GHM12_56EDG
SC[20,]=c(rep(0,19),rep(1,21)) # 1059GHM12_57EDG
SC[21,]=c(rep(0,22),rep(1,18)) # 633GHM12_56FG
SC[22,]=c(rep(0,22),rep(1,18)) # 1089GHM12_57FC
SC[23,]=c(rep(0,22),rep(1,18)) # 1111GHM12_57FC
SC[24,]=c(rep(0,24),rep(1,16)) # 693GDB12_56GH
SC[25,]=c(rep(0,24),rep(1,16)) # 723GHM12_56GH
SC[26,]=c(rep(0,27),rep(1,13)) # 1246NS12_57GH
SC[27,]=c(rep(0,27),rep(1,13)) # 1276EMC12_57GH
SC[28,]=c(rep(0,27),rep(1,13)) # 809DGB12_56GH
SC[29,]=c(rep(0,33),rep(1,7)) # 929GHM12_56HW
SC[30,]=c(rep(0,33),rep(1,7)) # 989GDB12_56HW
SC[31,]=c(rep(0,33),rep(1,7)) # 1546GHM12_57HW
SC[32,]=c(rep(0,33),rep(1,7)) # 1191GDB12_56HW
SC[33,]=c(rep(0,33),rep(1,7)) # 1189GDB12_56HW
SC[34,]=c(rep(0,33),rep(1,7)) # 1612GDB12_57KM
SC[35,]=c(rep(0,35),rep(1,5)) # 1636GDB12_57LM
SC[36,]=c(rep(0,35),rep(1,5)) # 1647GDB12_57_1M
SC[37,]=c(rep(0,38),rep(1,2)) # HS1443SQ12_56JP
SC[38,]=c(rep(0,38),rep(1,2)) # HS1484SQ12_56JP
SC[39,]=c(rep(0,38),rep(1,2)) # HS1513SQ12_56JP
SC[40,]=c(rep(0,40)) # HS1578SQ12_56KL
SC[41,]=c(rep(0,40)) # HS1603SQ12_56KL

# Isotonic Regression taking into account the Strati Matrix
IsotonicDistorsion = IsotonicCurve(StratiConstraints = SC, Holley_NC_14C_OSL,interactive = FALSE)

# Save the Strati matrix as a .pdf figure
p_strati <- recordPlot()
pdf("Holley_NorthCoast_StratiMatrix.pdf", width = 10, height = 8)
replayPlot(p_strati)
dev.off()

# Visualize the computed ages
knitr::kable(IsotonicDistorsion$Ages )

# Plot and save the Isotonic Regression ages as a .pdf
AgeWith_IR <-  plot_Ages(object = IsotonicDistorsion, plot_mode = "density")
p <- recordPlot(AgeWith_IR)
pdf("Holley_NorthCoast_Ages_IR.pdf", width = 10, height = 8)
replayPlot(p)
dev.off()

# Plot the Isotonic Curve
HolleyPlots = PlotIsotonicCurve(object = IsotonicDistorsion, level = .95)

# Save the Harris diagram as a .png figure
HolleyPlots$DAG
ggsave("Holley_NorthCoast_DAG.png", 
       plot = HolleyPlots$DAG, 
       width = 10, height = 8,
       dpi = 300)

# Save the Isotonic curve as a .png figure
HolleyPlots$ribbon
ggsave("Holley_NorthCoast_IsotonicCurve.png", 
       plot = HolleyPlots$ribbon, 
       width = 10, height = 8,
       dpi = 300)

# Save ages in a .csv file
write.csv(IsotonicDistorsion$Ages,
          file = "Holley_NorthC_14C-OSL_IsotonicRegression.csv",
          row.names = FALSE)

# Save MCMC values in a .csv file
write.csv(as.matrix(IsotonicDistorsion[["Sampling"]][[1]]),
          file = "Holley_NorthC_IsotonicRegression_MCMC.csv",
          row.names = FALSE)

## CRAMB PROFILE - C14 and OSL ##

# Set the working directory
setwd("C:/Users/Amelie/Documents/Holley shelter/Cramb Profile/14C_OSL")

# STEP 1: Bayesian Computation of unconstrained ages #

# Radiocarbon ages and luminescence De (Q + F)
Holley_CP_14C_OSL_Meas=c(36050, 32580, 34090, 35070, 96.1, 103.1, 38190, 35850, 36440, 34980, 36990, 101.4, 98.9)

# Luminescence De err (Q + F)
Holley_CP_14C_OSL_dosesErr=c(3.5, 5.2, 3.9, 6.1)

# Luminescence Dr (Q + F)
Holley_CP_14C_OSL_doseRates=c(2.69, 2.73, 2.32, 2.38)

# Radiocarbon err
C14_Er=c(740, 410, 530, 780, 760, 590, 680, 540, 1829)

# Create a list 
dt=list(M = Holley_CP_14C_OSL_Meas, sD = Holley_CP_14C_OSL_dosesErr, ddot = Holley_CP_14C_OSL_doseRates,sigmaC14Cal=C14_Er)

# Assign OSL and C14 data (0 = 14C, 1 = OSL)
OSLorC14=c(0,0,0,0,1,1,0,0,0,0,0,1,1)

# Sample Names in stratigraphical order (14C and OSL)
SampleNames=c("109SQ10", "130SQ10", "170SQ10", "173SQ10", "HRS5-Q","HRS5-F","376SQ10", "388SQ10",
              "415SQ10", "445SQ10", "455SQ10","HRS6-Q","HRS6-F")
nb_sample=length(SampleNames)

# Assign systematic errors in %
sigma_s =  c(s_betaK = 0.11,
             s_betaU = 0.11,
             s_betaTh = 0.11,
             s_gammaK = 0,
             s_gammaU = 0,
             s_gammaTh = 0,
             s_gammaDR = 0.05,
             s_CAL = 0.04,
             s_intDR = 0)

# Path to the .csv file with corresponding Dr values
PathThetaInput=c("C:/Users/Amelie/Documents/Holley shelter/Cramb Profile/14C_OSL/")
#ThetaInput = paste(PathThetaInput,"Holley_CP_matrix.csv",sep="")

# Read the original matrix from the .csv file
ThetaInputMatrix <- read.csv("Holley_CP_matrix.csv", sep=",")

# Create a matrix accounting for shared systematic uncertainties on the dose rate components
Theta=create_ThetaMatrix(input=ThetaInputMatrix,output_file = NULL, sigma_s=sigma_s)

# Check the matrix values
print(Theta)

# Model Bayesian ages without stratigraphical constraints
Holley_CP_14C_OSL=Compute_AgeS_DC14(DATA=dt, Nb_sample = nb_sample, SampleNames = SampleNames,
                                    encoding = OSLorC14, ThetaMatrix = Theta, prior = "Unconstrained", 
                                    CalibrationCurve = "C:/Users/Amelie/Documents/Holley shelter/SHCal20.csv", jags_method = "rjags",
                                    PriorAge =rep(c(1,60),nb_sample), Iter = 10000, burnin = 5000, t = 5)

#### Optional steps
# Save the output as an .RData file to be re-used if necessary
save(Holley_CP_14C_OSL, file = "MCMC_Holley_CP_14C-OSL_NoStratiNoCov.RData")

# Load the ouput to re-run the next steps (Isotonic Distorsion) if necessary
load("MCMC_Holley_CP_14C-OSL_NoStratiNoCov.RData")
#### End Optional steps

# Plot and save the unconstrained ages
AgeWith_NoConstraints = plot_Ages(Holley_CP_14C_OSL, plot_mode = "density")
p <- recordPlot()
pdf("Holley_CP_Ages_unconstrained.pdf", width = 10, height = 8)
replayPlot(p)
dev.off()

# STEP 2: Isotonic Regression on Constrained Ages #
### NOTE: Sample 195SQ10 outlier, not taken into account

# Define a matrix with stratigraphic constraints
SC2 = matrix(data=0,ncol=13,nrow=14)
SC2[1,]=rep(1,13)
SC2[2,]=c(0,0,0,0,1,1,1,1,1,1,1,1,1) # 109SQ10
SC2[3,]=c(0,0,0,0,1,1,1,1,1,1,1,1,1) # 130SQ10
SC2[4,]=c(0,0,0,0,1,1,1,1,1,1,1,1,1) # 170SQ10
SC2[5,]=c(0,0,0,0,1,1,1,1,1,1,1,1,1) # 173SQ10
SC2[6,]=c(0,0,0,0,0,0,0,0,0,0,0,1,1) # HRS5-Q
SC2[7,]=c(0,0,0,0,0,0,0,0,0,0,0,1,1) # HRS5-F
SC2[8,]=c(0,0,0,0,0,0,0,0,0,0,0,1,1) # 376SQ10
SC2[9,]=c(0,0,0,0,0,0,0,0,0,0,0,1,1) # 388SQ10
SC2[10,]=c(0,0,0,0,0,0,0,0,0,0,0,1,1) # 415SQ10
SC2[11,]=c(0,0,0,0,0,0,0,0,0,0,0,1,1) # 445SQ10
SC2[12,]=c(0,0,0,0,0,0,0,0,0,0,0,1,1) # 455SQ10
SC2[13,]=c(0,0,0,0,0,0,0,0,0,0,0,0,0) # HRS6-Q
SC2[14,]=c(0,0,0,0,0,0,0,0,0,0,0,0,0) # HRS6-F

# Isotonic Regression taking into account the Strati Matrix
Cramb_IsotonicDistorsion = IsotonicCurve(StratiConstraints = SC2, Holley_CP_14C_OSL,interactive = FALSE)

# Save the Strati matrix as a .pdf figure
p_strati <- recordPlot()
pdf("Holley_CP_StratiMatrix.pdf", width = 10, height = 8)
replayPlot(p_strati)
dev.off()

# Visualize the computed ages
knitr::kable(Cramb_IsotonicDistorsion$Ages )

# Plot and save the Isotonic Regression ages as a .pdf
AgeWith_IR <-  plot_Ages(object = Cramb_IsotonicDistorsion, plot_mode = "density")
p <- recordPlot(AgeWith_IR)
pdf("Holley_CP_Ages_IR.pdf", width = 10, height = 8)
replayPlot(p)
dev.off()

# Plot the Isotonic Curve
Holley_CrambPlots = PlotIsotonicCurve(object = Cramb_IsotonicDistorsion, level = .95)

# Save the Harris diagram as a .png figure
Holley_CrambPlots$DAG
ggsave("Holley_CP_DAG.png", 
       plot = Holley_CrambPlots$DAG, 
       width = 10, height = 8,
       dpi = 300)

# Save the Isotonic curve as a .png figure
Holley_CrambPlots$ribbon
ggsave("Holley_CP_IsotonicCurve.png", 
       plot = Holley_CrambPlots$ribbon, 
       width = 10, height = 8,
       dpi = 300)

# Save ages in a .csv file
write.csv(Cramb_IsotonicDistorsion$Ages,
          file = "Holley_Cramb_C14-OSL_IsotonicRegression.csv",
          row.names = FALSE)

# Save MCMC values in a .csv file
write.csv(as.matrix(Cramb_IsotonicDistorsion[["Sampling"]][[1]]),
          file = "Holley_Cramb_IsotonicRegression_MCMC.csv",
          row.names = FALSE)
