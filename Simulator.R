####################################################
###       Basic Simulator Model                   ##
###       Kendra Holt & Brooke Davis              ##
###         Edited by Lauren Weir                 ##
####################################################

# Clear environment
rm(list=ls())

#setwd("C:/Users/WeirL/Documents/GitHub/Simple_Simulator")

# Load required packages
library(dplyr)
library(rlist)

# Source functions, which are in base folder
FuncPath <- paste("R/Functions.R", sep="/") # main functions for simulation routine
pFuncPath <- paste("R/plottingFunctions.R", sep="/") # functions to make plots
hFuncPath <- paste("R/helperFunctions.R", sep="/") # helper functions that are called from the simulation routine
sFuncPath <- paste( "R/setupFunctions.R", sep="/") # functions that help set-up the simulation routine (initialze arrays, read-in data)
statsFuncPath <- paste( "R/statsFunctions.R", sep="/") # functions that help set-up the simulation routine (initialze arrays, read-in data)


source(FuncPath)
source(pFuncPath)
source(hFuncPath)
source(sFuncPath)
source(statsFuncPath)



#################################################################
# DU 2 Harrison
#################################################################
#setwd("C:/Users/WeirL/Documents/GitHub/Simple_Simulator/DU2")

setwd("C:/gitHub/Simple_Simulator/DU2")

# Create DataOut directory if it doesn't yet exist
outputDir <- paste0("DataOut")
if (file.exists(outputDir) == FALSE){
  dir.create(outputDir)
}

# Create Figure directory if it doesn't yet exist
figuresDir <- paste0("Figures")
if (file.exists(figuresDir) == FALSE){
  dir.create(figuresDir)
}
########################

#======================================================================================
# Temporary code for testing posterior sampling

Prod<-"TimeVarying"
Fish <- "Fish0"
Names <- NULL
nSims <- 10000
Depensatory <- "Dep"
P_Scenario <- "0p_12yrs"

# Scenario 1: Sample from posterior estimated without LN Bias Correction; Projections without bias correction
Scenario <- 1
Run <- "DU2_postNoBC_projNoBC"
Names[Scenario] <- paste(Scenario, Prod, Fish, P_Scenario, Run, sep="_")

BlobHar <- Init.Blob(SR = Prod, Years = 2015:2031, HRS = Fish, nSims = nSims,
                     Name=Names[Scenario],  BM="WSP", Initialization ="Escapement", EV_Type="1", Prod_Scenario = P_Scenario,
                     CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0, AutoCorr=F,
                     LNormBias=F, samplePosterior = T, inputLNormBias=F)

# Load data to blob
BlobHar$Data <- Read.Data(BlobHar,FolderPath="DataIn")
# Run simulations
BlobHar$Sims <- Run.CN.MSE.Sim(BlobHar)
# Save blob
Save.Blob(BlobHar)

# Scenario 2: Sample from posterior estimated with LN Bias Correction; Projections with bias correction
Scenario <- 2
Run <- "DU2_postBC_projBC"
Names[Scenario] <- paste(Scenario, Prod, Fish, P_Scenario, Run, sep="_")

BlobHar <- Init.Blob(SR = Prod, Years = 2015:2031, HRS = Fish, nSims = nSims,
                     Name=Names[Scenario],  BM="WSP", Initialization ="Escapement", EV_Type="1", Prod_Scenario = P_Scenario,
                     CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0, AutoCorr=F,
                     LNormBias=T, samplePosterior = T, inputLNormBias=T)

# Load data to blob
BlobHar$Data <- Read.Data(BlobHar,FolderPath="DataIn")
# Run simulations
BlobHar$Sims <- Run.CN.MSE.Sim(BlobHar)
# Save blob
Save.Blob(BlobHar)

# # Scenario 3: Sample from posterior estimated without LN Bias Correction; Projections with bias correction
 Scenario <- 3
Run <- "DU2_postNoBC_projBC"
Names[Scenario] <- paste(Scenario, Prod, Fish, P_Scenario, Run, sep="_")

BlobHar <- Init.Blob(SR = Prod, Years = 2015:2031, HRS = Fish, nSims = nSims,
                     Name=Names[Scenario],  BM="WSP", Initialization ="Escapement", EV_Type="1", Prod_Scenario = P_Scenario,
                     CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0, AutoCorr=F,
                     LNormBias=T, samplePosterior = T, inputLNormBias=F)


# Load data to blob
BlobHar$Data <- Read.Data(BlobHar,FolderPath="DataIn")
# Run simulations
BlobHar$Sims <- Run.CN.MSE.Sim(BlobHar)
# Save blob
Save.Blob(BlobHar)



# Scenario 4: Use median estimated without LN Bias Correction; Projections without bias correction
Scenario <- 4
Run <- "DU2_medNoBC_projNoBC"
Names[Scenario] <- paste(Scenario, Prod, Fish, P_Scenario, Run, sep="_")

BlobHar <- Init.Blob(SR = Prod, Years = 2015:2031, HRS = Fish, nSims = nSims,
                     Name=Names[Scenario],  BM="WSP", Initialization ="Escapement", EV_Type="1", Prod_Scenario = P_Scenario,
                     CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0, AutoCorr=F,
                     LNormBias=F, samplePosterior = F, inputLNormBias=F)

# Load data to blob
BlobHar$Data <- Read.Data(BlobHar,FolderPath="DataIn")
# Run simulations
BlobHar$Sims <- Run.CN.MSE.Sim(BlobHar)
# Save blob
Save.Blob(BlobHar)

# Scenario 5: Use median estimated with LN Bias Correction; Projections with bias correction
Scenario <- 5
Run <- "DU2_medBC_projBC"
Names[Scenario] <- paste(Scenario, Prod, Fish, P_Scenario, Run, sep="_")

BlobHar <- Init.Blob(SR = Prod, Years = 2015:2031, HRS = Fish, nSims = nSims,
                     Name=Names[Scenario],  BM="WSP", Initialization ="Escapement", EV_Type="1", Prod_Scenario = P_Scenario,
                     CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0, AutoCorr=F,
                     LNormBias=T, samplePosterior = F, inputLNormBias=T)

# Load data to blob
BlobHar$Data <- Read.Data(BlobHar,FolderPath="DataIn")
# Run simulations
BlobHar$Sims <- Run.CN.MSE.Sim(BlobHar)
# Save blob
Save.Blob(BlobHar)

# # Scenario 6: Use medians estimated without LN Bias Correction; Projections with bias correction
Scenario <- 6
Run <- "DU2_medNoBC_projBC"
Names[Scenario] <- paste(Scenario, Prod, Fish, P_Scenario, Run, sep="_")

BlobHar <- Init.Blob(SR = Prod, Years = 2015:2031, HRS = Fish, nSims = nSims,
                     Name=Names[Scenario],  BM="WSP", Initialization ="Escapement", EV_Type="1", Prod_Scenario = P_Scenario,
                     CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0, AutoCorr=F,
                     LNormBias=T, samplePosterior = F, inputLNormBias=F)


# Load data to blob
BlobHar$Data <- Read.Data(BlobHar,FolderPath="DataIn")
# Run simulations
BlobHar$Sims <- Run.CN.MSE.Sim(BlobHar)
# Save blob
Save.Blob(BlobHar)




# Compare scenarios
EscapePlotsPercentiles_Compare(c(Names[3],Names[2], Names[1]), PlotName="Harrison - Bias correction comparison - posterior - percent", 
                               Stock = "Harrison", Legend_Names = c("PostNoBC_projBC", "PostBC_projBC","PostNoBC_projNoBC"))
EscapePlots_Compare(c(Names[3],Names[2],Names[1]), PlotName="Harrison - Bias correction comparison - posterior", 
                    Stock = "Harrison", Legend_Names = c("PostNoBC_projBC", "PostBC_projBC","PostNoBC_projNoBC"))

EscapePlotsPercentiles_Compare(c(Names[6],Names[5], Names[4]), PlotName="Harrison - Bias correction comparison - median - percent", 
                               Stock = "Harrison", Legend_Names = c("MedNoBC_projBC", "MedBC_projBC","MedNoBC_projNoBC"))
EscapePlots_Compare(c(Names[6],Names[5],Names[4]), PlotName="Harrison - Bias correction comparison - median", 
                    Stock = "Harrison", Legend_Names = c("MedNoBC_projBC", "MedBC_projBC","MedNoBC_projNoBC"))

#checkNSims (Names, perfMetric = "EscapeChange") 

# End of temporary code for posterior method  ===========================================



#####################################################

Prod<-c("TimeVarying")#AutoCorr, Simple, TimeVarying, Survival
Fish <- c("Fish0")#, "Fish+10", "Fish10","Fish20","Fish30","Fish40","Fish50","Fish60","Fish70","Fish80","Fish90","Fish100")
Names <- NULL
Scenario <- 1
nSims <- 100
Run <- "DU2"
Depensatory <- c("Dep")
P_Scenario <- c("0p_12yrs")#, "10p_12yrs", "-10p_12yrs", "20p_12yrs", "-20p_12yrs", "30p_12yrs", "-30p_12yrs", "40p_12yrs", "-40p_12yrs", "50p_12yrs", "-50p_12yrs")

# Loop over OMs (basic Ricker vs others)
for(i in 1:length(Prod)){
  # Loop over fishing scenarios
  for(j in 1:length(Fish)){
    for (p in 1:length(P_Scenario)){
    # Store name for this scenario
    Names[Scenario] <- paste(Scenario, Prod[i], Fish[j], P_Scenario[p], Run, sep="_")
    print(Names[Scenario])
    # Initialize "blob"
      BlobHar <- Init.Blob(SR = Prod[i], Years = 2015:2031, HRS = Fish[j], nSims = nSims,
                          Name=Names[Scenario],  BM="WSP", Initialization ="Escapement", EV_Type="1", Prod_Scenario = P_Scenario[p],
                          CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0, AutoCorr=F,
                          LNormBias=T, samplePosterior = T)
    
      # BlobHar <- Init.Blob(SR = Prod[i], Years = 2015:2031, HRS = Fish[j], nSims = nSims,
      #                      Name=Names[Scenario],  BM="WSP", Initialization ="Escapement", EV_Type="1", Prod_Scenario = P_Scenario[p],
      #                      CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0, AutoCorr=F,
      #                      LNormBias=F, samplePosterior = F)
    
    # Load data to blob
    BlobHar$Data <- Read.Data(BlobHar,FolderPath="DataIn")
    # Run simulations
    BlobHar$Sims <- Run.CN.MSE.Sim(BlobHar)
    # Save blob
    Save.Blob(BlobHar)
    # increment scenario
    Scenario <- Scenario + 1
    } # end productivity scenario loop
  } # end fishery scenario loop
  
} # end OM loop

# save list of names for this run
saveRDS(Names, paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

#Need to read in Escapment Lead in and Names for plotting
Esc_LeadIn <- read.csv("DataIn/Escape_LeadIn.csv")
Names <- readRDS(paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

#Can pick which Names you want
EscapePlots.Quant(Names[1], PlotName = "Scenario3_LNPar_withBCF", nr=1, nc=1) 
#EscapePlots.Quant(Names[1], PlotName = "Scenario1_noLN_withoutBCF", nr=1, nc=1)
#EscapePlots(Names[1], PlotName = "TV_Current_2020to2031Mean_withBCFs", nr=1, nc=1) 


Names <- c("Scenario1_OrgPar_noBCF","Scenario2_LNPar_withBCF","Scenario3_OrgPar_withBCF")
labs<- c("Scenario 1 - Original Analysis", "Scenario 2 - Correction Factor in both Est and Proj","Scenario 3 - Original Param with BCF")
EscapePlots_4Quant(Names, PlotName = "Harrison - Various BCF Scenarios", labels=labs, nr=1, nc=3)

savemeans(Names, Run)
savemedians(Names, Run)

RecoveryResults(Names, Run)
#RecoveryResultsgeoMean(Names, Run)

write_HR_atage_sim(Names, Run)
write_totalHR_sim(Names, Run)

# Excel files are too large to open properly
HR_atage <- read.csv("DataOut/HRatAge_SimResults_DU2.csv")
HR_total <- read.csv("DataOut/TotalCatchRate_SimResults_DU2.csv")
results<-read.csv("DataOut/Recovery_Results_DU2.csv")

#################################################################
# Examples for the other DUs with different productivity Values
#################################################################

#################################################################
# DU 8 Portage
#################################################################
setwd("C:/Users/WeirL/Documents/GitHub/Simple_Simulator/DU8")


# Create DataOut directory if it doesn't yet exist
outputDir <- paste0("DataOut")
if (file.exists(outputDir) == FALSE){
  dir.create(outputDir)
}

# Create Figure directory if it doesn't yet exist
figuresDir <- paste0("Figures")
if (file.exists(figuresDir) == FALSE){
  dir.create(figuresDir)
}
########################

Prod<-c(0,50,55,60,65,70,75,80,85,90)#, "Recursive_Bayes_5yr" ) # Add in , "Time_Varying" if being used 
Fish <- c("Fish0")#, "Fish25_5yrs", "Fish50_5yrs","Fish75_5yrs")#, "PT_Fish50", "NoFishing")
Names <- NULL
Scenario <- 1
nSims <- 10000
Run <- "DU8"
Depensatory <- c("Dep")

# Loop over OMs (basic Ricker vs others)
for(i in 1:length(Prod)){
  # Loop over fishing scenarios
  for(j in 1:length(Fish)){
    # Store name for this scenario
    Names[Scenario] <- paste(Scenario, Prod[i], Fish[j], Run, sep="_")
    # Initialize "blob"
    Blob8 <- Init.Blob(SR = Prod[i], Years = 2013:2033, HRS = Fish[j], nSims = nSims,
                         Name=Names[Scenario],  BM="WSP", Initialization ="Escapement",EV_Type="1",
                         Prod_Scenario = 1, CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0)
    # Load data to blob
    Blob8$Data <- Read.Data(Blob8)
    # Run simulations
    Blob8$Sims <- Run.CN.MSE.Sim(Blob8)
    # Save blob
    Save.Blob(Blob8) 
    # increment scenario
    Scenario <- Scenario + 1
  } # end fishery scenario loop
  
} # end OM loop

# Save performance measures related to escapement change (writes csv called "Results")

# save list of names for this run
saveRDS(Names, paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

#Need to read in Escapment Lead in and Names for plotting
Esc_LeadIn <- read.csv("DataIn/Escape_LeadIn.csv")
Names <- readRDS(paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

EscapePlots_4Quant(Names[2:5], PlotName = "Current_DU8_50-65alphaQuant", alpha=c(50,55,60, 65),nr=2, nc=2)
EscapePlots_4Quant(Names[6:9], PlotName = "Current_DU8_70-85alphaQuant", alpha=c(70,75,80, 85),nr=2, nc=2)

savemeans(Names, Run)

RecoveryResults(Names, Run)

write_HR_atage_sim(Names, Run)
write_totalHR_sim(Names, Run)


#################################################################
# DU 9 MFR-Springs
#################################################################
setwd("C:/Users/WeirL/Documents/GitHub/Simple_Simulator/DU9")


# Create DataOut directory if it doesn't yet exist
outputDir <- paste0("DataOut")
if (file.exists(outputDir) == FALSE){
  dir.create(outputDir)
}

# Create Figure directory if it doesn't yet exist
figuresDir <- paste0("Figures")
if (file.exists(figuresDir) == FALSE){
  dir.create(figuresDir)
}
########################

Prod<-c(0,50,55,60,65,70,75,80,85,90) 
Fish <- c("Fish0")
Names <- NULL
Scenario <- 1
nSims <- 10000
Run <- "DU9"
Depensatory <- c("Dep")
BigB <- c(1)

# Loop over OMs (basic Ricker vs others)
for(i in 1:length(Prod)){
  # Loop over fishing scenarios
  for(j in 1:length(Fish)){
    
    for(bb in 1:length(BigB)){
    # Store name for this scenario
    Names[Scenario] <- paste(Scenario, Prod[i], Fish[j], Run, bb, sep="_")
    # Initialize "blob"
    Blob9 <- Init.Blob(SR = Prod[i], Years = 2013:2033, HRS = Fish[j], nSims = nSims,
                       Name=Names[Scenario],  BM="WSP", Initialization ="Escapement",EV_Type="1",
                       Prod_Scenario = 1, CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=bb)
    # Load data to blob
    Blob9$Data <- Read.Data(Blob9)
    # Run simulations
    Blob9$Sims <- Run.CN.MSE.Sim(Blob9)
    # Save blob
    Save.Blob(Blob9) 
    # increment scenario
    Scenario <- Scenario + 1
    } # End Big Bar loop
  } # end fishery scenario loop
  
} # end OM loop


# save list of names for this run
saveRDS(Names, paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

#Need to read in Escapment Lead in and Names for plotting
Esc_LeadIn <- read.csv("DataIn/Escape_LeadIn.csv")
Names <- readRDS(paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

EscapePlots_4Quant(Names[2:5], PlotName = "Current_DU9_50-65alpha_BB1QT", alpha=c(50,55,60, 65),nr=2, nc=2)
EscapePlots_4Quant(Names[6:9], PlotName = "Current_DU9_70-85alpha_BB1QT", alpha=c(70,75,80, 85),nr=2, nc=2)

savemeans(Names, Run)

RecoveryResults(Names, Run)

write_HR_atage_sim(Names, Run)
write_totalHR_sim(Names, Run)



#################################################################
# DU 10 MFR-Summers
#################################################################
setwd("C:/Users/WeirL/Documents/GitHub/Simple_Simulator/DU10")


# Create DataOut directory if it doesn't yet exist
outputDir <- paste0("DataOut")
if (file.exists(outputDir) == FALSE){
  dir.create(outputDir)
}

# Create Figure directory if it doesn't yet exist
figuresDir <- paste0("Figures")
if (file.exists(figuresDir) == FALSE){
  dir.create(figuresDir)
}
########################

Prod<-c(0,50,55,60,65,70,75,80,85,90)
Fish <- c("Fish0")
Names <- NULL
Scenario <- 1
nSims <- 10000
Run <- "DU10"
Depensatory <- c("Dep")
BigB <- c(1)

# Loop over OMs (basic Ricker vs others)
for(i in 1:length(Prod)){
  # Loop over fishing scenarios
  for(j in 1:length(Fish)){
    for(bb in 1:length(BigB)){
    # Store name for this scenario
    Names[Scenario] <- paste(Scenario, Prod[i], Fish[j], Run,bb, sep="_")
    # Initialize "blob"
    Blob10 <- Init.Blob(SR = Prod[i], Years = 2013:2033, HRS = Fish[j], nSims = nSims,
                       Name=Names[Scenario],  BM="WSP", Initialization ="Escapement",EV_Type="1",
                       Prod_Scenario = 1, CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=bb)
    # Load data to blob
    Blob10$Data <- Read.Data(Blob10)
    # Run simulations
    Blob10$Sims <- Run.CN.MSE.Sim(Blob10)
    # Save blob
    Save.Blob(Blob10) 
    # increment scenario
    Scenario <- Scenario + 1
    } # End BigBAr loop
  } # end fishery scenario loop
  
} # end OM loop

# save list of names for this run
saveRDS(Names, paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

#Need to read in Escapment Lead in and Names for plotting
Esc_LeadIn <- read.csv("DataIn/Escape_LeadIn.csv")
Names <- readRDS(paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

EscapePlots_4Quant(Names[2:5], PlotName = "Current_DU10_50-65alpha_BB1QT", alpha=c(50,55,60, 65),nr=2, nc=2)
EscapePlots_4Quant(Names[6:9], PlotName = "Current_DU10_70-85alpha_BB1QT", alpha=c(70,75,80, 85),nr=2, nc=2)

savemeans(Names, Run)

RecoveryResults(Names, Run)

write_HR_atage_sim(Names, Run)
write_totalHR_sim(Names, Run)


#################################################################
# DU 11 UFR-Springs
#################################################################
setwd("C:/Users/WeirL/Documents/GitHub/Simple_Simulator/DU11")


# Create DataOut directory if it doesn't yet exist
outputDir <- paste0("DataOut")
if (file.exists(outputDir) == FALSE){
  dir.create(outputDir)
}

# Create Figure directory if it doesn't yet exist
figuresDir <- paste0("Figures")
if (file.exists(figuresDir) == FALSE){
  dir.create(figuresDir)
}
########################

Prod<-c(0,50,55,60,65,70,75,80,85,90)
Fish <- c("Fish0")
Names <- NULL
Scenario <- 1
nSims <- 10000
Run <- "DU11"
Depensatory <- c("Dep")
BigB <- c(1)

# Loop over OMs (basic Ricker vs others)
for(i in 1:length(Prod)){
  # Loop over fishing scenarios
  for(j in 1:length(Fish)){
    
    for(bb in 1:length(BigB)){
      # Store name for this scenario
      Names[Scenario] <- paste(Scenario, Prod[i], Fish[j], Run, bb, sep="_")
      # Initialize "blob"
      Blob11 <- Init.Blob(SR = Prod[i], Years = 2013:2033, HRS = Fish[j], nSims = nSims,
                         Name=Names[Scenario],  BM="WSP", Initialization ="Escapement",EV_Type="1",
                         Prod_Scenario = 1, CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=bb)
      # Load data to blob
      Blob11$Data <- Read.Data(Blob11)
      # Run simulations
      Blob11$Sims <- Run.CN.MSE.Sim(Blob11)
      # Save blob
      Save.Blob(Blob11) 
      # increment scenario
      Scenario <- Scenario + 1
    } # End Big Bar loop
  } # end fishery scenario loop
  
} # end OM loop

# save list of names for this run
saveRDS(Names, paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

#Need to read in Escapment Lead in and Names for plotting
Esc_LeadIn <- read.csv("DataIn/Escape_LeadIn.csv")
Names <- readRDS(paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

EscapePlots_4Quant(Names[2:5], PlotName = "Current_DU11_50-65alpha_BB1QT", alpha=c(50,55,60, 65),nr=2, nc=2)
EscapePlots_4Quant(Names[6:9], PlotName = "Current_DU11_70-85alpha_BB1QT", alpha=c(70,75,80, 85),nr=2, nc=2)

savemeans(Names, Run)

RecoveryResults(Names, Run)

write_HR_atage_sim(Names, Run)
write_totalHR_sim(Names, Run)



#################################################################
# DU 17 NTh-Summers
#################################################################
setwd("C:/Users/WeirL/Documents/GitHub/Simple_Simulator/DU17")


# Create DataOut directory if it doesn't yet exist
outputDir <- paste0("DataOut")
if (file.exists(outputDir) == FALSE){
  dir.create(outputDir)
}

# Create Figure directory if it doesn't yet exist
figuresDir <- paste0("Figures")
if (file.exists(figuresDir) == FALSE){
  dir.create(figuresDir)
}
########################

Prod<-c(0,50,55,60,65,70,75,80,85,90) 
Fish <- c("Fish0")
Names <- NULL
Scenario <- 1
nSims <- 10000
Run <- "DU17"
Depensatory <- c("Dep")

# Loop over OMs (basic Ricker vs others)
for(i in 1:length(Prod)){
  # Loop over fishing scenarios
  for(j in 1:length(Fish)){
    # Store name for this scenario
    Names[Scenario] <- paste(Scenario, Prod[i], Fish[j], Run, sep="_")
    # Initialize "blob"
    Blob17 <- Init.Blob(SR = Prod[i], Years = 2013:2033, HRS = Fish[j], nSims = nSims,
                       Name=Names[Scenario],  BM="WSP", Initialization ="Escapement",EV_Type="1",
                       Prod_Scenario = 1, CC_Scenario = 1, Smolts_Scenario = 1 ,  Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0)
    # Load data to blob
    Blob17$Data <- Read.Data(Blob17)
    # Run simulations
    Blob17$Sims <- Run.CN.MSE.Sim(Blob17)
    # Save blob
    Save.Blob(Blob17) 
    # increment scenario
    Scenario <- Scenario + 1
  } # end fishery scenario loop
  
} # end OM loop

# save list of names for this run
saveRDS(Names, paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

#Need to read in Escapment Lead in and Names for plotting
Esc_LeadIn <- read.csv("DataIn/Escape_LeadIn.csv")
Names <- readRDS(paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))

EscapePlots_4Quant(Names[2:5], PlotName = "Current_DU17_50-65alphaQT", alpha=c(50,55,60, 65),nr=2, nc=2)
EscapePlots_4Quant(Names[6:9], PlotName = "Current_DU17_70-85alphaQT", alpha=c(70,75,80, 85),nr=2, nc=2)

savemeans(Names, Run)

RecoveryResults(Names, Run)

write_HR_atage_sim(Names, Run)
write_totalHR_sim(Names, Run)




