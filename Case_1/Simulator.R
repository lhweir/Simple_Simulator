####################################################
###       Basic Simulator Model                   ##
###       Kendra Holt & Brooke Davis              ##
####################################################

# Clear environment
rm(list=ls())

# Load required packages
library(dplyr)

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

setwd("Case_1")

# Create DataOut directory if it doesn't yet exist
outputDir <- "DataOut"
if (file.exists(outputDir) == FALSE){
  dir.create(outputDir)
}

# Create Figure directory if it doesn't yet exist
figuresDir <- "Figures"
if (file.exists(figuresDir) == FALSE){
  dir.create(figuresDir)
}


#**********************************************************************
# Scenarios for April meeting

Prod<-c("Basic_Ricker", "Recursive_Bayes_5yr")
Fish <- c("Current", "Fish50", "PT_Fish50", "NoFishing")
Names <- NULL
Scenario <- 1
nSims <- 2
Run <- "ER_Review_Test"

# Loop over OMs (basic Ricker vs recursive bayes)
for(i in 1:length(Prod)){
  # Loop over fishing scenarios
  for(j in 1:length(Fish)){
    # Store name for this scenario
    Names[Scenario] <- paste(Scenario, Prod[i], Fish[j], Run, sep="_")
    # Initialize "blob"
    Blob <- Init.Blob(SR = Prod[i], Years = 2011:2065, HRS = Fish[j], nSims = nSims,
                  Name=Names[Scenario],  BM="Basic_Ricker", Initialization ="Escapement",EV_Type="1",
                  Prod_Scenario=1, CC_Scenario=1, Smolts_Scenario = 1 ,  Exclude_Jacks=T)
    # Load data to blob
    Blob$Data <- Read.Data(Blob)
    # Run simulations
    Blob$Sims <- Run.CN.MSE.Sim(Blob)
    # Save blob
    Save.Blob(Blob)
    # increment scenario
    Scenario <- Scenario+1
  } # end fishery scenario loop

  # Now scenario with increased productivity
  Names[Scenario] <-  paste(Scenario, Prod[i], "Prod25", Run, sep="_")
  Blob <- Init.Blob(SR = Prod[i], Years = 2011:2065, HRS = "Current", nSims = nSims,
            Name=Names[Scenario],  BM="Basic_Ricker", Initialization ="Escapement",EV_Type="1",
            Prod_Scenario="25p_10yrs", CC_Scenario=1, Smolts_Scenario = 1 ,  Exclude_Jacks=T)

  Blob$Data <- Read.Data(Blob)
  Blob$Sims <- Run.CN.MSE.Sim(Blob)
  Save.Blob(Blob)
  Scenario <- Scenario+1

  # Now increase CC
  Names[Scenario] <-  paste(Scenario, Prod[i], "CC25", Run, sep="_")
  Blob <- Init.Blob(SR = Prod[i], Years = 2011:2065, HRS = "Current", nSims = nSims,
                    Name=Names[Scenario], BM="Basic_Ricker", Initialization ="Escapement",EV_Type="1",
                    Prod_Scenario=1, CC_Scenario="25p_10yrs", Smolts_Scenario = 1 , Exclude_Jacks=T)
  Blob$Data <- Read.Data(Blob)
  Blob$Sims <- Run.CN.MSE.Sim(Blob)
  Save.Blob(Blob)
  Scenario <- Scenario+1

} # end OM loop

# save list of names for this run
saveRDS(Names, paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))


####################################################
##  Call summary & plot functions
######################################################

# Save performance measures related to escapement change (writes csv called "Results")
writePerfMeas_EscChange(Names, Run)

# Plot 4-panel escapement
# Read in escapement "lead-in" data
Esc_LeadIn <- read.csv("DataIn/Escape_LeadIn.csv")
# Load names if not already loaded
Names <- readRDS(paste("DataOut/Scenario_Names_", Run, ".RDS", sep=""))
# call plot function
Escape_4_Panel(Names=Names, Esc_LeadIn, PlotName=paste("Escape_4Panel_", Run , sep=""), Plot_Spawners = T, Smax=F)

# Catch 4-panel plot for PT and terminal
Catch_4_Panel(Names=Names, PlotName="Catch_4Panel_PT", Mods2Plot_1 = c(1,2,3,4), Mods2Plot_2 = c(1,5,6), Mods2Plot_3 = c(7,8,9,10),
              Mods2Plot_4 = c(7,11,12), Scenarios = "Fishing_Habitat", Fishery="PT")


Catch_4_Panel(Names=Names, PlotName="Catch_4Panel_Term", Mods2Plot_1 = c(1,2,3,4), Mods2Plot_2 = c(1,5,6), Mods2Plot_3 = c(7,8,9,10),
              Mods2Plot_4 = c(7,11,12), Scenarios = "Fishing_Habitat", Fishery="Term")


# Dot and line plots showing average catch over the time series for all scenarios
plotCatch(Names, PlotName = "Catch_Barplots")


# Can also plot individual cohort abundance and escapement trajectories
# can give these functions 1 scenario, or multiple, but multiple gets messy fast
# don't need to give legend names, but names are usually messy

# Plot cohort abundance by age and stock for two sample runs
CohortPlots(Names[1:2], PlotName = "Sample_Cohort_Abundance",
            Legend_Names = c('Base Curr. Fishing', "Base 50% Fishing"))

# Same for escapement
EscapePlots(Names[1], PlotName = "Escapement_Example")






