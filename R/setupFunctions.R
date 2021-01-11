#########################################################
##     Helper Functions for Chinook MSE Simulator      ##
##            Brooke Davis & Kendra Holt               ##
##              Edited by Lauren Weir                  ##
#########################################################


# ========================================================================
# Functions included in this file:
# ========================================================================
# Init.Blob
# Read.Data
# Save.Blob
# ========================================================================


#######################################
#  Init.Blob
######################################
# Purpose: Initiatizes a list to save simulation results for a single scenario
# Arguments: Name = scenario name
#           SR = type of SR model fit used to model productivity ("Basic_Ricker" or "Recursive_Bayes")
#           BM = type of SR model fit used to calculate benchmarks ("Basic_Ricker" or "Recursive_Bayes")
#           Initialization = method used to initialize abundance ("Escapement" or "Cohort")
#           EV_Type = method used to calculate annual EV recruitment scalars.  If 1, EV scalars are not used.
#           HRS = scenario used to specify annual harvest rate scalars ("Current", "Fish50", "PT_Fish50", "NoFishing")
#           nSims = number of simulation replicates
#           minVuln = minimum vulnerability (should move later ...)
#           Prod_Scenario = Productivity scenario (how alpha changes over time)
#           CC_Scenario = Carrying Capacity scenario (how Smax changes over time)
#           Depensatory_Effects = If TRUE, depensatory effects are included
#           BigBar = If Big Bar = 1, then effects from Big Bar are included
#           AutoCorr = If TRUE, depensatory effects are included
#           projectLNormBias = If TRUE, a lognormal bias correction is applied to projected recruitment 
#           inputLNormBias = If TRUE, the SR parameters input from csv file were estimated using a lognormal bias correction factor
###########################################

# Initialize simulation run with options
Init.Blob <- function(Name, SR, BM, Initialization, Years, HRS, EV_Type, nSims, Prod_Scenario,
                       CC_Scenario,  Smolts_Scenario, Exclude_Jacks=T, Depensatory_Effects=F, BigBar=0, 
                       AutoCorr, projectLNormBias, samplePosterior=F, inputLNormBias=F) {
  Blob <- list()
  Blob$Options <- list()
  # Add options
  Blob$Options$Name <- Name
  Blob$Options$SR <- SR
  Blob$Options$BM <- BM
  Blob$Options$Initialization <- Initialization
  Blob$Options$Years <- Years
  Blob$Options$HRS <- HRS
  Blob$Options$EV_Type = EV_Type
  Blob$Options$nSims <- nSims
  Blob$Options$Prod_Scenario <- Prod_Scenario
  Blob$Options$CC_Scenario <- CC_Scenario
  Blob$Options$Smolts_Scenario <- Smolts_Scenario
  Blob$Options$Exclude_Jacks <- Exclude_Jacks
  Blob$Options$Depensatory_Effects <- Depensatory_Effects
  Blob$Options$BigBar <- BigBar
  Blob$Options$AutoCorr <- AutoCorr
  Blob$Options$LNormBias <- projectLNormBias
  Blob$Options$samplePosterior <- samplePosterior
  Blob$Options$inputLNormBias <- inputLNormBias
  Blob
}



#######################################
#  Read.Data
######################################
# Purpose:  Reads in all csv data files
# Arguments: Blob = master scenario list to which data and model results are saved to
#            FolderPath = spaecifies name of folder that data input files are stored in
###########################################
Read.Data <- function(Blob, FolderPath="DataIn"){

  # Create list of Data to be output
  Data <- list()

  #*********************************
  # Add NY to Data
  Data$NY <- length(Blob$Options$Years)

  #**********************************
  # Get info about periods, regions
  Period_Tab <- read.csv(paste(FolderPath, "/Periods.csv", sep=""))
  # Number of periods
  Data$NP <- length(unique(Period_Tab$Period_Name))
  #Same for regions
  Data$Region_Tab <- read.csv(paste(FolderPath, "/Fishery_Regions.csv", sep=""))
  Data$NR <- length(unique(Data$Region_Tab$Region_Name))

  #**********************
  #Read in stock info
  StocksInfo <- read.csv(paste(FolderPath, "/Stocks_Info.csv", sep=""))
  Data$MaxAge <- max(StocksInfo$MaxAge)
  Data$NAges <- Data$MaxAge-1 # don't keep track of age 1's except for cohort
  Data$Ages <- 2:Data$MaxAge
  Data$Stocks <- unique(StocksInfo$StockID)
  Data$NS <- length(unique(Data$Stocks))

  # Add SR Params and BM's according to options
  if(Data$Stocks != "Harrison") {
  SR_File <- read.csv(paste(FolderPath, "/SR_Params/HabModel_SR.csv", sep=""))
  SR_Dat <- SR_File[SR_File$Reduction == Blob$Options$SR, ]}

  if(Data$Stocks == "Harrison") {
    if (Blob$Options$inputLNormBias == TRUE) {
        SR_Dat <- read.csv(paste(FolderPath, "/SR_Params/", Blob$Options$SR, "_LN_SR.csv", sep="")) }
    else {
        SR_Dat <- read.csv(paste(FolderPath, "/SR_Params/", Blob$Options$SR, "_SR.csv", sep=""))
    }
  }
  
  BM_Dat <- read.csv(paste(FolderPath, "/Benchmarks/", Blob$Options$BM, "_BM.csv", sep=""))

  # Convert to StockID fields to characters to avoid warning message:
  StocksInfo[,1]<-as.character(StocksInfo[,1])
  BM_Dat[,1]<-as.character(BM_Dat[,1])
  SR_Dat[,1]<-as.character(SR_Dat[,1])

  # Now add to StocksInfo file
  Data$StocksInfo <-  StocksInfo %>% left_join( SR_Dat, by="StockID") %>% left_join( BM_Dat, by="StockID")

  # Add joint posterior from SR estimation if full posterior is required
  if (Blob$Options$samplePosterior == TRUE) {
    if (Blob$Options$inputLNormBias == TRUE) {
      Data$SR_Post <- read.csv(paste(FolderPath, "/SR_Params/", Blob$Options$SR, "_LN_SR_fullPosterior.csv", sep=""))
    }
    else {
      Data$SR_Post <- read.csv(paste(FolderPath, "/SR_Params/", Blob$Options$SR, "_SR_fullPosterior.csv", sep=""))
    }
  }
  
  # Now read in multipliers for productivity/CC scenario and create data frame
  Prod_Mults <- data.frame(Year = Blob$Options$Year, Mult=1)
  CC_Mults <- data.frame(Year = Blob$Options$Year, Mult=1)
  # If not scenario where mult=1, fill in
  if(Blob$Options$Prod_Scenario != 1){
    Prod_DF <- read.csv(paste(FolderPath,"/Productivity_Change/",Blob$Options$Prod_Scenario, ".csv", sep=""))
    # Merge with Prod_Mults
    Prod_Mults$Mult[which(Prod_Mults$Year %in% Prod_DF$Year)] <- Prod_DF$Mult
  }
  if(Blob$Options$CC_Scenario != 1){
    CC_DF <- read.csv(paste(FolderPath,"/CC_Change/",Blob$Options$CC_Scenario, ".csv", sep=""))
    # Merge with CC_Mults
    CC_Mults$Mult[which(CC_Mults$Year %in% CC_DF$Year)] <- CC_DF$Mult
  }
  # add to list for output
  Data$Prod_Mults <- Prod_Mults
  Data$CC_Mults <- CC_Mults

  #*****************************
  # Read in fisheries info
  Data$FisheryInfo <- read.csv(paste(FolderPath, "/Fisheries.csv", sep=""))

  #Number of fisheries
  Data$NF <- length(unique(Data$FisheryInfo$FisheryID))

  #*******************************
  # Read in stock initialization data

  if(Blob$Options$Initialization == "Cohort"){
    # Base year cohort abundance by age
    Data$BaseCohort <- read.csv(paste(FolderPath, "/Starting_Cohort.csv", sep=""))
  } else if (Blob$Options$Initialization == "Escapement") {
    #**** read in escapement file insteaed
    Data$Escapement <- read.csv(paste(FolderPath, "/Escapement_Init.csv", sep=""))
  }

  #********************************
  # Exploitation Rates
  Data$Base_ER <- read.csv(paste(FolderPath, "/Effort_Based_Total_Mortality.csv", sep=""))

  #**********************************************
  #Regional Distribution coefficients
  Data$RDist <- read.csv(paste(FolderPath, "/Dist_Table.csv", sep=""))

  #*******************************************
  # Terminal harvest rate by stock, fishery, age
  Data$TermHR <- read.csv(paste(FolderPath, "/Stock_TermHR_Total_Mortality.csv", sep=""))

  #**********************************
  # Fishery harvest Rate scalars
  ##### Depends on fishing scenario option use
  # Not sure if this is best layout but will work for now? Maybe built in way to set HRS instead
  # of reading in from file?
  Data$HRS <- read.csv(paste(FolderPath, "/HRS/", Blob$Options$HRS, ".csv", sep=""))

  #*****************************
  # Group cohort survival
  Data$Grp_Surv <- read.csv(paste(FolderPath, "/Group_Survival.csv", sep=""))

  #********************************
  # Maturation by stock, age -- period in which they mature is in StocksInfo
  Data$Maturation <- read.csv(paste(FolderPath, "/Stock_MatRate.csv", sep=""))

  #********************************
  # Which file read in will depend on which type of EV we are using
  # Not sure if will go back to using EV's but will leave as option
  if(Blob$Options$EV_Type != 1){
    Data$EVs <- read.csv(paste(FolderPath, "/" , Blob$Options$EV_Type, ".csv", sep=""))
  }

  #**********************************
  #Smolts scenario
  if(Blob$Options$Smolts_Scenario != 1){
     Data$Smolts <- read.csv(paste(FolderPath, "/Smolt_Scenarios/", Blob$Options$Smolts_Scenario, ".csv", sep=""))
  }

  #*********************************
  # Set up the escapement survival if BigBar is TRUE
  if(Blob$Options$BigBar != 0){
    Data$BigBar <- read.csv(paste(FolderPath, "/BigBar/BigBar", Blob$Options$BigBar,".csv", sep=""))
  }
  
  # Return data list to be added to blob
  Data
} # end read.data function

###################
##  Save Blob    ##
###################

Save.Blob <- function(Blob){
  saveRDS(Blob, paste("DataOut/",Blob$Options$Name, ".rds", sep=""))
}

