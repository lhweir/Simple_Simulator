#############################################################################
#   CN MSE Functions that create summary tables of performance statistics
#  Brooke Davis & Kendra Holt
#################################################################################

# ==================================================================================================
# Functions included in this file:
# ================================================================================================
# writePerfMeas_EscChange   - Calculates, summarizes, and writes csv with performance measures
#                             related to changes in escapement


# ==================================================================================================

writePerfMeas_EscChange<-function(Names, Run) {
  
  # want to extract table of median change in escapment over 10 years (2016 to 2025) and 20 years (2016 to 2035)
  # also save all blobs
  Blobs <- list()
  Escape <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    Escape[[mm]] <- Blobs[[mm]]$Sims$Escape
  }
  
  Esc_Init <- Blobs[[1]]$Data$Escapement
  
  # For each stock, model, get median escapements 
  # also get sims of 2022:2025 and 2032:2035
  Stocks <- Blobs[[1]]$Data$Stocks
  NY <- Blobs[[1]]$Data$NY
  Years <- Blobs[[1]]$Options$Years
  Esc_Means <- list()
  Esc_SDs <- list()
  Esc_22_25 <- list()
  Esc_32_35 <- list()
  EscChange_2025 <- list()
  EscChange_2035 <- list()
  AvgCatch <- list()
  AvgCatch_Init <- list()
  for(mm in 1:length(Names)){
    Esc_Means[[mm]] <- list()
    Esc_SDs[[mm]] <- list()
    Esc_22_25[[mm]] <- list()
    Esc_32_35[[mm]] <- list()
    EscChange_2025[[mm]] <- list()
    EscChange_2035[[mm]] <- list()
    nSims <- Blobs[[mm]]$Options$nSims
    
    # Get catch organized
    Term_Catch <- Blobs[[mm]]$Sims$Term_Catch
    Term_Yearly <- lapply(1:nSims, function(x) apply(Term_Catch[[x]], 1, sum))
    PT_Catch <- Blobs[[mm]]$Sims$PT_Catch
    PT_Yearly <-  lapply(1:nSims, function(x) apply(PT_Catch[[x]], 1, sum))
    # Now add together
    Total_Catch <- lapply(1:nSims, function(x) PT_Yearly[[x]] + Term_Yearly[[x]])
    # For each sim get avg of first 5 years, average of all years (minus initialization years)
    AvgCatch[[mm]] <- unlist(lapply(1:nSims, function(x) mean(Total_Catch[[x]][6:55])))
    AvgCatch_Init[[mm]] <- unlist(lapply(1:nSims, function(x) mean(Total_Catch[[x]][6:10])))
    
    for(ss in 1:length(Stocks)){
      
      # Get initial avg escapement from data inputs
      Esc_Start <- mean(Esc_Init$Abundance[which(Esc_Init$StockID==Stocks[ss] & Esc_Init$Year %in% 2012:2015)])
      StockDat <- lapply(Escape[[mm]], "[", ,ss , ) 
      # sum over stocks
      StockDat_Summed <- lapply(1:length(StockDat), function(x) apply(StockDat[[x]], 1, sum) )
      # sum over ages
      Esc_Means[[mm]][[ss]] <-  sapply(1:NY, function(x) mean(sapply(StockDat_Summed, "[", x) ))
      Esc_SDs[[mm]][[ss]] <- sapply(1:NY, function(x) sd(sapply(StockDat_Summed, "[", x) ))
      # extract 100 sims of averages
      Esc_22_25[[mm]][[ss]] <-  sapply(1:nSims, function(x) mean((StockDat_Summed[[x]][which(Years %in% 2022:2025)])) )
      Esc_32_35[[mm]][[ss]]<-  sapply(1:nSims, function(x) mean((StockDat_Summed[[x]][which(Years %in% 2032:2035)])) )
      
      EscChange_2025[[mm]][[ss]]<-Esc_22_25[[mm]][[ss]]/Esc_Start
      EscChange_2035[[mm]][[ss]]<-Esc_32_35[[mm]][[ss]]/Esc_Start
      
    } # end stock loop
  } # end mod loop
  
  # Make data frame summarizing differences
  ChangeDF <- data.frame(Stock=character(), Model=character(), EscChange_10y=numeric(), EscChange_20y=numeric(),
                         Pr2025=numeric(), Pr3035=numeric())
  
  for(ss in 1:length(Stocks)){
    # Get initial avg escapement from data inputs
    Esc_Start <- Esc_Init$Abundance[which(Esc_Init$StockID==Stocks[ss] & Esc_Init$Year %in% 2012:2015)]
    for(mm in 1:length(Names)){
      # number of sims where this is greater than init
      NewRow <- data.frame(Stock=Stocks[ss], Model=Names[mm],
                           EscChange_10y = median (EscChange_2025[[mm]][[ss]]),
                           EscChange_20y = median (EscChange_2035[[mm]][[ss]]),
                           Pr2025=sum(Esc_22_25[[mm]][[ss]] > mean(Esc_Start))/nSims,
                           Pr3035=sum(Esc_32_35[[mm]][[ss]] > mean(Esc_Start))/nSims)
      ChangeDF <- rbind(ChangeDF, NewRow)
    } # end model loop
  } # end stock loop
  
  write.csv(ChangeDF,paste("DataOut/Results_",Run ,".csv", sep=""))
  
}

