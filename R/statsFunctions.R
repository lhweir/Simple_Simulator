#############################################################################
#   CN MSE Functions that create summary tables of performance statistics
#  Brooke Davis & Kendra Holt
#################################################################################

# ==================================================================================================
# Functions included in this file:
# ================================================================================================
# per.change                - Calculates percent change for the recovery results function
# writePerfMeas_EscChange   - Calculates, summarizes, and writes csv with performance measures
#                             related to changes in escapement
# write_HR_atage            - Calculates harvest rates at age for every year and simulation
# write_totalHR             - Calculates total HR for each year and simulation
# savemeans                 - Saves the mean for each model run
# RecoveryResults           - Save the p.change, abund, and evaluates whether the recovery targets are met. 
#                             Saves a 1 if each target is met 
# ==================================================================================================

# Calculates percent change
per.change<-function(data.in,na.rm=TRUE){
  
  # Log the data, but have to convert zeros to NA
  if(length(data.in$spn[data.in$spn == 0] != 0)) {data.in$spn[data.in$spn == 0] <-NA}
  data.in$spn <- log(data.in$spn)
  
  if(na.rm){data.in <- na.omit(data.in)}
  n<-dim(data.in)[[1]]-1
  lm<-lm(data.in[,"spn"]~data.in[,"yr"])
  pchange<-(exp(lm$coefficients[[2]]*n)-1)
  return(pchange)
}   

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
  Esc_27_30 <- list()
  #Esc_32_35 <- list()
  EscChange_2030 <- list()
  #EscChange_2035 <- list()
  AvgCatch <- list()
  AvgCatch_Init <- list()
  for(mm in 1:length(Names)){
    Esc_Means[[mm]] <- list()
    Esc_SDs[[mm]] <- list()
    Esc_27_30[[mm]] <- list()
    #Esc_32_35[[mm]] <- list()
    EscChange_2030[[mm]] <- list()
    #EscChange_2035[[mm]] <- list()
    nSims <- Blobs[[mm]]$Options$nSims
    
    # Get catch organized
    Term_Catch <- Blobs[[mm]]$Sims$Term_Catch
    Term_Yearly <- lapply(1:nSims, function(x) apply(Term_Catch[[x]], 1, sum))
    PT_Catch <- Blobs[[mm]]$Sims$PT_Catch
    PT_Yearly <-  lapply(1:nSims, function(x) apply(PT_Catch[[x]], 1, sum))
    # Now add together
    Total_Catch <- lapply(1:nSims, function(x) PT_Yearly[[x]] + Term_Yearly[[x]])
    # For each sim get avg of first 5 years, average of all years (minus initialization years)
    AvgCatch[[mm]] <- unlist(lapply(1:nSims, function(x) mean(Total_Catch[[x]][6:17])))
    AvgCatch_Init[[mm]] <- unlist(lapply(1:nSims, function(x) mean(Total_Catch[[x]][6:10])))
    
    for(ss in 1:length(Stocks)){
      
      # Get initial avg escapement from data inputs
      Esc_Start <- mean(Esc_Init$Abundance[which(Esc_Init$StockID==Stocks[ss] & Esc_Init$Year %in% 2014:2018)])
      StockDat <- lapply(Escape[[mm]], "[", ,ss , ) 
      # sum over stocks
      StockDat_Summed <- lapply(1:length(StockDat), function(x) apply(StockDat[[x]], 1, sum) )
      # sum over ages
      Esc_Means[[mm]][[ss]] <-  sapply(1:NY, function(x) mean(sapply(StockDat_Summed, "[", x) ))
      Esc_SDs[[mm]][[ss]] <- sapply(1:NY, function(x) sd(sapply(StockDat_Summed, "[", x) ))
      # extract 100 sims of averages
      Esc_27_30[[mm]][[ss]] <-  sapply(1:nSims, function(x) mean((StockDat_Summed[[x]][which(Years %in% 2027:2030)])) )
      #Esc_32_35[[mm]][[ss]]<-  sapply(1:nSims, function(x) mean((StockDat_Summed[[x]][which(Years %in% 2032:2035)])) )
      
      EscChange_2030[[mm]][[ss]]<-Esc_27_30[[mm]][[ss]]/Esc_Start
      #EscChange_2035[[mm]][[ss]]<-Esc_32_35[[mm]][[ss]]/Esc_Start
      
    } # end stock loop
  } # end mod loop
  
  # Make data frame summarizing differences
  ChangeDF <- data.frame(Stock=character(), Model=character(), EscChange_12y=numeric(), Pr2030=numeric())
  
  for(ss in 1:length(Stocks)){
    # Get initial avg escapement from data inputs
    Esc_Start <- Esc_Init$Abundance[which(Esc_Init$StockID==Stocks[ss] & Esc_Init$Year %in% 2012:2015)]
    for(mm in 1:length(Names)){
      # number of sims where this is greater than init
      NewRow <- data.frame(Stock=Stocks[ss], Model=Names[mm],
                           EscChange_12y = median (EscChange_2030[[mm]][[ss]]),
                           #EscChange_20y = median (EscChange_2035[[mm]][[ss]]),
                           Pr2030=sum(Esc_27_30[[mm]][[ss]] > mean(Esc_Start))/nSims)
                           #Pr3035=sum(Esc_32_35[[mm]][[ss]] > mean(Esc_Start))/nSims)
      ChangeDF <- rbind(ChangeDF, NewRow)
    } # end model loop
  } # end stock loop
  
  write.csv(ChangeDF,paste("DataOut/Results_",Run ,".csv", sep=""))
  
}


##########################################################
# Functions to calculate the HR in each year and each sim
##########################################################
write_HR_atage_year <-function(Names, Run) {
  
  # want to extract table of harvest rates used in the model
  # also save all blobs
  Blobs <- list()
  PT_Catch <- list()
  Term_Catch <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    PT_Catch[[mm]] <- Blobs[[mm]]$Sims$PT_Rates
    Term_Catch[[mm]] <- Blobs[[mm]]$Sims$Term_Rates
  }
  
  
  # For each stock, model, sim and Age get mean harvest rates 
  Stock <- Blobs[[1]]$Data$Stocks
  NY <- Blobs[[1]]$Data$NY
  Age <- Blobs[[1]]$Data$Ages
  Years <- Blobs[[1]]$Options$Years
  nSims <- Blobs[[1]]$Options$nSims
  
  HarvRates <- list()
  rr<-1
  
  for (mm in 1:length(Names)){
    
    print(Names[[mm]])
    
    for(ss in 1:nSims){
      
      # Make data frame summarizing differences
      HarvRatesDF <- data.frame(Stock=character(), Model=character(), Fish.Type = character(), Year = numeric(), Sim=numeric(), PT_Age2=numeric(),PT_Age3=numeric(),
                                PT_Age4=numeric(), PT_Age5=numeric(),T_Age2=numeric(),T_Age3=numeric(),T_Age4=numeric(), T_Age5=numeric())
      
      for (ff in 1:2){
        
        for (yy in 1:length(Years)) {
          
          ifelse(ff == 1, f.type <- "CA", f.type <- "US")
          
          if (ff == 1) {
            t2 = ifelse(Stock=="Harrison",Term_Catch[[mm]][[ss]][yy,,ff,1,1], Term_Catch[[mm]][[ss]][yy,,ff,1,2])
            t3 = ifelse(Stock=="Harrison",Term_Catch[[mm]][[ss]][yy,,ff,1,2], Term_Catch[[mm]][[ss]][yy,,ff,1,3])
            t4 = ifelse(Stock=="Harrison",Term_Catch[[mm]][[ss]][yy,,ff,1,3], Term_Catch[[mm]][[ss]][yy,,ff,1,4])
            t5 = ifelse(Stock=="Harrison",Term_Catch[[mm]][[ss]][yy,,ff,1,4], Term_Catch[[mm]][[ss]][yy,,ff,1,5]) }
          
          if (ff == 2) {
            t2 <- ifelse(Stock=="Harrison",Term_Catch[[mm]][[ss]][yy,,ff,1,1],0) # ifelse(ff==2 & Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][yy,,ff,1,1]),0),
            t3 <- ifelse(Stock=="Harrison",Term_Catch[[mm]][[ss]][yy,,ff,1,2],0)
            t4 <- ifelse(Stock=="Harrison",Term_Catch[[mm]][[ss]][yy,,ff,1,3],0)
            t5 <- ifelse(Stock=="Harrison",Term_Catch[[mm]][[ss]][yy,,ff,1,4],0) }
          
          
          NewRow <- data.frame(Stock=Stock, Model=Names[mm], Sim = ss, Fish.Type = f.type, Year = yy,
                               PT_Age2 = ifelse(Stock=="Harrison", PT_Catch[[mm]][[ss]][yy,,ff,1,1], PT_Catch[[mm]][[ss]][yy,,ff,1,2]),
                               PT_Age3 = ifelse(Stock=="Harrison", PT_Catch[[mm]][[ss]][yy,,ff,1,2], PT_Catch[[mm]][[ss]][yy,,ff,1,3]),
                               PT_Age4 = ifelse(Stock=="Harrison", PT_Catch[[mm]][[ss]][yy,,ff,1,3], PT_Catch[[mm]][[ss]][yy,,ff,1,4]),
                               PT_Age5 = ifelse(Stock=="Harrison", PT_Catch[[mm]][[ss]][yy,,ff,1,4], PT_Catch[[mm]][[ss]][yy,,ff,1,5]),
                               T_Age2 = t2,
                               T_Age3 = t3,
                               T_Age4 = t4,
                               T_Age5 = t5 )
          HarvRatesDF <- rbind(HarvRatesDF, NewRow)
          
        } # end year loop
        
      } # end fishery loop
      
      HarvRates[[rr]]<- HarvRatesDF
      rr <- rr+1
    } # end sim loop
    rr <- rr+1
    
  } # end model loop
  
  HRatAge <- bind_rows(HarvRates)
  
  write.csv(HRatAge,paste("DataOut/HRatAge_Results_",Run ,".csv", sep=""))
  
}

write_totalHR_year <-function(Names, Run) {
  
  # want to extract table of harvest rates used in the model
  # also save all blobs
  Blobs <- list()
  PT_Catch <- list()
  Term_Catch <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
  }
  
  
  # For each stock, model, sim and Age get mean harvest rates 
  Stock <- Blobs[[1]]$Data$Stocks
  NY <- Blobs[[1]]$Data$NY
  Age <- Blobs[[1]]$Data$Ages
  Years <- Blobs[[1]]$Options$Years
  nSims <- Blobs[[1]]$Options$nSims
  
  
  # Make list
  TotalHR <- list()
  rr<-1                        
  
  for (mm in 1:length(Names)){
    print(Names[[mm]])
    for(ss in 1:nSims){
      
      # Make data frame summarizing differences
      TotalHR_DF <- data.frame(Stock=character(), Model=character(), Sim=numeric(), Year= numeric(), Total_Catch=numeric(),Total_Escape=numeric(), Total_HR =numeric())
        
        for (yy in 1:length(Years)){
          
          totalcatch <- sum(Blobs[[mm]]$Sims$PT_Catch[[ss]][yy,,1,1,]) +sum(Blobs[[mm]]$Sims$PT_Catch[[ss]][yy,,2,1,])
          totalescape <- sum(Blobs[[mm]]$Sims$Escape[[ss]][yy,,])
        
        
        NewRow <- data.frame(Stock=Stock, Model=Names[mm], Sim = ss, Year = yy, 
                             Total_Catch = totalcatch,
                             Total_Escape = totalescape,
                             Total_HR = totalcatch/(totalcatch+totalescape))
        
        TotalHR_DF <- rbind(TotalHR_DF, NewRow)
      
       } #end year loop
      
      TotalHR[[rr]] <- TotalHR_DF
      rr<-rr+1
    } # end sim loop
    rr<- rr+1
  } # end model loop
  
  TotalCatchRate <- bind_rows(TotalHR)
  
  write.csv(TotalCatchRate,paste("DataOut/TotalCatchRate_Results_",Run ,".csv", sep=""))
  
}
##########################################################

##########################################################
# Functions to calculate the avergae HR in each sim
##########################################################

write_HR_atage_sim <-function(Names, Run) {
  
  # want to extract table of harvest rates used in the model
  # also save all blobs
  Blobs <- list()
  PT_Catch <- list()
  Term_Catch <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    PT_Catch[[mm]] <- Blobs[[mm]]$Sims$PT_Rates
    Term_Catch[[mm]] <- Blobs[[mm]]$Sims$Term_Rates
  }
  
  
  # For each stock, model, sim and Age get mean harvest rates 
  Stock <- Blobs[[1]]$Data$Stocks
  NY <- Blobs[[1]]$Data$NY
  Age <- Blobs[[1]]$Data$Ages
  Years <- Blobs[[1]]$Options$Years
  nSims <- Blobs[[1]]$Options$nSims
  start.yr <- length(Blobs[[1]]$Data$Escapement$Year)+1
  
  HarvRates <- list()
  rr<-1
  
  for (mm in 1:length(Names)){
    
    print(Names[[mm]])
    
    for(ss in 1:nSims){
      
      # Make data frame summarizing differences
      HarvRatesDF <- data.frame(Stock=character(), Model=character(), Fish.Type = character(), Sim=numeric(), PT_Age2=numeric(),PT_Age3=numeric(),
                                PT_Age4=numeric(), PT_Age5=numeric(),T_Age2=numeric(),T_Age3=numeric(),T_Age4=numeric(), T_Age5=numeric())
      
      for (ff in 1:2){
        
        
        ifelse(ff == 1, f.type <- "CA", f.type <- "US")
        
        if (ff == 1) {
          t2 = ifelse(Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,1]), mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,2]))
          t3 = ifelse(Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,2]), mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,3]))
          t4 = ifelse(Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,3]), mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,4]))
          t5 = ifelse(Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,4]), mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,5])) }
        
        if (ff == 2) {
          t2 <- ifelse(Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,1]),0) # ifelse(ff==2 & Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][yy,,ff,1,1]),0),
          t3 <- ifelse(Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,2]),0)
          t4 <- ifelse(Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,3]),0)
          t5 <- ifelse(Stock=="Harrison",mean(Term_Catch[[mm]][[ss]][start.yr:NY,,ff,1,4]),0) }
        
        
        NewRow <- data.frame(Stock=Stock, Model=Names[mm], Sim = ss, Fish.Type = f.type,
                             PT_Age2 = ifelse(Stock=="Harrison", mean(PT_Catch[[mm]][[ss]][start.yr:NY,,ff,1,1]), mean(PT_Catch[[mm]][[ss]][start.yr:NY,,ff,1,2])),
                             PT_Age3 = ifelse(Stock=="Harrison", mean(PT_Catch[[mm]][[ss]][start.yr:NY,,ff,1,2]), mean(PT_Catch[[mm]][[ss]][start.yr:NY,,ff,1,3])),
                             PT_Age4 = ifelse(Stock=="Harrison", mean(PT_Catch[[mm]][[ss]][start.yr:NY,,ff,1,3]), mean(PT_Catch[[mm]][[ss]][start.yr:NY,,ff,1,4])),
                             PT_Age5 = ifelse(Stock=="Harrison", mean(PT_Catch[[mm]][[ss]][start.yr:NY,,ff,1,4]), mean(PT_Catch[[mm]][[ss]][start.yr:NY,,ff,1,5])),
                             T_Age2 = t2,
                             T_Age3 = t3,
                             T_Age4 = t4,
                             T_Age5 = t5 )
        HarvRatesDF <- rbind(HarvRatesDF, NewRow)
        
      } # end fishery loop
      
      HarvRates[[rr]]<- HarvRatesDF
      rr <- rr+1
    } # end sim loop
    rr <- rr+1
    
  } # end model loop
  
  HRatAge <- bind_rows(HarvRates)
  
  write.csv(HRatAge,paste("DataOut/HRatAge_SimResults_",Run ,".csv", sep=""))
  
}

write_totalHR_sim <-function(Names, Run) {
  
  # want to extract table of harvest rates used in the model
  # also save all blobs
  Blobs <- list()
  PT_Catch <- list()
  Term_Catch <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
  }
  
  
  # For each stock, model, sim and Age get mean harvest rates 
  Stock <- Blobs[[1]]$Data$Stocks
  NY <- Blobs[[1]]$Data$NY
  Age <- Blobs[[1]]$Data$Ages
  Years <- Blobs[[1]]$Options$Years
  nSims <- Blobs[[1]]$Options$nSims
  start.yr <- length(Blobs[[1]]$Data$Escapement$Year)+1
  
  # Make list
  TotalHR <- list()
  rr<-1                        
  
  for (mm in 1:length(Names)){
    print(Names[[mm]])
    
    for(ss in 1:nSims){
      
      # Make data frame summarizing differences
      TotalHR_DF <- data.frame(Stock=character(), Model=character(), Sim=numeric(), Year= numeric(), Total_Catch=numeric(),Total_Escape=numeric(), Total_HR =numeric())
      
      for (yy in start.yr:length(Years)){
        # Need to add term catch in here too
        totalcatch <- sum(Blobs[[mm]]$Sims$PT_Catch[[ss]][yy,,1,1,]) + sum(Blobs[[mm]]$Sims$PT_Catch[[ss]][yy,,2,1,]) + sum(Blobs[[mm]]$Sims$Term_Catch[[ss]][yy,,])
        totalescape <- sum(Blobs[[mm]]$Sims$Escape[[ss]][yy,,])
        
        
        NewRow <- data.frame(Stock=Stock, Model=Names[mm], Sim = ss, Year = yy, 
                             Total_Catch = totalcatch,
                             Total_Escape = totalescape,
                             Total_HR = totalcatch/(totalcatch+totalescape))
        
        TotalHR_DF <- rbind(TotalHR_DF, NewRow)
        
      } #end year loop
      
      AvgtotHR_DF <- data.frame(Stock=TotalHR_DF$Stock[1], Model=TotalHR_DF$Model[1], Sim=TotalHR_DF$Sim[1], 
                                Avg_Catch=mean(TotalHR_DF$Total_Catch),Avg_Escape=mean(TotalHR_DF$Total_Escape), Avg_HR =mean(TotalHR_DF$Total_HR, na.rm=TRUE))
      
      TotalHR[[rr]] <- AvgtotHR_DF
      rr<-rr+1
    } # end sim loop
    rr<- rr+1
  } # end model loop
  
  TotalCatchRate <- bind_rows(TotalHR)
  
  write.csv(TotalCatchRate,paste("DataOut/TotalCatchRate_SimResults_",Run ,".csv", sep=""))
  
}
##########################################################

##########################################################
# Functions to calculate means
##########################################################

savemeans <- function(Names, Run) {
  
  Blobs <- list()
  Escape <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    Escape[[mm]] <- Blobs[[mm]]$Sims$Escape
  }
  
  # For each stock, model, sim and Age get mean harvest rates 
  Stock <- Blobs[[1]]$Data$Stocks
  NY <- Blobs[[1]]$Data$NY
  Age <- Blobs[[1]]$Data$Ages
  Years <- Blobs[[1]]$Options$Years
  nSims <- Blobs[[1]]$Options$nSims
  
  # Make data frame summarizing differences
  #MeansDF <- data.frame(Stock=character(), Model=character(), Escape=numeric())
  
  MyMean <- list()
  #MySD <- list()
    # first store all model data so can get ylims
    for(mm in 1:length(Names)){
      StockDat <- lapply(Escape[[mm]], "[", , 1 , ) 
      MyDat <- lapply(1:length(StockDat), function(x) apply(StockDat[[x]], 1, sum) )
      MyMean[[mm]] <-  sapply(1:NY, function(x) mean(sapply(MyDat, "[", x) ))
      #MySD[[mm]] <- sapply(1:Data$NY, function(x) sd(sapply(MyDat, "[", x) ))
    } # end mod loop
  
  means.table<- list.cbind(MyMean)
  colnames(means.table)<-Names
  Means <- cbind(Years,means.table)
  
  
  write.csv(Means,paste("DataOut/Means_",Run ,".csv", sep=""))

}


savemedians <- function(Names, Run) {
  
  Blobs <- list()
  Escape <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    Escape[[mm]] <- Blobs[[mm]]$Sims$Escape
  }
  
  # For each stock, model, sim and Age get mean harvest rates 
  Stock <- Blobs[[1]]$Data$Stocks
  NY <- Blobs[[1]]$Data$NY
  Age <- Blobs[[1]]$Data$Ages
  Years <- Blobs[[1]]$Options$Years
  nSims <- Blobs[[1]]$Options$nSims
  
  # Make data frame summarizing differences
  #MeansDF <- data.frame(Stock=character(), Model=character(), Escape=numeric())
  
  MyMedian <- list()
  #MySD <- list()
  
  for(mm in 1:length(Names)){
    StockDat <- lapply(Escape[[mm]], "[", , 1 , ) 
    MyDat <- lapply(1:length(StockDat), function(x) apply(StockDat[[x]], 1, sum) )
    MyMedian[[mm]] <-  sapply(1:NY, function(x) median(sapply(MyDat, "[", x) ))
    #MySD[[mm]] <- sapply(1:Data$NY, function(x) sd(sapply(MyDat, "[", x) ))
  } # end mod loop
  
  medians.table<- list.cbind(MyMedian)
  colnames(medians.table)<-Names
  Medians <- cbind(Years,medians.table)
  
  
  write.csv(Medians,paste("DataOut/Medians_",Run ,".csv", sep=""))
  
}

##########################################################

##########################################################
# Function to calculate results compared to recovery targets
##########################################################

RecoveryResults <-function(Names, Run) { 
  
  # want to extract table of harvest rates used in the model
  # also save all blobs
  Blobs <- list()
  PT_Catch <- list()
  Term_Catch <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
  }
  
  
  # For each stock, model, sim and Age get mean harvest rates 
  Stock <- Blobs[[1]]$Data$Stocks
  NY <- Blobs[[1]]$Data$NY
  Age <- Blobs[[1]]$Data$Ages
  Years <- Blobs[[1]]$Options$Years
  nSims <- Blobs[[1]]$Options$nSims
  inityrs <- Blobs[[1]]$Data$Escapement$Year
  
  #to get start year for abundance to be average over
  a.yrs <- length(Years) - ifelse(Stock == "Harrison" , 3 ,4)
  #to get the start of the % change
  pc.yrs<- length(Years)-length(inityrs) -1
  
  Results <- list()
  rr<-1
  
  for (mm in 1:length(Names)){
    print(Names[mm])
    for(ss in 1:nSims){
      
      #spn <- log(rowSums((Blobs[[mm]]$Sims$Escape[[ss]][c(length(inityrs)+1):length(Years),,])))
      #yr <- 0:pc.yrs
      data.in <- data.frame(spn =rowSums(Blobs[[mm]]$Sims$Escape[[ss]][c(length(inityrs)+1):length(Years),,]), yr = 0:pc.yrs)
      
      p.change <- per.change(data.in)
         
      abund <- mean(rowSums((Blobs[[mm]]$Sims$Escape[[ss]][a.yrs:length(Years),,])))
         
      low_pc <- ifelse(p.change >= -0.3,1,0)
      low_esc <- ifelse(abund >= Blobs[[mm]]$Data$StocksInfo$Low_Esc,1,0)
      up_pc <- ifelse(p.change < -0.3,0,1)
      up_esc <- ifelse(abund >= Blobs[[mm]]$Data$StocksInfo$Up_Esc,1,0)
      Low_RT <- ifelse(low_esc ==1 & low_pc ==1, 1, 0)
      Upp_RT <- ifelse(up_esc ==1 & up_pc ==1, 1,0)
          
      ResultsDF <- data.frame(Stock=Stock, Model=Names[mm], FishScenario = Blobs[[mm]]$Options$HRS , ProdScenario=Blobs[[mm]]$Options$Prod_Scenario , Sim = ss, Abundance = abund, perc.Change = p.change, Low_esc_met = low_esc, Low_pc_met = low_pc,
                               Low_RT_met = Low_RT, Upp_esc_met = up_esc, Upp_pc_met = up_pc, Upp_RT_met= Upp_RT)
          
      Results[[rr]] <- ResultsDF
      rr<- rr+1
    } # end sim loop
    rr<- rr+1
  } # end model loop
  
  Recovery_Results <- bind_rows(Results)
  
  write.csv(Recovery_Results,paste("DataOut/Recovery_Results_",Run ,".csv", sep=""))
  
}


##########################################################

##########################################################
# Functions to calculate the average HR in each model
##########################################################

write_HR_atage_avg <-function(Stock="DU2",Run) {
  
  HR_atage <- read.csv("DataOut/HRatAge_SimResults_DU2.csv")
     
  # For each stock, model, sim and Age get mean harvest rates 
  Models <- unique(HR_atage$Model)
  fish <- unique(HR_atage$Fish.Type)
  
  # Make data frame summarizing differences
  HarvRatesDF <- data.frame(Stock=character(), Model=character(), Fish.Type = character(), PT_Age2=numeric(),PT_Age3=numeric(),
                            PT_Age4=numeric(), PT_Age5=numeric(),T_Age2=numeric(),T_Age3=numeric(),T_Age4=numeric(), T_Age5=numeric())
  
  for (mm in 1:length(Models)){
    
      for (ff in 1:length(fish)){
        
        NewRow <- data.frame(Stock=Stock, Model=Models[mm], Fish.Type = fish[ff],
                             PT_Age2 = mean(HR_atage$PT_Age2[HR_atage$Model == Models[mm] & HR_atage$Fish.Type == fish[ff]]),
                             PT_Age3 = mean(HR_atage$PT_Age3[HR_atage$Model == Models[mm] & HR_atage$Fish.Type == fish[ff]]),
                             PT_Age4 = mean(HR_atage$PT_Age4[HR_atage$Model == Models[mm] & HR_atage$Fish.Type == fish[ff]]),
                             PT_Age5 = mean(HR_atage$PT_Age5[HR_atage$Model == Models[mm] & HR_atage$Fish.Type == fish[ff]]),
                             T_Age2 = mean(HR_atage$T_Age2[HR_atage$Model == Models[mm] & HR_atage$Fish.Type == fish[ff]]),
                             T_Age3 = mean(HR_atage$T_Age3[HR_atage$Model == Models[mm] & HR_atage$Fish.Type == fish[ff]]),
                             T_Age4 = mean(HR_atage$T_Age4[HR_atage$Model == Models[mm] & HR_atage$Fish.Type == fish[ff]]),
                             T_Age5 = mean(HR_atage$T_Age5[HR_atage$Model == Models[mm] & HR_atage$Fish.Type == fish[ff]]) )
        
        HarvRatesDF <- rbind(HarvRatesDF, NewRow)
        
      } # end fishery loop
    
  } # end model loop
  
  write.csv(HarvRatesDF,paste("DataOut/HRatAge_Avg_",Run ,".csv", sep=""))
  
}

write_totalHR_avg <-function(Stock,Run) {
  
  HR_total <- read.csv("DataOut/TotalCatchRate_SimResults_DU2.csv")
  
  # For each stock, model, sim and Age get mean harvest rates 
  Models <- unique(HR_total$Model)

  # Make data frame summarizing differences
  TotalHR_DF <- data.frame(Stock=character(), Model=character(), Avg_Total_HR =numeric())
  
  for (mm in 1:length(Models)){
    
        NewRow <- data.frame(Stock=Stock, Model=Models[mm], Avg_Total_HR = mean(HR_total$Avg_HR[HR_total$Model==Models[mm]]))
        
        TotalHR_DF <- rbind(TotalHR_DF, NewRow)
        
  } # end model loop
  
  write.csv(TotalHR_DF,paste("DataOut/AvgTotalCatch_", Run ,".csv", sep=""))
  
}

##########################################################
