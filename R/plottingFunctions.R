################################
#   CN MSE Plotting Functions
#  Brooke Davis & Kendra Holt
################################

# ==================================================================================================
# Plot functions included in this file:
# ================================================================================================
# CohortPlots ()              - plots cohort abundance trajectories at age and by stock fr one or more scenarios  
# EscapePlots ()              - plots escapement trajectories by stock for one or more scenarios
# Escape_4_Panel ()           - plots escapement trajectories across scenarios
# PlotEscapeChange()          - dot and line plots showing the ratio of change in escapement (future:recent)
# PlotCatch()                 - dot and line plots showing average Pre-terminal, Terminal, and Total catch (averaged over all years) 
# checkNSims ()               - boxplots to check number of simulation replicates needed to get stable results
# Catch_4_Panel()             - plots catch over time, across scenarios
# AI_4_Panel()                - plots AABM abundance index (AI) over time, across scenarios
# ==================================================================================================



###########################################################
#  CohortPlots
#######################################################
# Purpose: Plots cohort abundance trajectory at age and by stock for one or more scenarios 
# Arguments: Names = names of scenarios to be plotted; each name must have an rds file with that name saved in DataOut folder
#           PlotName = name that will be used when saving plots as a pdf file in Figures folder             
################################################################

CohortPlots <- function(Names, PlotName, Legend_Names=NA){
 
  # want to extract table of median change in escapment over 10 years (2016 to 2025) and 20 years (2016 to 2035)
  # also save all blobs
  Blobs <- list()
  Escape <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    Escape[[mm]] <- Blobs[[mm]]$Sims$Escape
  }
  
  #Colors -- need to plot less than 7 scenarios or not enough colors
  cols <- c("#0000ff", "#b22222", "#7E3FCC", "#006400", "#FFC725", "#808080", "#EF29E4")
  TGrey <- "#80808050"
  
  # if legend names given, use these, otherwise use scenario names
  if(is.na(Legend_Names[1])){
    Legend_Names <- Names
  }
  
  CohortDat <- AbundDat <- list()
  for(mm in 1:length(Blobs)){
    CohortDat[[mm]] <- Blobs[[mm]]$Sims$Cohort
    AbundDat[[mm]] <-  Blobs[[mm]]$Sims$N
  }
  
  # Get require info from first blob, assume they are the same
  Opts <- Blobs[[1]]$Options
  Data <- Blobs[[1]]$Data
    
    # For each Stock, age plot Cohort over time of each
    pdf(paste("Figures/",PlotName,".pdf", sep=""))
    par(mfrow=c(2,3), oma=c(3,3,2,2), mar=c(2,2,2,2))
    for(ss in 1:length(Data$Stocks)){
      for(aa in 1:Data$MaxAge){
        MyMean <- list()
        MySD <- list()
        for(mm in 1:length(Blobs)){
          if(aa == 1){
            MyDat <- lapply(CohortDat[[mm]], "[", , ss, aa)
          } else {
            # Data needs to come from N
            SA_Dat <- lapply(AbundDat[[mm]], "[", ,1, , ss, aa-1)
            # now need to sum over regions
            if(Data$NR > 1 ){
              MyDat <- lapply(1:length(SA_Dat), function(x) apply(SA_Dat[[x]], 1, sum) )
            } else {
              MyDat <- SA_Dat
            }
          } # end else age >=2
          # Get means and SDs
          MyMean[[mm]] <-  sapply(1:Data$NY, function(x) mean(sapply(MyDat, "[", x) ))
          MySD[[mm]] <- sapply(1:Data$NY, function(x) sd(sapply(MyDat, "[", x) ))
        } # end model loop -- now that have all data, set ylims
        
        if(Opts$nSims >1 ){
          ylims <- c( min(c(unlist(MyMean) - unlist(MySD))), max(c( unlist(MyMean) + unlist(MySD) )))
        } else {
          ylims <- c( min(c( unlist(MyMean) )), max(c( unlist(MyMean) )))
        } # end nsims if
        # remove initialization years if applicable
        if(Opts$Initialization == "Escapement"){
          Years_To_Plot <- Opts$Years[-(1:Data$MaxAge)]
          MyMean <- lapply(MyMean, "[", which(Opts$Years %in% Years_To_Plot))
          MySD <- lapply(MySD, "[", which(Opts$Years %in% Years_To_Plot))
        } else {
          Years_To_Plot <- Opts$Years
        }
        
        # Plot first model, just make black
        plot(Years_To_Plot, MyMean[[1]], type="l", ylim=ylims)
        polygon(y=c(MyMean[[1]] - MySD[[1]], rev(MyMean[[1]] + MySD[[1]])), 
                x=c(Years_To_Plot, rev(Years_To_Plot)), col=TGrey, border="black")
        
         if(length(Names) > 1){
          for(mm in 2:length(Names)){
            lines(Years_To_Plot, MyMean[[mm]], col= cols[mm-1])
            # Add tranparent error bars
            polygon(y=c(MyMean[[mm]] - MySD[[mm]], rev(MyMean[[mm]] + MySD[[mm]])), x=c(Years_To_Plot, rev(Years_To_Plot)), 
                    col=paste(cols[mm-1], "50", sep=""), border=cols[mm-1])
          } # end mod loop
        } # if more than one mod if
        mtext(side=3, text=paste( "Age", aa, sep=" "))
      } # end age loop
      
      # if max is age 6 put legend on last plot
      if(Data$MaxAge==6){
        legend("topright", col=c("black", cols[1:(length(Blobs)-1)]), legend=Legend_Names, bty="n")
      } else {
        # else if 5 is max age put on last black plotting region
        plot.new()
        legend("center", fill=c("black", cols[1:(length(Blobs)-1)]), legend=Legend_Names)
      }
      mtext(side=3, outer=T, text=Data$Stocks[ss])
    } # End stock loop
    dev.off()
} # end cohort plot function


###########################################################
#  EscapePlots
#######################################################
# Purpose: Plots escapement trajectory by stock for one or more scenarios 
#             (note: if using > 1 scenario, escapement trajectories for each scnenario will be shown in the same plot) 
# Arguments: Names = names of scenarios to be plotted; each name must have an rds file with that name saved in DataOut folder
#            PlotName = name that will be used when saving plots as a pdf file in Figures folder 
#            LeadIn = Indicates whether lead-in (pre-initiaization) escapements should be plotted (True or False)
#            Smax = Whether or not to add Smax to time-series
################################################################

EscapePlots <- function(Names, PlotName, LeadIn = T, Legend_Names=NA){
  
  # want to extract table of median change in escapment over 10 years (2016 to 2025) and 20 years (2016 to 2035)
  # also save all blobs
  Blobs <- list()
  Escape <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    Escape[[mm]] <- Blobs[[mm]]$Sims$Escape
  }

  
  #Colors -- need to plot less than 7 or not enough colors
  cols <- c("#0000ff", "#b22222", "#7E3FCC", "#006400", "#FFC725", "#808080", "#EF29E4")
  TGrey <- "#80808050"
  
  # if legend names given, use these, otherwise use scenario names
  if(is.na(Legend_Names[1])){
    Legend_Names <- Names
  }
  
  Escape <- AbundDat <- list()
  for(mm in 1:length(Blobs)){
    Escape[[mm]] <- Blobs[[mm]]$Sims$Escape
  }
  
  # also read in "leadin" escapement data
  if(LeadIn == T){
    Esc_LeadIn <- read.csv("DataIn/Escape_LeadIn.csv")
  }
  
  # Get require info from first blob, assume they are the same
  Opts <- Blobs[[1]]$Options
  Data <- Blobs[[1]]$Data
  
  pdf(paste("Figures/", PlotName ,".pdf", sep=""))
  par(mfrow=c(2,2), oma=c(2,2,2,1), mar=c(2,2,2,1))
  
  for(ss in 1:length(Data$Stocks)){
    
    MyMean <- list()
    MySD <- list()
    # first store all model data so can get ylims
    for(mm in 1:length(Blobs)){
      StockDat <- lapply(Escape[[mm]], "[", ,ss , ) 
      MyDat <- lapply(1:length(StockDat), function(x) apply(StockDat[[x]], 1, sum) )
      MyMean[[mm]] <-  sapply(1:Data$NY, function(x) mean(sapply(MyDat, "[", x) ))
      MySD[[mm]] <- sapply(1:Data$NY, function(x) sd(sapply(MyDat, "[", x) ))
    } # end mod loop
    # remove initialization years if applicable
    if(Opts$Initialization == "Escapement"){
      Years_To_Plot <- Opts$Years[-(1:Data$MaxAge)]
      MyMean <- lapply(MyMean, "[", which(Opts$Years %in% Years_To_Plot))
      MySD <- lapply(MySD, "[", which(Opts$Years %in% Years_To_Plot))
    } else {
      Years_To_Plot <- Opts$Years
    }
    
    # Get ylims
    if(Opts$nSims >1 ){
      if(LeadIn ==F){
        ylims <- c( 0, max(c( unlist(MyMean) + unlist(MySD) )))
      } else {
        Edat <- Esc_LeadIn[which(as.character(Esc_LeadIn$StockID)==as.character(Data$Stocks[ss])),]
        ylims <- c( 0 , max(c( unlist(MyMean) + unlist(MySD), na.omit(Edat$Escape) )))
      }
    } else {
      ylims <-  c( min(c( unlist(MyMean) )), max(c( unlist(MyMean) )))
    }  
    # get escapament init data
    Edat_Init <- Data$Escapement[which(Data$Escapement$StockID==Data$Stocks[ss]),]
    # Plot
    if(LeadIn==T){
      plot(Edat$Year, Edat$Escape, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Edat$Year), 
                                                                                max(Opts$Years)), lwd=2)
      # now add initialization data
      lines(Edat_Init$Year, Edat_Init$Abundance, col="black", lwd=2)
    } else {
      plot(Edat_Init$Year, Edat_Init$Abundance, col="black", type="l", ylim=ylims, lwd=2, 
           xlim=c(min(Opts$Years), max(Opts$Years)))
    }
    # Now add model projections
    lines(Years_To_Plot, MyMean[[1]], col="darkgrey", type="l", lwd=2.5)
    polygon(y=c(MyMean[[1]] - MySD[[1]], rev(MyMean[[1]] + MySD[[1]])), 
            x=c(Years_To_Plot, rev(Years_To_Plot)), col=TGrey, border="darkgrey")
    
    if(length(Blobs) >= 2){
      for(mm in 2:length(Blobs)){
        lines(Years_To_Plot, MyMean[[mm]], col=cols[mm-1], lwd=2.5)
        # end mods loop
        # Add tranparent error bars
        polygon(y=c(MyMean[[mm]] - MySD[[mm]], rev(MyMean[[mm]] + MySD[[mm]])), 
                x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(cols[mm-1], 50, sep=""), border=cols[mm-1])
      } # end mod loop
    }
    # label stock
    mtext(side=3, text=Data$Stocks[ss])
    # if last on page do overall labels
    if(ss %% 4 == 0){
      mtext(side=1, text="Year", outer=T, line=0.5)
      mtext(side=2, text = "Escapement", outer=T, line=0.5)
      #mtext(side=3, text="Escapement and Targets", outer=T, line=0.5)
      if(length(Blobs) >= 2){
        legend("topright", fill=c("black", cols[1:(length(Blobs)-1)]), legend=Legend_Names)
      } # and more than one model if
    } # end final plot of page if
    
  } # end Stock loop
  dev.off()
} # End Escape Plot function



#*###########################################################
#  Escape_4_Panel
#######################################################
# Purpose: Plots escapement trajectories for the 2 OMs and 6 MPs shown at April 2018 meetings
# Arguments: Names = names of scenarios to be plotted; each name must have an rds file with that name saved in DataOut folder
#            Esc_LeadIn = Indicates whether lead-in (pre-initiaization) escapements should be plotted (True or False)
#            PlotName = name that will be used when saving plots as a pdf file in Figures folder 
################################################################
# Use Esc_medians and Esc_SDs from compiled outputs (indexed by mm,ss) 

Escape_4_Panel <- function(Names, Esc_LeadIn, PlotName="Escape_4Panel", Smax = T, Plot_Spawners=F,
                           Mods2Plot_1 = c(1,2,3,4), Mods2Plot_2 = c(1,5,6), Mods2Plot_3 = c(7,8,9,10), 
                           Mods2Plot_4 = c(7,11,12), Scenarios = "Fishing_Habitat", BiasNum) {
  
  
  # want to extract table of median change in escapment over 10 years (2016 to 2025) and 20 years (2016 to 2035)
  # also save all blobs
  Blobs <- list()
  Escape <- list()
  Spawners <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    Escape[[mm]] <- Blobs[[mm]]$Sims$Escape
    if(Plot_Spawners==T){
      Spawners[[mm]] <- Blobs[[mm]]$Sims$Spawners
    }
  }
  
  cols <- c("#0000ff", "#b22222", "#006400", "#FF9937","#7F3DCC",  "#808080")
  #  Blue, Green, Red, Orange, Purple, grey
  TGrey <- "#80808050"

  
  pdf(paste("Figures/", PlotName ,".pdf", sep=""))
  par(mfrow=c(2,2), oma=c(2,3,2,1), mar=c(2,2,2,1))

  Stocks <- Blobs[[1]]$Data$StocksInfo$StockID
  Years <- Blobs[[1]]$Options$Years
  NY <- Blobs[[1]]$Data$NY
  MaxAge <- Blobs[[1]]$Data$MaxAge
  Escapement <- Blobs[[1]]$Data$Escapement
  
  StocksInfo_List <- list()
  for(mm in 1:length(Blobs)){
    StocksInfo_List[[mm]] <- Blobs[[1]]$Data$StocksInfo
  }
  
  
  for(ss in 1:length(Stocks)){
    MyMean <- list()
    MySD <- list()
    # first store all model data so can get ylims
    for(mm in 1:length(Blobs)){
      if(Plot_Spawners==T){
        MyDat <- lapply(Spawners[[mm]], "[", ,ss  ) 
      } else {
        StockDat <- lapply(Escape[[mm]], "[", ,ss , ) 
        MyDat <- lapply(1:length(StockDat), function(x) apply(StockDat[[x]], 1, sum) )
      }
      MyMean[[mm]] <-  sapply(1:NY, function(x) mean(sapply(MyDat, "[", x) ))
      MySD[[mm]] <- sapply(1:NY, function(x) sd(sapply(MyDat, "[", x) ))
    } # end mod loop
    # remove initialization years
      Years_To_Plot <- Years[-(1:MaxAge)]
      MyMean <- lapply(MyMean, "[", which(Years %in% Years_To_Plot))
      MySD <- lapply(MySD, "[", which(Years %in% Years_To_Plot))
    # Get ylims
        Edat <- Esc_LeadIn[which(Esc_LeadIn$StockID==Stocks[ss] & Esc_LeadIn$Year < min(Years)), ]
        # need wider bounds for some stocks to make room for legend
      if(Stocks[ss] == "FRsp5.2") {
        ylims <- c( 0 , max(c( unlist(MyMean) + unlist(MySD), Edat$Escape )*1.25, na.rm=T))
      } else if(Stocks[ss]=="FRsu5.2") {
        ylims <- c( 0 , max(c( unlist(MyMean) + unlist(MySD), Edat$Escape )*1.4, na.rm=T))
      } else {
        ylims <- c( 0 , max(c( unlist(MyMean) + unlist(MySD), Edat$Escape ), na.rm=T))
      }
    # get escapament init data
    Edat_Init <- Escapement[which(Escapement$StockID==Stocks[ss]),]
    # Plot
    #**************************************************
    # Panel 1 -- Basic Ricker with 4 fishing scenarios OR smolt release scenarios
    Cols <- c("#000000", cols[1:3])
    # start with leadin and inits
    plot(Edat$Year, Edat$Escape, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Edat$Year), max(Years)), lwd=2, ann=F)
    # now add initialization data
    lines(Edat_Init$Year, Edat_Init$Abundance, col="black", lwd=2)
    # Now add model fits
      for(mm in 1:length(Mods2Plot_1)){
        lines(Years_To_Plot, MyMean[[Mods2Plot_1[mm]]], col=Cols[mm], lwd=2)
        # end mods loop
        # Add tranparent error bars
        polygon(y=c(MyMean[[Mods2Plot_1[mm]]] - MySD[[Mods2Plot_1[mm]]], rev(MyMean[[Mods2Plot_1[mm]]] + MySD[[Mods2Plot_1[mm]]])), 
                x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
      } # end mod loop
    # add Lines for capacity, Smax
    if(Smax == T & StocksInfo_List[[mm]]$Type[which(StocksInfo_List[[mm]]$StockID==Stocks[ss])] == "Natural"){
      # extract alpha and beta
      # only works if natural stock
      alpha <- StocksInfo_List[[Mods2Plot_1[mm]]]$Ricker_A[which(StocksInfo_List[[Mods2Plot_1[mm]]]$StockID==Stocks[ss])]
      beta <- StocksInfo_List[[Mods2Plot_1[mm]]]$Ricker_B[which(StocksInfo_List[[Mods2Plot_1[mm]]]$StockID==Stocks[ss])]
      abline(h=1/beta, lty=2, col="orange", lwd=2)
      abline(h=1/beta*log(alpha), lty=3, col="orange", lwd=2)
    }
    
    if(Smax==T){
     if(Scenarios=="Fishing_Habitat"){
      legend("topleft", col=c(Cols, "white", "orange", "orange"), 
             legend=c("Current", "50% Curr.", "50% Curr. PT", "No Fishing", "", "Smax", "Capacity"), title="Fishing Scenario",
             bty="n", cex=0.8, lwd=2, lty=c(1,1,1,1,1,2,3))
     } else if(Scenarios=="Bias_Smolts"){
       legend("topleft", col=c(Cols, "white", "orange", "orange"), 
              legend=c("Current", "150% Curr.", "150% Curr. 2026", "", "Smax", "Capacity"), title="Smolt Release Scenario",
              bty="n", cex=0.8, lwd=2, lty=c(1,1,1,1,2,3))
     }
    } else {
      if(Scenarios=="Fishing_Habitat"){
        legend("topleft", col=c(Cols), 
               legend=c("Current", "50% Curr.", "50% Curr. PT", "No Fishing"), title="Fishing Scenario",
               bty="n", cex=0.8, lwd=2, lty=c(1,1,1,1))
      } else if(Scenarios=="Bias_Smolts"){
        legend("topleft", col=c(Cols), 
               legend=c("Current", "150% Curr.", "150% Curr. 2026"), title="Smolt Release Scenario",
               bty="n", cex=0.8, lwd=2, lty=c(1,1,1))
      }
    }
    mtext(side=2, line=3, text="Historical Productivity")
    
    #**************************************************
    # Panel 2 -- Basic Ricker with 2 Productivity/CC scenarios OR AI bias scenarios
    Cols <- c("#000000", cols[4:5])
    # start with leadin and inits
    plot(Edat$Year, Edat$Escape, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Edat$Year), max(Years)), lwd=2, ann=F)
    # now add initialization data
    lines(Edat_Init$Year, Edat_Init$Abundance, col="black", lwd=2)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_2)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_2[mm]]], col=Cols[mm], lwd=2)
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_2[mm]]] - MySD[[Mods2Plot_2[mm]]], rev(MyMean[[Mods2Plot_2[mm]]] + MySD[[Mods2Plot_2[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
      # This is a little more complicated for this because Smax/capacity will change
      # will use final year
      if(Smax == T & StocksInfo_List[[mm]]$Type[which(StocksInfo_List[[mm]]$StockID==Stocks[ss])] == "Natural"){
        # extract alpha and beta
        alpha <- StocksInfo_List[[Mods2Plot_2[mm]]]$Ricker_A[which(StocksInfo_List[[Mods2Plot_2[mm]]]$StockID==Stocks[ss])] * Blobs[[Mods2Plot_2[mm]]]$Data$Prod_Mults$Mult[length(Years)]
        beta <- StocksInfo_List[[Mods2Plot_2[mm]]]$Ricker_B[which(StocksInfo_List[[Mods2Plot_2[mm]]]$StockID==Stocks[ss])] / Blobs[[Mods2Plot_2[mm]]]$Data$CC_Mults$Mult[length(Years)]
        abline(h=1/beta, lty=2, col=Cols[mm], lwd=2)
        abline(h=1/beta*log(alpha), lty=3, col=Cols[mm], lwd=2)
      }
    } # end mod loop
   
    if(Smax == T){
      if(Scenarios=="Fishing_Habitat"){
        legend("topleft", col=c(Cols, "white", "grey", "grey"), legend=c("No Change", "Prod. 25%", "CC 25%", "", "Smax", "Capacity"  ), 
           title="Habitat Scenario", bty="n", cex=0.8, lwd=2, lty=c(1,1,1,1,2,3))
      } else if(Scenarios=="Bias_Smolts"){
        legend("topleft", col=c(Cols, "white", "grey", "grey"), legend=c("No Bias", "+15% AI Bias", "-15% Ai Bias", "", "Smax", "Capacity"  ), 
               title="AI Bias Scenarios", bty="n", cex=0.8, lwd=2, lty=c(1,1,1,1,2,3))
      }
    } else { #if Smax=F
      if(Scenarios=="Fishing_Habitat"){
       legend("topleft", col=Cols, legend=c("No Change", "Prod. 25%", "CC 25%"  ), 
             title="Habitat Scenario", bty="n", cex=0.8, lwd=2, lty=c(1,1,1))
      } else if(Scenarios=="Bias_Smolts"){
        legend("topleft", col=Cols, legend=c("No Bias", paste("+", BiasNum ,"% AI Bias", sep=""), paste("-", BiasNum, "% Ai Bias", sep="")), 
               title="AI Bias Scenario", bty="n", cex=0.8, lwd=2, lty=c(1,1,1))
      }
    }
    
    
    #**************************************************
    # Panel 3 --  R.B with 4 fishing scenarios
    Cols <- c("#000000", cols[1:3])
    # start with leadin and inits
    plot(Edat$Year, Edat$Escape, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Edat$Year), max(Years)), lwd=2, ann=F)
    # now add initialization data
    lines(Edat_Init$Year, Edat_Init$Abundance, col="black", lwd=2)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_3)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_3[mm]]], col=Cols[mm], lwd=2)
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_3[mm]]] - MySD[[Mods2Plot_3[mm]]], rev(MyMean[[Mods2Plot_3[mm]]] + MySD[[Mods2Plot_3[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
    } # end mod loop
    if(Smax == T  & Smax == T & StocksInfo_List[[mm]]$Type[which(StocksInfo_List[[mm]]$StockID==Stocks[ss])] == "Natural"){
      # extract alpha and beta
      alpha <- StocksInfo_List[[Mods2Plot_3[mm]]]$Ricker_A[which(StocksInfo_List[[Mods2Plot_3[mm]]]$StockID==Stocks[ss])]
      beta <- StocksInfo_List[[Mods2Plot_3[mm]]]$Ricker_B[which(StocksInfo_List[[Mods2Plot_3[mm]]]$StockID==Stocks[ss])]
      abline(h=1/beta, lty=2, col="orange", lwd=2)
      abline(h=1/beta*log(alpha), lty=3, col="orange", lwd=2)
    }
    mtext(side=2, line=3, text="Recent Productivity")
    
    #***************************************************
    # Panel 4 -- R.B with 2 productivity/CC scenarios
    Cols <- c("#000000", cols[4:5])
    # start with leadin and inits
    plot(Edat$Year, Edat$Escape, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Edat$Year), max(Years)), lwd=2, ann=F)
    # now add initialization data
    lines(Edat_Init$Year, Edat_Init$Abundance, col="black", lwd=2)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_4)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_4[mm]]], col=Cols[mm], lwd=2)
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_4[mm]]] - MySD[[Mods2Plot_4[mm]]], rev(MyMean[[Mods2Plot_4[mm]]] + MySD[[Mods2Plot_4[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
      # add Smax/capacity lines
      if(Smax == T  & Smax == T & StocksInfo_List[[mm]]$Type[which(StocksInfo_List[[mm]]$StockID==Stocks[ss])] == "Natural"){
        # extract alpha and beta
        alpha <- StocksInfo_List[[Mods2Plot_4[mm]]]$Ricker_A[which(StocksInfo_List[[Mods2Plot_4[mm]]]$StockID==Stocks[ss])] * Blobs[[Mods2Plot_4[mm]]]$Data$Prod_Mults$Mult[length(Years)]
        beta <- StocksInfo_List[[Mods2Plot_4[mm]]]$Ricker_B[which(StocksInfo_List[[Mods2Plot_4[mm]]]$StockID==Stocks[ss])] / Blobs[[Mods2Plot_4[mm]]]$Data$CC_Mults$Mult[length(Years)]
        abline(h=1/beta, lty=2, col=Cols[mm], lwd=2)
        abline(h=1/beta*log(alpha), lty=3, col=Cols[mm], lwd=2)
      }
    } # end mod loop
   
    # label stock
    mtext(side=3, text=Stocks[ss], outer=T)
    # add overall labels
      mtext(side=1, text="Year", outer=T, line=0.5)
      mtext(side=2, text = "Escapement", outer=T)
  } # end Stock loop
  dev.off()
  
} # End Escape_4_Panel function


#*###########################################################
#  plotEscapeChange
#######################################################
# Purpose: Barplots showing the ratio of change in escapement (future:recent) for the 2 OMs and 6 MPs shown at April 2018 meetings
# Arguments: Names = names of scenarios to be plotted; each name must have an rds file with that name saved in DataOut folder
#            PlotName = name that will be used when saving plots as a pdf file in Figures folder 
# Note: Lots of hardwiring in this function; will need to be updated to be more generic
################################################################
plotEscapeChange<-function(Names, PlotName) {
 
  # want to extract table of median change in escapment over 10 years (2016 to 2025) and 20 years (2016 to 2035)
  # also save all blobs
  Blobs <- list()
  Escape <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    Escape[[mm]] <- Blobs[[mm]]$Sims$Escape
  }
  
  Esc_Init <- Blob$Data$Escapement

  pdf(paste("Figures/", PlotName ,".pdf", sep=""))
  par(mfrow=c(2,1), oma=c(2,3,2,1), mar=c(2,2,2,1))
  
  # For each stock, model, get median escapements 
  # also get sims of 2022:2025 and 2032:2035
  Stocks <- Blob$Data$Stocks
  NY <- Blob$Data$NY
  Years <- Blob$Options$Years
  Esc_Medians <- list()
  Esc_SDs <- list()
  Esc_22_25 <- list()
  Esc_32_35 <- list()
  EscChange_2025 <- list()
  EscChange_2035 <- list()
  
  for(mm in 1:length(Names)){
    Esc_Medians[[mm]] <- list()
    Esc_SDs[[mm]] <- list()
    Esc_22_25[[mm]] <- list()
    Esc_32_35[[mm]] <- list()
    
    EscChange_2025[[mm]] <- list()
    EscChange_2035[[mm]] <- list()

    for(ss in 1:length(Stocks)){
      
      # Get initial avg escapement from data inputs
      Esc_Start <- mean(Esc_Init$Abundance[which(Esc_Init$StockID==Stocks[ss] & Esc_Init$Year %in% 2012:2015)])
      
      StockDat <- lapply(Escape[[mm]], "[", ,ss , ) 
      
      # sum over ages
      StockDat_Summed <- lapply(1:length(StockDat), function(x) apply(StockDat[[x]], 1, sum) )
      Esc_Medians[[mm]][[ss]] <-  sapply(1:NY, function(x) mean(sapply(StockDat_Summed, "[", x) ))
      Esc_SDs[[mm]][[ss]] <- sapply(1:NY, function(x) sd(sapply(StockDat_Summed, "[", x) ))
      # extract 100 sims of averages
      Esc_22_25[[mm]][[ss]] <-  sapply(1:nSims, function(x) mean((StockDat_Summed[[x]][which(Years %in% 2022:2025)])) )
      Esc_32_35[[mm]][[ss]]<-  sapply(1:nSims, function(x) mean((StockDat_Summed[[x]][which(Years %in% 2032:2035)])) )
      
      EscChange_2025[[mm]][[ss]]<-Esc_22_25[[mm]][[ss]]/Esc_Start
      EscChange_2035[[mm]][[ss]]<-Esc_32_35[[mm]][[ss]]/Esc_Start
      
      
      # TO FIX: need to fix this so not hardwired!
      ifelse(mm<=6, om<-"Average", om<-"Recent")
      if (mm == 1 | mm == 7) mp<-"CurrentF"
      if (mm == 2 | mm == 8) mp<-"Fish50"    
      if (mm == 3 | mm == 9) mp<-"PT_Fish50"
      if (mm == 4 | mm == 10) mp<-"NoFishing"        
      if (mm == 5 | mm == 11) mp<-"IncProd"
      if (mm == 6 | mm == 12) mp<-"IncCap"
      
      EscChange_10y_med<-median(EscChange_2025[[mm]][[ss]])
      EscChange_20y_med<-median(EscChange_2035[[mm]][[ss]])
      
      EscChange_10y_10quant<-quantile(EscChange_2025[[mm]][[ss]],0.1)[[1]]
      EscChange_20y_10quant<-quantile(EscChange_2035[[mm]][[ss]],0.1)[[1]]
      
      EscChange_10y_90quant<-quantile(EscChange_2025[[mm]][[ss]],0.9)[[1]]
      EscChange_20y_90quant<-quantile(EscChange_2035[[mm]][[ss]],0.9)[[1]]
      
      # Initiate data frame
      if (ss == 1 & mm == 1) {
        
        
        df<-data.frame(Stocks[ss],om, mp, EscChange_10y_med, EscChange_20y_med, EscChange_10y_10quant, EscChange_20y_10quant, 
                       EscChange_10y_90quant, EscChange_20y_90quant)
        
      } else {
        tmp<-data.frame(Stocks[ss],om, mp, EscChange_10y_med, EscChange_20y_med, EscChange_10y_10quant, EscChange_20y_10quant, 
                        EscChange_10y_90quant, EscChange_20y_90quant)
        df<-rbind(df,tmp)
      }
      
    } # end stock loop
  } # end mod loop
  
  names(df)<-c("Stock" ,"OM", "MP", "med_10y", "med_20y", "lower_10y", "lower_20y", "upper_10y", "upper_20y")
  
  # Remove hatchery stocks
  Hatchery_Stocks <- Blob$Data$StocksInfo$StockID[which(Blob$Data$StocksInfo$Type=="Hatchery")]
  df <- df[which(df$Stock %notin% Hatchery_Stocks),]
  
  # Plot Average Productivity
  df.ave<-df[which(df$OM=="Average"),]
  # Plot Recent Productivity
  df.rec<-df[which(df$OM=="Recent"),]
  
  mpList <- c("CurrentF", "Fish50", "PT_Fish50", "NoFishing", "IncProd", "IncCap")
  stockList <- unique(df$Stock)
  
  # Create plot for 10-year medians
  
  colList<-c("black","red", "blue","steelblue2","purple","forestgreen")
  
  # Plot OM for Average Productivity
  xlim<-c(1,(length(stockList)* (length(mpList)+1)) )
  ylim<-range(df.ave$lower_10y*0.95,df.ave$upper_10y*1.1)
  
  plot(xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )
  
  for (ss in 1:length(stockList)) {
    
    ind <- (1:length(mpList))+(length(mpList)*(ss-1))+(ss-1)
  
    points(ind,df.ave$med_10y[which(df.ave$Stock==stockList[ss])],pch=19,col=colList)
    arrows(ind,df.ave$lower_10y[which(df.ave$Stock==stockList[ss])],ind,df.ave$upper_10y[which(df.ave$Stock==stockList[ss])],
      length=0, col=colList) 
    axis(side=2)
    
    abline(v=ind[length(mpList)]+1,lty=3, lwd=2)

    box()
    
  } # end stock loop
  
  mtext(side=3, "Average Productivity")
  
  # need legend
  legend("topright", bg = "white", col=colList, legend=mpList, pch=19, cex=0.8)
  
  
  # Plot OM for Recent Productivity
  ylim<-range(df.rec$lower_10y*0.95,df.rec$upper_10y*1.1)
  
  plot(xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )
  
  for (ss in 1:length(stockList)) {
    
    ind<-(1:length(mpList))+(length(mpList)*(ss-1))+(ss-1)
    
    points(ind,df.rec$med_10y[which(df.rec$Stock==stockList[ss])],pch=19,col=colList)
    arrows(ind,df.rec$lower_10y[which(df.rec$Stock==stockList[ss])],ind,df.rec$upper_10y[which(df.rec$Stock==stockList[ss])],
           length=0, col=colList) 
    axis(side=2)
    
    abline(v=ind[6]+1,lty=3, lwd=2)
    
    box()
    
  } # end stock loop
  
  mtext(side=1, at=3 + (length(mpList)+1)*c(0:(length(stockList)-1)), text=stockList, las=2)
  mtext(side=3, "Recent Productivity")
  
  mtext(side=2, outer=T, text="Relative change in Escapement")
  dev.off()
  
} # end plot escape change function


#*###########################################################
#  plotCatch
#######################################################
# Purpose: Dot and Line showing average Pre-terminal, Terminal, and Total catch (averaged over all years) 
#                for the 2 OMs and 6 MPs shown at April 2018 meetings
# Arguments: Names = names of scenarios to be plotted; each name must have an rds file with that name saved in DataOut folder
#            PlotName = name that will be used when saving plots as a pdf file in Figures folder 
# Note: Lots of hardwiring in this function; will need to be updated to be more generic
################################################################
plotCatch<-function(Names, PlotName,  
                    mpNames = c("CurrentF", "Fish50", "PT_Fish50", "NoFishing", "Productivity", "Capacity")) {

  
  # want to extract table of median change in escapment over 10 years (2016 to 2025) and 20 years (2016 to 2035)
  # also save all blobs
  Blobs <- list()
  Escape <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
  }
  
  PT_Catch_Medians <- NULL
  PT_Catch_SDs <- NULL
  Term_Catch_Medians <- NULL
  Term_Catch_SDs <- NULL
  Total_Catch_Medians <- NULL
  Total_Catch_SDs <- NULL
  
  NY <- Blobs[[1]]$Data$NY
  # which years to exclude?
  MaxAge <- Blobs[[1]]$Data$MaxAge

  # extract catches across all stocks and ages from each scenario
  for(mm in 1:length(Names)){
   
    # Extract catch from Blob
    Term_Catch_Sims <-  Blobs[[mm]]$Sims$Term_Catch
    PT_Catch_Sims <- Blobs[[mm]]$Sims$PT_Catch
    
    # sum over ages, then stocks ( how to do in one step?)
    Term_Catch_Summed_1 <- lapply(1:length(Term_Catch_Sims), function(x) apply(Term_Catch_Sims[[x]], 1:2,sum) )
    Term_Catch_Summed_2 <- lapply(1:length(Term_Catch_Summed_1), function(x) apply(Term_Catch_Summed_1[[x]], 1,sum) )
    # Now get mean across years
    Term_Catch_Means <- unlist(lapply(1:length(Term_Catch_Summed_2), function(x) mean(Term_Catch_Summed_2[[x]][(MaxAge+1):NY])))
    
    PT_Catch_Summed_1 <- lapply(1:length(PT_Catch_Sims), function(x) apply(PT_Catch_Sims[[x]], 1:2,sum) )
    PT_Catch_Summed_2 <- lapply(1:length(PT_Catch_Summed_1), function(x) apply(PT_Catch_Summed_1[[x]], 1,sum) )
    PT_Catch_Means <- unlist(lapply(1:length(PT_Catch_Summed_2), function(x) mean(PT_Catch_Summed_2[[x]][(MaxAge+1):NY])))

    
    # Now summarize as Medians and SDs
    Term_Catch_Medians[mm] <-  median(Term_Catch_Means)
    Term_Catch_SDs[mm] <- sd(Term_Catch_Means)
    
    PT_Catch_Medians[mm] <- median (PT_Catch_Means)
    PT_Catch_SDs[mm] <- sd(PT_Catch_Means)
    
    Total_Catch_Medians[mm] <-  median(PT_Catch_Means + Term_Catch_Means )
    Total_Catch_SDs[mm] <- sd(PT_Catch_Means + Term_Catch_Means )
    
  }
  
  pdf(paste("Figures/", PlotName ,".pdf", sep=""))
  
  par(mfrow=c(3,1),mar=c(1,3,1,1), oma=c(10,3,2,2))
  
  #colList<-c("black","red", "blue","steelblue2","purple","forestgreen","black","red", "blue","steelblue2","purple","forestgreen")
  
  # First plot, PT_Catch across scenerios
  # get ylims
  ylims <- c(0, max(PT_Catch_Medians+PT_Catch_SDs)*1.25)
  plot(1:length(Names), PT_Catch_Medians, pch=19, cex=1.2, ann=F, axes=F, xlim=c(0,(length(Names)+1)), ylim=ylims)
  for(mm in 1:length(Names)){
    arrows(mm,  (PT_Catch_Medians[mm] - PT_Catch_SDs[mm]), mm, (PT_Catch_Medians[mm] + PT_Catch_SDs[mm]),
           length=0) 
  }
  axis(2, at= c(0, 50000, 100000), labels=c("0", "50,000", "100,000"))
  abline(v= mean(c (0,(length(Names)+1) )))
  box()
  axis(1, at=1:length(Names), labels=FALSE)
  #text(x=1:12, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
       #labels=mpNames, srt=45, adj=1, xpd=TRUE)
  mtext(side=3, "Pre-Terminal Catch")
  text(x=c(3, (length(Names)/2 + 4)), y=par()$usr[4], labels=c("Average Productivity", "Recent Productivity"), pos=1, cex=1.2)
  
  # Now terminal
  # First plot, PT_Catch across scenerios
  # get ylims
  ylims <- c(0, max(Term_Catch_Medians+Term_Catch_SDs)*1.25)
  plot(1:length(Names), Term_Catch_Medians, pch=19, cex=1.2, ann=F, axes=F, xlim=c(0,(length(Names)+1)), ylim=ylims)
  for(mm in 1:length(Names)){
    arrows(mm,  (Term_Catch_Medians[mm] - Term_Catch_SDs[mm]), mm, (Term_Catch_Medians[mm] + Term_Catch_SDs[mm]),
           length=0) 
  }
  axis(2, at= c(0, 50000, 100000), labels=c("0", "50,000", "100,000"))
  abline(v= mean(c (0,(length(Names)+1) )))
  box()
  axis(1, at=1:length(Names), labels=FALSE)
  #text(x=1:12, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
       #labels=mpNames, srt=45, adj=1, xpd=TRUE)
  mtext(side=3, "Terminal Catch")
  
  
  # Now total Catch
  # First plot, PT_Catch across scenerios
  # get ylims
  ylims <- c(0, max(Total_Catch_Medians+Total_Catch_SDs)*1.25)
  plot(1:length(Names), Total_Catch_Medians, pch=19, cex=1.2, ann=F, axes=F, xlim=c(0,(length(Names)+1)), ylim=ylims)
  for(mm in 1:length(Names)){
    arrows(mm,  (Total_Catch_Medians[mm] - Total_Catch_SDs[mm]), mm, (Total_Catch_Medians[mm] + Total_Catch_SDs[mm]),
           length=0) 
  }
  axis(2, at= c(0, 100000, 200000), labels=c("0", "100,000", "200,000"))
  abline(v= mean(c (0,(length(Names)+1) )))
  box()
  axis(1, at=1:length(Names), labels=FALSE)
  text(x=1:length(Names), y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
       labels=mpNames, srt=45, adj=1, xpd=NA)
  mtext(side=3, "Total Catch")
 
  mtext(side=1, outer=T, "Scenario", line=5)
  mtext(side=2, outer=T, "Average Catch")
  
 dev.off()
  
} # end plotcatch function


#*###########################################################
#  checkNSims
#######################################################
# Purpose: Check number of simulation replicates needed to get stable results
#   
# Arguments: Names = names of scenarios to be plotted; each name must have an rds file with that name saved in DataOut folder
#            perfMetric = must specify which performance measuree to plot; at present, two options available: EscapeChange or AveCatch
# Note: # --> More testing is needed to make sure this works for current set-up
################################################################
checkNSims<-function(Names, perfMetric) {
  
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
      Esc_Start <- Esc_Init$Abundance[which(Esc_Init$StockID==Stocks[ss] & Esc_Init$Year %in% 2012:2015)]
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
  
  # For each stock, boxplot of 10 and 20 year escapement change for each # of sims, basic model
  Mods <- seq(1:length(Names))
  cols <- c("black", "blue", "darkgreen")
  
  
  if (perfMetric == "EscapeChange") {

    pdf("Figures/Escape_Change_Sims.pdf")
    par(mfrow=c(2,1), oma=c(3,2,2,2), mar=c(2,2,0,0))
    indx <- 1
    plot(0, ann=F, axes=F, xlim=c(0,34), ylim=c(0,10), type="n")
  
    for(ss in 1:length(Stocks)){
      for(mm in 1:length(Mods)){
        boxplot( EscChange_2025[[Mods[mm]]][[ss]], at=indx, add=T, axes=F, border=cols[mm])
        indx <- indx+1
      }
      indx <- indx+2
    } # end stocks loop
    axis(2)
    axis(1, at=seq(2,32, by=5), labels=F)
    legend("topleft", col=cols, lwd=2, legend=c("100", "500", "1000"), title="# of Sims", bty="n", cex=0.75)
    mtext(side=3, text="Init vs 2022:2025")
    # now 2035
    indx <- 1
    plot(0, ann=F, axes=F, xlim=c(0,34), ylim=c(0,10), type="n")
    for(ss in 1:length(Stocks)){
      for(mm in 1:length(Mods)){
        boxplot( EscChange_2035[[Mods[mm]]][[ss]], at=indx, add=T, axes=F, border=cols[mm])
        indx <- indx+1
      }
      indx <- indx+2
    } # end stocks loop
    axis(2)
    axis(1, at=seq(2,32, by=5), labels=F)
    # add labels
    text(x=seq(2,32, by=5), y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
       labels=Stocks, srt=45, adj=1, xpd=NA, cex=0.8)
    mtext(side=3, text="Init vs 2032:2035")
    mtext(side=1, outer=T, text="Stock", line=1.5)
    mtext(side=2, outer=T, text="Change in Escapement")
  
  dev.off() 
  } # end of if perMetric == "EscapeChange"
  
  
  if (perfMetric == "AveCatch") {
    # Same for avg catch
    pdf("Figures/Avg_Catch_Sims.pdf")
    par(mfrow=c(2,1), oma=c(2,2,2,5), mar=c(4,2,0,0))
    # get ylims for both
    MaxCatch <- max(c(unlist(AvgCatch), unlist(AvgCatch_Init)))
    indx <- 1
    plot(0, ann=F, axes=F, xlim=c(0,4), ylim=c(50000,MaxCatch), type="n")
    for(mm in 1:length(Mods)){
      boxplot( AvgCatch_Init[[Mods[mm]]], at=indx, add=T, axes=F, border=cols[mm])
      indx <- indx+1
    }
    axis(2)
    axis(1, at=1:3 , labels=c("100", "500", "1000"))
    box()
    mtext(side=3, text="Avg Catch 2016:2020")
    # now 2035
    indx <- 1
    plot(0, ann=F, axes=F, xlim=c(0,4), ylim=c(50000,MaxCatch), type="n")
    for(mm in 1:length(Mods)){
      boxplot( AvgCatch[[Mods[mm]]], at=indx, add=T, axes=F, border=cols[mm])
      indx <- indx+1
    }
    axis(2)
    axis(1, at=1:3 , labels=c("100", "500", "1000"))
    box()
    # add labels
    mtext(side=3, text="Avg Catch all Years")
    mtext(side=1, outer=T, text="Number of Sims")
    mtext(side=2, outer=T, text="Average Catch")
  
  
    dev.off()
  
  } # End of If perfMetric "AveCatch"

}

#*###########################################################
#  Catch_4_Panel
#######################################################
# Purpose: Plots catch trajectories across OMs and MP's
#Arguments: Names = names of scenarios to be plotted; each name must have an rds file with that name saved in DataOut folder
#            PlotName = name that will be used when saving plots as a pdf file in Figures folder 
#            Mods2Plot_ = which models are plotted in each panel, added for flexibility
#            Scenarios = which set of scenarios, added for felixbility
#            Fishery = which fishery are we plotting? PT or Term
################################################################

Catch_4_Panel <- function(Names, PlotName="Catch_4Panel",
                           Mods2Plot_1 = c(1,2,3,4), Mods2Plot_2 = c(1,5,6), Mods2Plot_3 = c(7,8,9,10), 
                           Mods2Plot_4 = c(7,11,12), Scenarios = "Fishing_Habitat", Fishery="PT", BiasNum) {
  
  # Read in blobs
  Blobs <- list()
  Catch <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    if(Fishery=="PT"){
      Catch[[mm]] <- Blobs[[mm]]$Sims$PT_Catch
    } else if(Fishery=="Term"){
       Catch[[mm]] <- Blobs[[mm]]$Sims$Term_Catch
    }
  }
  
  cols <- c("#0000ff", "#b22222", "#006400", "#FF9937","#7F3DCC",  "#808080")
  #  Blue, Green, Red, Orange, Purple, grey
  TGrey <- "#80808050"
  
  pdf(paste("Figures/", PlotName ,".pdf", sep=""))
  par(mfrow=c(2,2), oma=c(2,3,2,1), mar=c(2,2,2,1))
  
  Stocks <- Blobs[[1]]$Data$StocksInfo$StockID
  Years <- Blobs[[1]]$Options$Years
  NY <- Blobs[[1]]$Data$NY
  MaxAge <- Blobs[[1]]$Data$MaxAge
  
  for(ss in 1:length(Stocks)){
    MyMean <- list()
    MySD <- list()
    # first store all model data so can get ylims
    for(mm in 1:length(Blobs)){
      if(Fishery=="PT"){
        # reduce to list of array NY, NF_PT, NAges
        StockDat <- lapply(Catch[[mm]], "[",,,,ss, ) 
      } else if(Fishery=="Term"){
        # reduce to list of array NY, NAges
        StockDat <- lapply(Catch[[mm]], "[",,ss, ) 
      }
      # Now sum over fisheries and ages, not have list of 100 vectors of length NY
      MyDat <- lapply(1:length(StockDat), function(x) apply(StockDat[[x]], 1, sum) )
      # Not get mean and SD across model runs, for this stock
      MyMean[[mm]] <-  sapply(1:NY, function(x) mean(sapply(MyDat, "[", x) ))
      MySD[[mm]] <- sapply(1:NY, function(x) sd(sapply(MyDat, "[", x) ))
    } # end mod loop
    # remove initialization years
    Years_To_Plot <- Years[-(1:MaxAge)]
    MyMean <- lapply(MyMean, "[", which(Years %in% Years_To_Plot))
    MySD <- lapply(MySD, "[", which(Years %in% Years_To_Plot))
    # Get ylims
    # need wider bounds for some stocks to make room for legend
    if(Stocks[ss] %in% c("FRsp4.2", "FRsu4.1", "FRsu5.2", "FRL.N", "WCVI.H", "WC.N", "WC.H", "CBC", "USG" )){
      ylims <- c( 0 , max( unlist(MyMean) + unlist(MySD), na.rm=T)*1.4)
    } else {
      ylims <- c( 0 , max( unlist(MyMean) + unlist(MySD), na.rm=T))
    }

    # Plot
    #**************************************************
    # Panel 1 -- Basic Ricker with 4 fishing scenarios OR smolt release scenarios
    Cols <- c("#000000", cols[1:3])
    # start with leadin and inits
    plot(x=NA, y=NA, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Years_To_Plot), max(Years_To_Plot)), lwd=2, ann=F)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_1)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_1[mm]]], col=Cols[mm], lwd=2)
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_1[mm]]] - MySD[[Mods2Plot_1[mm]]], rev(MyMean[[Mods2Plot_1[mm]]] + MySD[[Mods2Plot_1[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
    } # end mod loop
    # add Lines for capacity, Smax

      if(Scenarios=="Fishing_Habitat"){
        legend("topleft", col=c(Cols), 
               legend=c("Current", "50% Curr.", "50% Curr. PT", "No Fishing"), title="Fishing Scenario",
               bty="n", cex=0.8, lwd=2, lty=c(1,1,1,1))
      } else if(Scenarios=="Bias_Smolts"){
        legend("topright", col=c(Cols), 
               legend=c("Current", "150% Curr.", "150% Curr. 2026"), title="Smolt Release Scenario",
               bty="n", cex=0.8, lwd=2, lty=c(1,1,1))
      }
    mtext(side=2, line=3, text="Historical Productivity")
    
    #**************************************************
    # Panel 2 -- Basic Ricker with 2 Productivity/CC scenarios OR AI bias scenarios
    Cols <- c("#000000", cols[4:5])
    # start with leadin and inits
    plot(x=NA, y=NA, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Years_To_Plot), max(Years_To_Plot)), lwd=2, ann=F)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_2)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_2[mm]]], col=Cols[mm], lwd=2)
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_2[mm]]] - MySD[[Mods2Plot_2[mm]]], rev(MyMean[[Mods2Plot_2[mm]]] + MySD[[Mods2Plot_2[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
    } # end mod loop

      if(Scenarios=="Fishing_Habitat"){
        legend("topleft", col=Cols, legend=c("No Change", "Prod. 25%", "CC 25%"  ), 
               title="Habitat Scenario", bty="n", cex=0.8, lwd=2, lty=c(1,1,1))
      } else if(Scenarios=="Bias_Smolts"){
        legend("topright", col=Cols, legend=c("No Bias", paste("+", BiasNum, "% AI Bias", sep=""), paste("-", BiasNum, "% Ai Bias", sep="")), 
               title="AI Bias Scenario", bty="n", cex=0.8, lwd=2, lty=c(1,1,1))
      }
    
    
    #**************************************************
    # Panel 3 --  R.B with 4 fishing scenarios
    Cols <- c("#000000", cols[1:3])
    # start with leadin and inits
    plot(x=NA, y=NA, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Years_To_Plot), max(Years_To_Plot)), lwd=2, ann=F)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_3)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_3[mm]]], col=Cols[mm], lwd=2)
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_3[mm]]] - MySD[[Mods2Plot_3[mm]]], rev(MyMean[[Mods2Plot_3[mm]]] + MySD[[Mods2Plot_3[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
    } # end mod loop
    mtext(side=2, line=3, text="Recent Productivity")
    
    #***************************************************
    # Panel 4 -- R.B with 2 productivity/CC scenarios
    Cols <- c("#000000", cols[4:5])
    # start with leadin and inits
    plot(x=NA, y=NA, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Years_To_Plot), max(Years_To_Plot)), lwd=2, ann=F)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_4)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_4[mm]]], col=Cols[mm], lwd=2)
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_4[mm]]] - MySD[[Mods2Plot_4[mm]]], rev(MyMean[[Mods2Plot_4[mm]]] + MySD[[Mods2Plot_4[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
      # add Smax/capacity lines
    } # end mod loop
    
    # label stock
    mtext(side=3, text=Stocks[ss], outer=T)
    # add overall labels
    mtext(side=1, text="Year", outer=T, line=0.5)
    mtext(side=2, text = paste(Fishery, "Catch"), outer=T)
  } # end Stock loop
  dev.off()
  
} # End catch 4 panel

#*###########################################################
#  AI_4_Panel
#######################################################
# Purpose: Plots catch trajectories across OMs and MP's
# Arguments: Names = names of scenarios to be plotted; each name must have an rds file with that name saved in DataOut folder
#            PlotName = name that will be used when saving plots as a pdf file in Figures folder 
#            Mods2Plot_ = which models are plotted in each panel, added for flexibility
#            Scenarios = which set of scenarios, added for flexibility
################################################################

AI_4_Panel <- function(Names, PlotName="AI_4Panel",
                          Mods2Plot_1 = c(1,2,3,4), Mods2Plot_2 = c(1,5,6), Mods2Plot_3 = c(7,8,9,10), 
                          Mods2Plot_4 = c(7,11,12), Scenarios = "Fishing_Habitat", BiasNums) {
  
  # Read in blobs
  Blobs <- list()
  AI <- list()
  for(mm in 1:length(Names)){
    Blobs[[mm]] <- readRDS(paste("DataOut/", Names[mm], ".rds", sep=""))
    AI[[mm]] <- Blobs[[mm]]$Sims$AI
  }
  
  Stocks <- Blobs[[1]]$Data$StocksInfo$StockID
  Years <- Blobs[[1]]$Options$Years
  NY <- Blobs[[1]]$Data$NY
  MaxAge <- Blobs[[1]]$Data$MaxAge
  
  cols <- c("#0000ff", "#b22222", "#006400", "#FF9937","#7F3DCC",  "#808080")
  #  Blue, Green, Red, Orange, Purple, grey
  TGrey <- "#80808050"
  
  MyMean <- list()
  MySD <- list()
  # first store all model data so can get ylims
  for(mm in 1:length(Blobs)){
    # Not get mean and SD across model runs
    MyMean[[mm]] <-  sapply(1:NY, function(x) mean(sapply(AI[[mm]], "[", x) ))
    MySD[[mm]] <- sapply(1:NY, function(x) sd(sapply(AI[[mm]], "[", x) ))
  } # end mod loop
  # remove initialization years
  Years_To_Plot <- Years[-(1:MaxAge)]
  MyMean <- lapply(MyMean, "[", which(Years %in% Years_To_Plot))
  MySD <- lapply(MySD, "[", which(Years %in% Years_To_Plot))
  
  pdf(paste("Figures/", PlotName ,".pdf", sep=""))
  par(mfrow=c(2,2), oma=c(2,3,2,1), mar=c(2,2,2,1))
    # Get ylims
      ylims <- c( 0 , max( unlist(MyMean) + unlist(MySD), na.rm=T))
    # Plot
    #**************************************************
    # Panel 1 -- Basic Ricker with 4 fishing scenarios OR smolt release scenarios
    Cols <- c("#000000", cols[1:3])
    # start with leadin and inits
    plot(x=NA, y=NA, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Years_To_Plot), max(Years_To_Plot)), lwd=2, ann=F)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_1)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_1[mm]]], col=Cols[mm], lwd=2)
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_1[mm]]] - MySD[[Mods2Plot_1[mm]]], rev(MyMean[[Mods2Plot_1[mm]]] + MySD[[Mods2Plot_1[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
    } # end mod loop
    # add Lines for capacity, Smax
    
    if(Scenarios=="Fishing_Habitat"){
      legend("topleft", col=c(Cols), 
             legend=c("Current", "50% Curr.", "50% Curr. PT", "No Fishing"), title="Fishing Scenario",
             bty="n", cex=0.8, lwd=2, lty=c(1,1,1,1))
    } else if(Scenarios=="Bias_Smolts"){
      legend("topright", col=c(Cols), 
             legend=c("Current", "150% Curr.", "150% Curr. 2026"), title="Smolt Release Scenario",
             bty="n", cex=0.8, lwd=2, lty=c(1,1,1))
    }
    mtext(side=2, line=3, text="Historical Productivity")
    
    #**************************************************
    # Panel 2 -- Basic Ricker with 2 Productivity/CC scenarios OR AI bias scenarios
    Cols <- c("#000000", rep(cols[4:6],each=2))
    ltys <- c(1,1,2,1,2,1,2)
    # start with leadin and inits
    plot(x=NA, y=NA, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Years_To_Plot), max(Years_To_Plot)), lwd=2, ann=F)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_2)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_2[mm]]], col=Cols[mm], lwd=2, lty=ltys[mm])
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_2[mm]]] - MySD[[Mods2Plot_2[mm]]], rev(MyMean[[Mods2Plot_2[mm]]] + MySD[[Mods2Plot_2[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
    } # end mod loop
    
    if(Scenarios=="Fishing_Habitat"){
      legend("topleft", col=Cols, legend=c("No Change", "Prod. 25%", "CC 25%"  ), 
             title="Habitat Scenario", bty="n", cex=0.8, lwd=2, lty=c(1,1,1))
    } else if(Scenarios=="Bias_Smolts"){
      legend("topright", col= c(Cols, Cols[2:3]), legend=c("No Bias", paste( "+", BiasNums ,"%",sep=""), paste( "-", BiasNums ,"%",sep="")  ), 
             title="AI Bias Scenario", bty="n", cex=0.8, lwd=2, lty=c(1,1,2,1,2,1,2))
    }
    
    
    #**************************************************
    # Panel 3 --  R.B with 4 fishing scenarios
    Cols <- c("#000000", cols[1:3])
    # start with leadin and inits
    plot(x=NA, y=NA, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Years_To_Plot), max(Years_To_Plot)), lwd=2, ann=F)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_3)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_3[mm]]], col=Cols[mm], lwd=2)
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_3[mm]]] - MySD[[Mods2Plot_3[mm]]], rev(MyMean[[Mods2Plot_3[mm]]] + MySD[[Mods2Plot_3[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
    } # end mod loop
    mtext(side=2, line=3, text="Recent Productivity")
    
    #***************************************************
    # Panel 4 -- R.B with 2 productivity/CC scenarios
    Cols <- c("#000000", rep(cols[4:6],each=2))
    # start with leadin and inits
    plot(x=NA, y=NA, col="darkgrey", type="l", ylim=ylims, xlim=c(min(Years_To_Plot), max(Years_To_Plot)), lwd=2, ann=F)
    # Now add model fits
    for(mm in 1:length(Mods2Plot_4)){
      lines(Years_To_Plot, MyMean[[Mods2Plot_4[mm]]], col=Cols[mm], lwd=2, lty=ltys[mm])
      # end mods loop
      # Add tranparent error bars
      polygon(y=c(MyMean[[Mods2Plot_4[mm]]] - MySD[[Mods2Plot_4[mm]]], rev(MyMean[[Mods2Plot_4[mm]]] + MySD[[Mods2Plot_4[mm]]])), 
              x=c(Years_To_Plot, rev(Years_To_Plot)), col= paste(Cols[mm], 50, sep=""), border=Cols[mm])
      # add Smax/capacity lines
    } # end mod loop
    
    # label stock
    mtext(side=3, text="AABM Abundance Index (AI)", outer=T)
    # add overall labels
    mtext(side=1, text="Year", outer=T, line=0.5)
    mtext(side=2, text = "AABM AI", outer=T)
  dev.off()
  
} # End catch 4 panel
