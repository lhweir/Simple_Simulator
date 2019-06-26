#########################################################
##  Functions for Chinook MSE Simulator version 1      ##
##            Brooke Davis & Kendra Holt               ##
#########################################################

# ========================================================================
# Functions included in this file:
# ========================================================================
# Run.CN.MSE.Sim ()        - main function that runs the simulation routine
# Distribute.Abundance ()  - distributes coastwide abundance among regions for a given time period in a given year
# Get.PT.Catch.Effort ()   - calculates catch, release mortality, and drop-off mortality for ISBM preterminal fisheries
# Get.Term.Catch ()        - calculates catch, release mortality, and drop-off mortality for terminal fisheries
# Get.Recruitment()        - takes escapement, returns age 1's (Natural, Hatchery, Enhanced)
# ========================================================================


###########################################################
#  Run.CN.Sim
#######################################################
# Purpose: Master function from which simulation replicates are run
# Arguments: Blob = list that contains all user-specified options and
#                     input data for a single scenario
################################################################
Run.CN.MSE.Sim <- function(Blob){

  # Make blob Data and options available in this environment so can access all Data Items
  Opts <- Blob$Options
  Data <- Blob$Data
  e <- environment()
  list2env(Opts, e)
  list2env(Data, e)

  # Create lists to store values in
  CohortList <- list()
  NList <- list()
  CatchList <- list()
  MatureList <- list()
  CTList <- list()
  EscList <- list()
  SpawnList <- list()
  TermFishList <- list()
  # PT relmort and dropoff (will only have this for AABM/Catch-controlled fisheries)
  RelMortList <- list()
  DropoffList <- list()
  # For ER fisheries only have Incidental mortality, not explicit Release and dropoff mort
  IncMortList <- list()
  TermIncMortList <- list()

  # extract survival for each stock
  #** note SR_1 never gets used
  Surv <- array(dim=c(NS, NAges+1))
  for(ss in 1:NS){
    grp <- StocksInfo$Surv_Grp[which(StocksInfo$StockID==Stocks[ss])]
    Surv[ss,] <- as.numeric(Grp_Surv[which(Grp_Surv$Grp==as.character(grp)),  paste("SR_", 1:MaxAge,sep="")])
  }

  # Start simulation loop
  for(sim in 1:nSims){

    # =================================================
    #Set up arrays to store values for each simulation
    # ==================================================

    # Stick to indexing order
    # Year | period | region | Fishery | stock | age

    #Cohort -- starting abundance before distribution,
    #represents "new" fish to model that need to be distributed
    # note Cohort is indexed by age (1-5), where N is indexed by adult age (2-MaxAge)
    Cohort <- array(dim=c(NY, NS, NAges+1))

    #abundance N -- indexes by stock, year, period, region, adult ages age (ages 2-MaxAge indexed as 1-(MaxAge-1))
    N <- array(dim=c(NY, NP, NR, NS, NAges))

    # have seperated for each step to make checks easier
    # after pre-terminal fishing
    N1 <- array(dim=c(NY, NP, NR, NS, NAges))
    # after incidental mortality
    N2 <- array(dim=c(NY, NP, NR, NS, NAges))
    # after maturity (only changes in period where stock matures)
    N3 <- array(dim=c(NY, NP, NR, NS, NAges))

    # FisheryID's (numeric) for pre-terminal fisheries
    ff_PT <- FisheryInfo$FisheryID[which(FisheryInfo$Type !="TERM")]
    # pre-terminal catch
    Catch <- array(dim=c( NY, NP, length(ff_PT), NS, NAges))

    # Keep terminal catch separate, since different indices
    Catch_Term <- array(dim=c(NY, NS, NAges))
    # also want to save by Terminal Fishery
    ff_Term <- FisheryInfo$FisheryID[which(FisheryInfo$Type=="TERM")]
    Term_Fisheries <- array(dim=c(NY, length(ff_Term)))

    # for Effort-based will only have Incidental mortality (ER PT (non-AABM), all terminal)
    IncMort <- array(dim=c(NY, NP, length(ff_PT), NS, NAges))
    IncMort_Term <- array(dim=c(NY, NS, NAges))
    # Releases (Pre-terminal AABM)
    Release <- array(dim=c(NY, NP, length(ff_PT), NS, NAges))
    # Release Mortality (Pre-terminal AABM)
    RelMort <- array(dim=c(NY, NP, length(ff_PT), NS, NAges))
    # DropOff Mortality (Pre-terminal AABM)
    Dropoff <- array(dim=c(NY, NP, length(ff_PT), NS, NAges))

    # Mature Fish
    Mature <- array(dim=c(NY, NR, NS, NAges))
    # Escapement (after terminal fishing)
    Escape <- array(dim=c(NY, NS, NAges))
    # Spawners
    Spawners <- array(dim=c(NY, NS))
    # Recruits
    Recruits <- array(dim=c(NY, NS))

    #=================================================================
    ## Start annual loops here
    #======================================================================
    # Cycle over years
    for(yy in 1:length(Years)){

      #===================================================
      ## Set Inital Cohort Values from inputs
      #===================================================

      ## Need to set up Cohort, if initializing from from known COhort values, just fill in year 1
      ## if from escapement will need to fill in zeros for years up to MaxAge
        if(Initialization=="Cohort"){
          if(yy==1){
            for(ss in 1:NS){
              for(aa in 2:MaxAge){
                # Ages in array will start from 2
                Cohort[yy,ss,aa] <- BaseCohort[which(BaseCohort$StockID==as.character(Stocks[ss]) & BaseCohort$Age==(aa))]
              } # end age loop
              # get age1 cohort abund from age 1
              Cohort[yy, ss, 1] <-  Cohort[yy, ss, 2] / Surv[ss,2]
            } # end stock loop
          } # end yy==1 if
        } else if(Initialization=="Escapement"){
          if(yy <= MaxAge){
            Cohort[yy,,yy:MaxAge] <- 0
          } # end if initialization period
        } # end escapement init else

      #=======================================
      # Start period loop
      #=======================================
      for(pp in 1:NP){

        #=======================================================================
        ## If Period 1 apply annual overwinter survival to last year's abundance
        #=======================================================================
        if (pp == 1) {

          # If this is the first period in year 2+, calculate cohort abundance as surviving adults from last year
          #   (i.e., N3 array from previous year * survival)
          if (yy > 1) {

            # Loop over stocks
            for(ss in 1:NS){
              #==============================================================================
              # If only one region, take straight from last year's post-maturation abundance
              #==============================================================================
                if (NR == 1) {
                  if (as.character(StocksInfo$Ocean_Stream[which(StocksInfo$StockID==Stocks[ss])]) == "Ocean") {
                    # Age 1 cohort abundance will come from Get.Recruitment function below
                    # Age 2 cohort abundance: need to fill using last year's age 1 cohort abundance and applying age 1:2 survival
                    Cohort[yy, ss, 2] <- Cohort[yy-1, ss, 1] * Surv[ss,2]
                    # Age 3+: apply survival to remaining at-sea abundance (after maturation occurs) in the
                    # last period of last year, from the single region
                    # Cohort(this year,  ages 3,4,5) = Post_maturity_Abund(last year, ages 2,3,4) * survival(ages: 2:3, 3:4, 4:5)
                    # K.Holt note: if ocean-type, want to get survival for ages 2:3, 3:4, and 4:5
                    Cohort[yy, ss, 3:MaxAge] <- N3[yy-1,NP,1,ss,1:3] * Surv[ss, 3:MaxAge]
                  } # end ocean if

                  if (as.character(StocksInfo$Ocean_Stream[which(StocksInfo$StockID==Stocks[ss])]) == "Stream") {
                    # Age 1 cohort abundance will be set to 0 below (fish are still in freshwater)
                    # Age 2 cohort abundance will come from Get.Recruitment function below
                    # Age 3+: apply survival to remaining at-sea abundance (after maturation occurs) in the
                    # last period of last year, from the single region
                    # if stream-type, want to get survival for ages 1:2, 2:3, 3:4 -- survival by ocean age
                    Cohort[yy, ss, 3:MaxAge] <- N3[yy-1,NP,1,ss,1:3] * Surv[ss,2:(MaxAge-1)]
                  } # end stream if

              #=====================================================
              #If more than one region need to sum N3 across regions
              #=====================================================
              } else if(NR > 1) {

                  if (as.character(StocksInfo$Ocean_Stream[which(StocksInfo$StockID==Stocks[ss])]) == "Ocean") {
                    # Age 1 cohort abundance will come from Get.Recruitment function below
                    # Age 2 cohort abundance: need to fill using last year's age 1 cohort abundance and applying age 1:2 survival
                    Cohort[yy, ss, 2] <- Cohort[yy-1, ss, 1] * S[2]
                    # Age 3+: apply survival to remaining at-sea abundance (after maturation occurs) in the last period of last year, summed over all regions
                    Cohort[yy, ss, 3:MaxAge] <- apply(N3[yy-1,NP,,ss,1:3], 2, sum) * S[3:MaxAge]
                  } # end oean if

                  if (as.character(StocksInfo$Ocean_Stream[which(StocksInfo$StockID==Stocks[ss])]) == "Stream") {
                    # Age 1 cohort abundance will be set to 0 below (fish are still in freshwater)
                    # Age 2 cohort abundance will come from Get.Recruitment function below
                    # Age 3+: apply survival to remaining at-sea abundance (after maturation occurs) in the last period of last year, summed over all regions
                    Cohort[yy, ss, 3:MaxAge] <- apply(N3[yy-1,NP,,ss,1:3], 2, sum) * S[2:(MaxAge-1)]
                  } # end stream if

                } # end more than 1 region if

            } # end of stock loop
          } # end of if yy > 1

          #===========================================================
          ## Distribute Cohort Abundandance Across Regions for Period 1
          #==============================================================
          N[yy,pp,,,] <- Distribute.Abundance(Abund=Cohort[yy,,2:MaxAge], RDist, yy, pp, Stocks, Ages, NR, NS, Region_Tab)

        #============================================================================================
        ## If Period > 1, need to redistribute abundance from last period to regions
        #=============================================================================================
        # Distribute abundance from last period, summed over regions (only matters if multiple periods)
        } else if(pp>1) {
          # get remaining fish from last period summed over regions
          N_remaining <- apply(N3[yy,pp-1,,,], 2:3, sum) # array dimensions NS X NAges
          # Cohort has age 1 index, need to add dummy so indexing aligns
          N[yy,pp,,,] <- Distribute.Abundance(Abund=N_remaining, RDist, yy, pp, Stocks, Ages, NR, NS, Region_Tab)
        } # end else if period > 1

        #=====================================================================
        ##  Preterminal Fishing
        #====================================================================
        #Loop over PT fisheries
        for(ff in ff_PT){

          # Extract fishery region
          region <- FisheryInfo$Region[which(FisheryInfo$FisheryID==ff)]

          #===========================================================
          # ER-based fisheries (also called ISBM/effort-controlled)
          #=============================================================
          ## If fishery is ISBM, calculate preterminal catch and incidental mortality
          # Also do this for initialization years for AABM fisheries, just to represent depletion without
          # being able to actually calculate AI's/TAC's etc.
          PT_Catch<-Get.PT.Catch.Effort(Years, Stocks, Ages, yy, pp,region, ff, N, Base_ER, Base_ER_Landed, HRS, FisheryInfo)

          # Now pull individual values from PT_Catch object
          Catch[yy, pp, ff, , ] <- round(PT_Catch$Catch)
          Release[yy,pp,ff, , ] <- round(PT_Catch$Release)
          RelMort[yy,pp,ff, , ] <- round(PT_Catch$RelMort)
          Dropoff[yy,pp,ff, , ] <- round(PT_Catch$Dropoff)
          IncMort[yy,pp,ff, , ] <- round(PT_Catch$IncMort)

        } # end of Preterminal fishery loop

        #==========================================
        # Remove PT mortality from population
        #==========================================
        # Go through regions, stocks ages, and remove fishing, incidental mort.
        for(rr in 1:NR){
          for(ss in 1:NS){
            for(aa in 1:NAges){
              # extract vector of pre-terminal fisheries in region
              ff_in_rr<-FisheryInfo$FisheryID[which(FisheryInfo$Region==rr & FisheryInfo$Type != "TERM")]
              # Remove fishing over all fisheries in region (ff_in_rr)
              N1[ yy, pp, rr, ss, aa] <- N[yy, pp, rr, ss, aa] - sum(Catch[yy, pp, ff_in_rr, ss, aa], na.rm=T)
              # Remove incidental mortality from fisheries in region (ff_in_rr)
              N2[ yy, pp, rr, ss, aa] <- N1[yy, pp, rr, ss, aa] - sum(IncMort[yy, pp, ff_in_rr, ss, aa], na.rm=T)
            } # end Age loop
          } # end stock loop
        } # end region loop

        #======================================================
        ## Maturation
        #======================================================
        # Check which stocks mature in this period, if any
        ss_mature <-  StocksInfo$StockID[which(StocksInfo$Maturation_Period==pp)]
        for(ss in 1:NS){
          if(Stocks[ss] %in% ss_mature){
            # which years do we have maturation rates for?
            DatYrs <- unique(Maturation$Year[which(Maturation$StockID==as.character(Stocks[ss]))])
            # Pull out most recent of years we have values for in data frame
            Mat_Tab <- Maturation[which(Maturation$Year==max(DatYrs[DatYrs<=Years[yy]]) & Maturation$StockID==as.character(Stocks[ss])),]
            # get randomized M vector for that stock
            Means <- sapply( Ages, function(x) Mat_Tab$Maturation[which(Mat_Tab$Age==x)] )
            SDs <- sapply( Ages, function(x) Mat_Tab$SD[which(Mat_Tab$Age==x)] )
            MatRates <- Trunc.Norm(mean=Means, SD=SDs, min=0, max=1)
            # for each region calc mature fish and remove from pop
            for(rr in 1:NR){
              Mature[yy, rr, ss, ] <-  N2[yy, pp, rr ,ss, ] * MatRates
              # Now remove these from population
              N3[yy, pp, rr, ss, ] <- N2[yy, pp, rr, ss,  ]   -  Mature[yy, rr,ss, ]
            } # end region loop
          # if stock not maturing in this period just move abund forward
          } else {
            N3[yy, pp, , ss, ] <- N2[yy, pp, , ss,  ]
          } # end not maturing stock else
        } # end stock loop
      #=====================
      } # End Period loop #
      #====================

      #====================================================
      ##  Terminal Fishing & Escapement & Recruitment  ##
      #====================================================
      # Which fisheries are terminal?
      ff_Term <- FisheryInfo$FisheryID[which(FisheryInfo$Type=="TERM")]
      #============================================
      ## Do individually for each stock
      #==============================================
      for(ss in 1:NS){

        # check that there is a terminal fishery associated with this stock
        if(Stocks[ss] %in% FisheryInfo$Term_Stock){
          # which terminal fisheries target this stock?
          ff_Stock <- FisheryInfo$FisheryID[which(FisheryInfo$Term_Stock==Stocks[ss])]
          # Create list to store terminal catch values in
          Term_Catch <- list()
          for(ff in 1:length(ff_Stock)){
            # extract HR for age classes as vector
            HR <- sapply ( Ages, function(x) TermHR[which(TermHR$FisheryID==ff_Stock[ff]), paste("HR_Age", x, sep="")])
            HR_Landed <- sapply ( Ages, function(x) TermHR_Landed[which(TermHR_Landed$FisheryID==ff_Stock[ff]), paste("HR_Age", x, sep="")])
            # Which years do we have Harvest rate scalars for?
            DatYrs <- HRS$Year[which(HRS$FisheryID==ff_Stock[ff])]
            # Get HRS from HRS file, most recent year
            HRScalar <- HRS$HRS[which(HRS$FisheryID==ff_Stock[ff]  & HRS$Year==max(DatYrs[DatYrs<=Years[yy]]))]

            # Get Catch and incidental mortality
            Term_Catch[[ff]] <- Get.Term.Catch(HR, HR_Landed, HRScalar,yy,ss, ff=ff_Stock[ff], NP, Mature)

            # also want to save catch for each fishery
            Term_Fisheries[yy, which(ff_Term==ff_Stock[ff])] <- sum(Term_Catch[[ff]]$Catch)
          } # end ff_Stock Loop

          # Extract results and fill catch arrays
          Catch_Term[yy, ss , ] <- rowSums(sapply(Term_Catch, function(x) x$Catch))
          IncMort_Term[yy,ss , ] <- rowSums(sapply(Term_Catch, function(x) x$IncMort))

        } else { # else if not terminal fishery for this stock
          Catch_Term[yy, ss, ] <- IncMort_Term[yy,ss , ]  <-  0
          Term_Fisheries[yy, which(ff_Term==ff_Stock[ff])] <- 0
        } # end terminal fishery if for that stock

        #==========================================================
        ##  Escapement
        #==========================================================
        # Calc. escapement by summing run by stock over regions
        # and removing catch, incidental mortality

        # Get total run by mature fish over regions
        # if only one region need different indices to take sum
        if(NR > 1){
          Total_Run <- apply(Mature[yy, ,ss,], 2, sum)
        } else if(NR==1){
          Total_Run <- Mature[yy,rr,ss,]
        }

        # Remove Terminal fishing mort from total run to get escapement
         Escape[yy, ss, ] <- Total_Run - Catch_Term[yy,ss, ] -  IncMort_Term[yy, ss, ]

        # If initializing by Escapement, and in initialization years, manually fill in Spawners
        # From input data
        if(Initialization=="Escapement" & yy <= MaxAge){
          Spawners[yy, ss] <- Escapement$Abundance[which(Escapement$Year==Years[yy] & Escapement$StockID==Stocks[ss])]
        } else {
          # if excluding Jacks, take out of spawner numbers
          if(Exclude_Jacks == F){
             Spawners[yy, ss] <- sum(Escape[yy, ss, ])
          } else {
            Spawners[yy, ss] <- sum(Escape[yy, ss, 2:NAges ])
          }
        }

        #========================================================
        #  Recruitment to Age 1
        #===========================================================

        # if isn't last year project forward to next year
        if(yy < length(Years)){
          # Get appropriate EV value
          # if drawing from 10 year mean, sd
          if(EV_Type=="10yrmean"){
            EV_mu <- EVs$Mean[which(EVs$StockID==Stocks[ss])]
            EV_sd <- EVs$SD[which(EVs$StockID==Stocks[ss])]
            EV <- rnorm(1, EV_mu, EV_sd)
          # or no EV being used
          } else if(EV_Type=="1"){
            EV <- 1
          } else if(EV_Type=="Annual"){
            # if don't have EV for that year use latest year
            DatYrs <- EVs$Year[which(EVs$StockID==as.character(Stocks[ss]))]
            EV <- EVs$EV[which(EVs$StockID==as.character(Stocks[ss]) & EVs$Year==max(DatYrs[DatYrs<=Years[yy+1]]))]
          }# end EV if's

          # If ocean-type, age 1 recruitment is predicted for year yy+1:
          if (as.character(StocksInfo$Ocean_Stream[which(StocksInfo$StockID==Stocks[ss])]) == "Ocean") {
            Cohort[yy+1, ss, 1] <- Get.Recruitment(Opts, Data, Years, Stocks, Ages, MaxAge, ss, yy, Var=T, Spawners, Maturation,
                                                   StocksInfo, Surv, Prod_Mults, CC_Mults, Smolts_Scenario) * EV

          # If stream-type, age 2 recruitment is predicted for year yy+2
          } else if (as.character(StocksInfo$Ocean_Stream[which(StocksInfo$StockID==Stocks[ss])]) == "Stream") {
            # Set age 1 cohort abundance to zero
            Cohort[yy+1, ss, 1]<- 0
            # Calculate recruitment to age 2 in year yy +2
            if (yy < (length(Years)-1)) {
              # currently not randomized -- using mean instead of S for survival
              Cohort[yy+2, ss, 2] <- Get.Recruitment(Opts, Data, Years, Stocks, Ages, MaxAge, ss, yy, Var=T, Spawners, Maturation,
                                                     StocksInfo, Surv, Prod_Mults, CC_Mults, Smolts_Scenario) * EV
            }
          } # end stream-type else
        } # end if not final year -- if final year nothing needs to happen

      #==================
      } # end stock loop
      #==================

    #===================
    } # End Year loop
    #===================

    # Put run values into lists
    CohortList[[sim]] <- Cohort
    NList[[sim]] <- N
    CatchList[[sim]] <- Catch
    MatureList[[sim]] <- Mature
    CTList[[sim]] <- Catch_Term
    EscList[[sim]] <- Escape
    SpawnList[[sim]] <- Spawners
    TermFishList[[sim]] <- Term_Fisheries
    RelMortList[[sim]] <- RelMort
    DropoffList[[sim]] <- Dropoff
    IncMortList[[sim]] <- IncMort
    TermIncMortList[[sim]] <- IncMort_Term


  #================
  }  #end sim loop
  #================

  #================================
  # Compile outputs
  #================================

  Outputs <- list()

  Outputs$Cohort <- CohortList
  Outputs$N <- NList
  Outputs$PT_Catch <- CatchList
  Outputs$Mature <- MatureList
  Outputs$Term_Catch <- CTList
  Outputs$Escape <- EscList
  Outputs$Spawners <- SpawnList
  # Same terminal catch data, just summarized differently
  Outputs$Term_Fisheries <- TermFishList
  Outputs$RelMort <- RelMortList
  Outputs$Dropoff <- DropoffList
  Outputs$IncMort <- IncMortList
  Outputs$IncMort_Term <- TermIncMortList

  Outputs

} # end DGM simulator Function


###########################################################
#  Distribute.Abundance
#######################################################
# Purpose: Allocates stock- and age-specific cohort abundance among regions for a given period
# Arguments: Abund = Abundance array to be distributed across regions
#            RDist = distribution coefficient data frame
#            yy = year
#            pp = period
#            Stocks = vector of stock names
#            Ages = vector of ages (2:MaxAge)
#            NR = number of regions
#            NS = number of stocks
#            Region_Tab = Data table with Region names and indices
################################################################
Distribute.Abundance <- function (Abund, RDist, yy, pp, Stocks, Ages, NR, NS, Region_Tab) {

  #Cohort is matrix with dimensions NS x NAges
  #Initiate N_out array to output
  N_out <- array(dim=c(NR, NS, length(Ages)))

  # Distribute cohort
    for(ss in 1:NS){
      for(aa in 1:length(Ages)){
        # get vector of distributions (length NR)
        Dist_Vec <- RDist[which(RDist$StockID==as.character(Stocks[ss]) & RDist$Period==pp & RDist$Age == Ages[aa]),
                      which(names(RDist) %in% Region_Tab$Region_Name)]
        # use this dist vector to distribute abund of this stock and age
        N_out[,ss,aa] <- rmultinom(n=1, size= Abund[ss,aa], prob=Dist_Vec)
      } # end age loop
    } # end stock loop
  return(N_out)
} # end distribute.cohort function


###########################################################
#  Get.PT.Catch.Effort
#######################################################
# Purpose: Calculates catch, release mortality, and drop-off mortality for effort-controlled preterminal fisheries
# Arguments: Years = Vector of years being simulated
#            Stocks = Vector of stock names
#            Ages = Vector of ages (2:MaxAge)
#            yy = current year
#            pp = current period
#            region = region fishery takes place in
#            ff = fishery for which catch is to be calculated
#            eta = sub-legal release mortality rate
#            theta = drop-off mortality rate
#            N = abundance array, with abundance distributed among regions. Index order = [year, period, region, stock, age]
#            Base_ER = Input data frame of base ER's
#            HRS = harvest rate scalars to scale base ER for different scenarios
#            FisheryInfo = Input data frame with Fishery details
################################################################
Get.PT.Catch.Effort<-function(Years, Stocks, Ages, yy, pp, region, ff, N, Base_ER, Base_ER_Landed, HRS, FisheryInfo) {

  # Extract fishery havest rate scalar (HRS, or delta)
  # -- first, need to find closest lower year with HRS specified via input files
  DatYrs <-  HRS$Year[which(HRS$FisheryID==ff)]
  # -- then, choose largest of the years from same or previous years
  delta <- HRS$HRS[which(HRS$FisheryID==ff & HRS$Year==max(DatYrs[DatYrs<=Years[yy]]))]

  # extract matrix of ER values for all stocks X ages
  ER <- as.matrix(Base_ER[which(Base_ER$Period==pp & Base_ER$FisheryID==ff & Base_ER$StockID %in% Stocks),
                           which(names(Base_ER) %in% paste("ER_Age", Ages, sep=""))]  *  delta)

  # also get ER for landed catch
  ER_Landed <- as.matrix(Base_ER_Landed[which(Base_ER_Landed$Period==pp & Base_ER_Landed$FisheryID==ff &
                                                Base_ER_Landed$StockID %in% Stocks),
                                 which(names(Base_ER_Landed) %in% paste("ER_Age", Ages, sep=""))]  *  delta)

  # Get total mortality
  TotalMort_Mat <- N[yy,pp,region,,] * ER
  # Get landed catch
  Catch_Mat <-N[yy,pp,region,,] * ER_Landed
  # Get Incidental Mortality
  Inc_Mat <- TotalMort_Mat - Catch_Mat

  # No longer have explicit releases, relmort, dropoff, but will leave structure to consistant with AI-based
  CatchValues <- list()
  CatchValues$Catch <- Catch_Mat
  CatchValues$Release <- NA
  CatchValues$RelMort <- NA
  CatchValues$Dropoff <- NA
  CatchValues$IncMort <- Inc_Mat

  return(CatchValues)

} # end Get.PT.Catch.Effort function


###########################################################
#  Get.Term.Catch
#######################################################
# Purpose: Calculates catch, release mortality, and drop-off mortality for stock-specific terminal fisheries
# Arguments:


################################################################

Get.Term.Catch <-  function(HR, HR_Landed, HRScalar,yy,ss, ff=ff_Stock[ff], NP, Mature){

  # Get catch simply by multiplying by harvest rates to mature summed across regions
  HR_Tot <- HR*HRScalar
  HR_Landed <- HR_Landed*HRScalar
  # only works if more than one period
  if(NP > 1){
    Catch_Vec <- apply(Mature[yy, ,ss,], 2, sum) * HR_Landed
    Total_Mort_Vec <- apply(Mature[yy, ,ss,], 2, sum) * HR_Tot
  } else if(NP==1) {
    Catch_Vec <- Mature[yy, ,ss,] * HR_Landed
    Total_Mort_Vec <- Mature[yy, ,ss,] * HR_Tot
  }

  # Prepare outputs
  CatchValues <- list()
  CatchValues$Catch <- round(Catch_Vec)
  CatchValues$IncMort <- round(Total_Mort_Vec - Catch_Vec)

  return(CatchValues)
}



###########################################################
#  Get.Recruitment
#######################################################
# Purpose: Calculates age 1 recruitment based on spawner abundance
# Arguments:


################################################################

Get.Recruitment <- function(Opts, Data, Years, Stocks, Ages, MaxAge, ss, yy, Var, Spawners, Maturation, StocksInfo, Surv,
                             Prod_Mults, CC_Mults, Smolts_Scenario){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Case 1: Natural stock
  #~~~~~~~~~~~~~~~~~~~~~~~~~~
  # If compute SR = 0 & Truncate = 1, then use true Ricker relationship with Ricker b parameter from the input file
  if (StocksInfo$Type[which(StocksInfo$StockID==Stocks[ss])] == "Natural") {
      a <- StocksInfo$Ricker_A[which(StocksInfo$StockID== as.character(Stocks[ss]))] * Prod_Mults$Mult[which(Prod_Mults$Year == Years[yy])]
      # extract b from input data
      b <- StocksInfo$Ricker_B[which(StocksInfo$StockID==Stocks[ss])] / CC_Mults$Mult[which(CC_Mults$Year == Years[yy])]
      # calculate recruitment
      R <- a*Spawners[yy,ss] * exp(-b*Spawners[yy,ss])
    # Add Variability
    if(Var==T){
      # get process error (precision) from StockInfo file
      TauR <- StocksInfo$TauR[which(StocksInfo$StockID==Stocks[ss])]
      SD_R <- 1/sqrt(TauR)
      # Draw random variability
      E_R <- rnorm(1, 0, SD_R)
      # Multiply Recruitment by this value
      R <- R*exp(E_R)
    }
    # Convert to age 1's
    # get maturation rates for most recent year -- allows for both annual and one mat rate
    DatYrs <- unique(Maturation$Year[which(Maturation$StockID==as.character(Stocks[ss]))])
    Mat_Tab <- Maturation[which(Maturation$Year==max(DatYrs[DatYrs<=Years[yy]])),]
    # get M vector for that stock, ages included in this sim
    M <- c( Mat_Tab$Maturation[which(Mat_Tab$StockID==as.character(Stocks[ss]) & Mat_Tab$Age %in% Ages)])

    # get survival rate
    # Survival rate depends if stream or ocean type
    Stock_Type <- StocksInfo$Ocean_Stream[which(StocksInfo$StockID==Stocks[ss])]

    Age1Surv <- Get.Age1.Surv(S=Surv[ss,2:MaxAge] , M, MaxAge, Stock_Type)

    age1Fish <- R / Age1Surv

  } # end natural stock

  #~~~~~~~~~~~~~~~~~~~~~~~~
  # Case 2: Hatchery Stock
  #~~~~~~~~~~~~~~~~~~~~~~~~
  if (StocksInfo$Type[which(StocksInfo$StockID==Stocks[ss])] == "Hatchery") {

    # apply smot release scenario
    if(Opts$Smolts_Scenario != 1){
      # extract smolt change
      # get most recent year
      DatYrs <- unique(Data$Smolts$Year)
      Smolts_mult <- Data$Smolts$Scalar[which(Data$Smolts$Year==max(DatYrs[DatYrs<=Years[yy]]) &
                                                Data$Smolts$Stock==as.character(Stocks[ss]) )]
    } else {
      Smolts_mult <- 1
    }

    # Note, by my calculations, this value would be 22596335 for re-calibrated WCVI.H (with A = 5.52)
    smTarg <- StocksInfo$SmRelTarget[which(StocksInfo$StockID==Stocks[ss])] * Smolts_mult

    # converts smolt releases into age-1 recruits
    smSurvRate <- StocksInfo$SmSurvRate[which(StocksInfo$StockID==Stocks[ss])]

    # Calculate required brood take to meet smolt target
    # read in sex ratio and smolts per female
    smPerFemBT<-StocksInfo$SmPerFemBT[which(StocksInfo$StockID==Stocks[ss])]
    BTSexRatio<-StocksInfo$BTSexRatio[which(StocksInfo$StockID==Stocks[ss])]
    #  BT.req <- (smTarg*smSurvRate) / exp(enhProd) old equation
    BT_req <- smTarg/(smPerFemBT*BTSexRatio)

    # (note that smoltSurvRate is needed to convert units from smolts (used for smoltTarg)
    # to ocean age-1 recruits
    #age1Fish  <- smTarg * smSurvRate
    age1Fish <- ifelse( Spawners[yy, ss] < BT_req , Spawners[yy,ss]*smPerFemBT*BTSexRatio*smSurvRate, smTarg * SmoltSurv)

  } # end hatchery

  #~~~~~~~~~~~~~~~~~~~~~~~~~
  # Case 3: Enhanaced Stock
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  if (StocksInfo$Type[which(StocksInfo$StockID==Stocks[ss])] == "Enhanced") {


    # For natural component:
    a <- StocksInfo$Ricker_A[which(StocksInfo$StockID== as.character(Stocks[ss]))] * Prod_Mults$Mult[which(Prod_Mults$Year == Years[yy])]
    # extract b from input data
    b <- StocksInfo$Ricker_B[which(StocksInfo$StockID==Stocks[ss])] / CC_Mults$Mult[which(CC_Mults$Year == Years[yy])]

    # apply smot release scenario
    if(Opts$Smolts_Scenario != 1){
      # extract smolt change
      # get most recent year
      DatYrs <- unique(Data$Smolts$Year)
      Smolts_mult <- Data$Smolts$Scalar[which(Data$Smolts$Year==max(DatYrs[DatYrs<=Years[yy]]) &
                                                Data$Smolts$Stock==as.character(Stocks[ss]) )]
    } else {
      Smolts_mult <- 1
    }

    # For hatchery component:
    smTarg <- StocksInfo$SmRelTarget[which(StocksInfo$StockID==Stocks[ss])] * Smolts_mult
    smSurvRate <- StocksInfo$SmSurvRate[which(StocksInfo$StockID==Stocks[ss])]
    smPerFemBT<-StocksInfo$SmPerFemBT[which(StocksInfo$StockID==Stocks[ss])]
    BTSexRatio<-StocksInfo$BTSexRatio[which(StocksInfo$StockID==Stocks[ss])]
    maxBT<-StocksInfo$maxBT[which(StocksInfo$StockID==Stocks[ss])]

    # adjust spawners for hatchery broodtake
    BT_req <- smTarg/(smPerFemBT*BTSexRatio)
    BT <- ifelse(BT_req <= (Spawners[yy,ss]*maxBT), BT_req, (Spawners[yy,ss]*maxBT))
    Spawners[yy,ss] <- Spawners[yy,ss]-BT

    # calculate recruitment
    R <- a*Spawners[yy,ss] * exp(-b*Spawners[yy,ss])
    # Add Variability
    if(Var==T){
      # get process error (precision) from StockInfo file
      TauR <- StocksInfo$TauR[which(StocksInfo$StockID==Stocks[ss])]
      SD_R <- 1/sqrt(TauR)
      # Draw random variability
      E_R <- rnorm(1, 0, SD_R)
      # Multiply Recruitment by this value
      R <- R*exp(E_R)
    }
    # Convert to age 1's
    # get maturation rates for most recent year -- allows for both annual and one mat rate
    DatYrs <- unique(Maturation$Year[which(Maturation$StockID==as.character(Stocks[ss]))])
    Mat_Tab <- Maturation[which(Maturation$Year==max(DatYrs[DatYrs<=Years[yy]])),]
    # get M vector for that stock, ages included in this sim
    M <- c( Mat_Tab$Maturation[which(Mat_Tab$StockID==as.character(Stocks[ss]) & Mat_Tab$Age %in% Ages)])

    # get survival rate
    # Survival rate depends if stream or ocean type
    Stock_Type <- StocksInfo$Ocean_Stream[which(StocksInfo$StockID==Stocks[ss])]

    Age1Surv <- Get.Age1.Surv(S=Surv[ss, 2:MaxAge] , M, MaxAge, Stock_Type)

    age1Fish_natural <- R / Age1Surv

    # For hatchery component:

    #age1Fish_hatchery  <- smTarg * smSurvRate
    # BD changed 10/31/2018
    age1Fish_hatchery  <- BT*smPerFemBT*BTSexRatio*smSurvRate


    # Combine natural and hatchery age 1 production ======================================
    age1Fish <- age1Fish_hatchery + age1Fish_natural

  } # end enhanced


  age1Fish

} # end Get.Recruitment function








