#########################################################
##     Helper Functions for Chinook MSE Simulator      ##
##            Brooke Davis & Kendra Holt               ##
#########################################################


# ========================================================================
# Functions included in this file:
# ========================================================================
# %notin%       - adaptation of dplyr %in% functionignore style guidelines of "."'s separating words for func names to matchdplyr %in% style
# Trunc.Norm    - truncated normal distribution used for randomizing parameters, can take single values or vectors
# Get.Age1.Surv - Age one to adult survival used in S-R calculation
# ========================================================================



###########
## notin ##
###########
`%notin%` <- function(x,y) !(x %in% y)



###################
## Trunc.Norm   ###
###################

Trunc.Norm <- function(mean, SD, min=NA, max=NA){
  xx <- rnorm(length(mean), mean, SD)
  if(is.na(min)==F){
    xx <- ifelse(xx < min, min, xx)
  }
  if(is.na(max)==F){
    xx <- ifelse(xx > max, max, xx)
  }
  xx
}


#####################
##  Get.Age1.Surv  ##
#####################

# Given survival and maturation, calc age 1 to adult
Get.Age1.Surv <- function(S, M, MaxAge, Stock_Type){
  # S is vector of survival S_12, S_23, S_34, S_45
  # M is vector of Maturations starting at age 2 -- M_2, M_3, M_4, M_5

  if(Stock_Type=="Ocean"){
    SEQ <- NULL
    SEQ[1] <- S[1]
    for(i in 2:(MaxAge-1)){
      SEQ[i] <- SEQ[i-1]*(1-M[i-1])*S[i]
    } # end age loop
  } else if(Stock_Type=="Stream"){
    # for stream-type stocks
    SEQ <- NULL
    # STream-type fish will not mature until total age 3
    SEQ[1] <- 0
    SEQ[2] <- S[1]
    for(i in 3:(MaxAge-1)){
      SEQ[i] <- SEQ[i-1]*(1-M[i-1])*S[i-1]
    } # end age loop

  } # end stream type else
  # multiply each element of SEQ vector with vector of maturation rates

   sum(SEQ*M)


} # end Get.Age1.Surv function






