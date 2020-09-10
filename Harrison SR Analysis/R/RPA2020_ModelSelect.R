#=================================================================
#Analysis of Harrison recruitment data
#Author: Catarina Wor and Brooke Davis
#Date: August 2020
#=================================================================

#load in required packages and libraries

library(ggplot2)
library(TMB)
library(reshape)
library(dplyr)
library(TMBhelper)

source("TMB_functions.R")

#read in simple data set
SR <- read.csv("../data/Harrison_simples_Apr18.csv")

# Compile and load all TMB models before load tmbstan -- sometimes causes problems
# Simple Ricker
#dyn.unload(dynlib("Ricker_simple"))
compile("../R/Ricker_simple.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="")
dyn.load(dynlib("Ricker_simple"))

# Autocorr model
#dyn.unload(dynlib("../R/Ricker_autocorr_ch"))
compile("../R/Ricker_autocorr_ch.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="",tracesweep = TRUE)
dyn.load(dynlib("../R/Ricker_autocorr_ch"))

# TV alpha model
#dyn.unload(dynlib("../R/Rickerkf_ratiovar"))
compile("../R/Rickerkf_ratiovar.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="",tracesweep = TRUE)
dyn.load(dynlib("../R/Rickerkf_ratiovar"))

# load this after compile, since can mess with compiling for some reason
library(tmbstan)

#=================================================================
#Simple Ricker model

iteracs=100000

# LM version
srm <- lm(log(SR$R/SR$S_adj)~ SR$S_adj)
a_srm <- srm$coefficients[1]
b_srm <- -srm$coefficients[2]

summary(srm)

par(mfrow=c(2,2))
plot(srm)

#==========================================================
# TMB simple Ricker
# Scale Data so that it is close to 0, which helps with convergence in TMB
Scale <- 100000

mydata<-list(obs_R=SR$R/Scale,obs_S=SR$S_adj/Scale)
parameters_simple <- list(
  a=(a_srm),
  logbeta = log(b_srm*Scale),
  logSigObs= log(.4)  )

# also bound beta so that Smax is between 1 Spawner and 4*maximum observed
Min_S_Obs <- 1/Scale
Max_S_Obs <- max(SR$S_adj)/Scale
Min_LogB <- log(1/(Max_S_Obs*4))
Max_LogB <- log(1/(Min_S_Obs))

# create vector of bounds
#           a      logbeta   logSigObs
lower <- c(-Inf,   Min_LogB, -Inf)
upper <- c(Inf,    Max_LogB,  Inf)


obj<-MakeADFun(mydata,parameters_simple,random=NULL,DLL="Ricker_simple")
newtonOption(obj, smartsearch=FALSE)
  
opt<-nlminb(obj$par,obj$fn,obj$gr,  lower=lower, upper = upper)

# pull out estimates from ML fit
# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
All_Ests

# Add column with upper and lower CIs
All_Ests <- All_Ests %>% mutate(Lower = Estimate-1.96*Std..Error, Upper = Estimate+1.96*Std..Error)
All_Ests[All_Ests$Param == "beta", c("Estimate", "Std..Error", "Lower", "Upper")] <- All_Ests[ All_Ests$Param == "beta", c("Estimate", "Std..Error", "Lower", "Upper")]/Scale
All_Ests[All_Ests$Param == "Smax", c("Estimate", "Std..Error", "Lower", "Upper")] <- All_Ests[ All_Ests$Param == "Smax", c("Estimate", "Std..Error", "Lower", "Upper")]*Scale                      

#diagnosticts
qqnorm(obj$report()$residuals)
abline(0,1)
# doesn't look great

# Look at residuals
SRdiagsimple <- SR
SRdiagsimple$residuals <- obj$report()$residuals
SRdiagsimple$std_resids <- obj$report()$std_resids
SRdiagsimple$model <- "Simple Ricker"
SRdiagsimple$alpha <- obj$report()$alpha

rp <- ggplot(SRdiagsimple)
rp <- rp + geom_line(aes(x=BroodYear,y=residuals),size=1.2)
rp <- rp + geom_point(aes(x=BroodYear,y=residuals, col=residuals),stroke=3)
rp <- rp + geom_hline(yintercept = 0)
rp <- rp + theme_bw(16)
rp <- rp + scale_color_viridis_c(end = 0.8)
rp
# clear downwards pattern

# also look at std resids

ggplot(SRdiagsimple) + 
  geom_line(aes(x=BroodYear,y=std_resids),size=1.2) +
  geom_point(aes(x=BroodYear,y=std_resids, col=std_resids),stroke=3) +
  geom_hline(yintercept = 0) +
  theme_bw(16) + 
  scale_color_viridis_c(end = 0.8)

# a few residuals close to 2

#MCMC
simpleB<-list(
  obj=obj,
  nchain=3,
  iter=iteracs,
  lowbd = lower,  
  hibd= upper
)

posterior_simple <- posteriorsdf(simpleB)

#posteriors of derived quantities 
simpdf <- posterior_simple$posteriors

# switch parameters to match notation (note that alpha in TMB models is really a=log(alpha))
a <- (simpdf$value[simpdf$parameters=="a"])
alpha <- exp(simpdf$value[simpdf$parameters=="a"])
beta <- exp(simpdf$value[simpdf$parameters=="logbeta"])/Scale
Smax <- 1/exp(simpdf$value[simpdf$parameters=="logbeta"])*Scale
sig <- exp(simpdf$value[simpdf$parameters=="logSigObs"])

deriv_posteriors<-data.frame(chains=rep(simpdf$chains[simpdf$parameters=="logbeta"],4),
                             parameters = rep(c("a","b","S[max]","sigma"),each=length(a)),
                             value = c(a,beta,Smax,sig)
                             )

# create figure with posteriors
pm <- ggplot(deriv_posteriors) 
pm <- pm + geom_density(aes(x=value, y=..scaled.., color=chains), size=1.2, stat="density")
pm <- pm + facet_wrap(~parameters, scales="free",
  labeller = label_parsed)
pm <- pm + theme_bw(10)+labs(colour = "Chain", x="Value", y="Density")
pm <- pm + scale_color_viridis_d(end = 0.8,option = "A") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

pm
ggsave("../Figures/posterior_simple_model_fit.pdf", plot=pm, width=10,height=7)

# create table with estimates
# MLE Estimates are already compiled in All_Ests
MLE_Ests <- All_Ests %>% filter(Param %in% c("a", "alpha", "beta", "Smax", "sig")) %>%
                                          arrange( factor(Param, levels = c("a", "alpha", "beta", "Smax", "sig")))

simpleparb<-data.frame(parameters=c("a", "alpha","b","Smax", "sigma"),
                       MLE= MLE_Ests$Estimate,
                       MLE_Low = MLE_Ests$Lower,
                       MLE_Upper = MLE_Ests$Upper,
                       median=c(median(a), 
                                median(alpha),
                                median(beta), 
                                median(Smax), 
                                median(sig)),
                       "Lower 95% CI"= c(quantile(a, probs=0.025),
                             quantile(alpha, probs=0.025), 
                             quantile(beta, probs=0.025), 
                             quantile(Smax, probs=0.025), 
                             quantile(sig, probs=0.025)),
                       "Upper 95% CI" =c(quantile(a, probs=0.975),
                              quantile(alpha, probs=0.975), 
                              quantile(beta, probs=0.975), 
                              quantile(Smax, probs=0.975), 
                              quantile(sig, probs=0.975)))

write.csv(simpleparb, "../DataOut/Simple_Ricker_Ests.csv")


#=============================================================================================================
#Autocorrelated recruitment model 

mydataar<-list(obs_R=SR$R/Scale,obs_S=SR$S_adj/Scale)

parameters_autocorr <- list(
  a=(a_srm),
  logbeta = log(b_srm*Scale),
  logSigObs= log(.4), 
  rho=0.3)

# Need a bit narrower bounds for this MCMC to converge
# also bound beta so that Smax is between 1 Spawner and 2*maximum observed
Min_S_Obs <- 1/Scale
Max_S_Obs <- max(SR$S_adj)/Scale
Min_LogB <- log(1/(Max_S_Obs*2))
Max_LogB <- log(1/(Min_S_Obs))

# create vector of bounds
#           a      logbeta   logSigObs  rho
lower <- c(-Inf,   Min_LogB, -Inf,       -Inf)
upper <- c(Inf,    Max_LogB,  Inf,       Inf)

objar <- MakeADFun(mydataar,parameters_autocorr,random=NULL,DLL="Ricker_autocorr_ch")
#newtonOption(objar, smartsearch=FALSE)
  
optar <- nlminb(objar$par,objar$fn,objar$gr, lower = lower, upper=upper)

# pull out estimates from ML fit
# Create Table of outputs
All_Ests_ar <- data.frame(summary(sdreport(objar)))
All_Ests_ar$Param <- row.names(All_Ests_ar)
All_Ests_ar

# Add column with upper and lower CIs
All_Ests_ar <- All_Ests_ar %>% mutate(Lower = Estimate-1.96*Std..Error, Upper = Estimate+1.96*Std..Error)
All_Ests_ar[All_Ests_ar$Param == "beta", c("Estimate", "Std..Error", "Lower", "Upper")] <- All_Ests_ar[ All_Ests_ar$Param == "beta", c("Estimate", "Std..Error", "Lower", "Upper")]/Scale
All_Ests_ar[All_Ests_ar$Param == "Smax", c("Estimate", "Std..Error", "Lower", "Upper")] <- All_Ests_ar[ All_Ests_ar$Param == "Smax", c("Estimate", "Std..Error", "Lower", "Upper")]*Scale                      

# compile estimates and residuals
SRdiagauto<-SR
SRdiagauto$residuals <- objar$report()$residuals
SRdiagauto$std_resids <- objar$report()$std_resids
SRdiagauto$alpha <- objar$report()$alpha
SRdiagauto$model <- "Ricker + Autocorr. "

# look at residuals
rpar <- ggplot(SRdiagauto) +
     geom_line(aes(x=BroodYear,y=residuals),size=1.2) +
     geom_point(aes(x=BroodYear,y=residuals, col=residuals),stroke=3) +
     geom_hline(yintercept = 0) +
    theme_bw(10) +
     scale_color_viridis_c(end = 0.8)
rpar

# look at standardized residuals
ggplot(SRdiagauto) +
  geom_line(aes(x=BroodYear,y=std_resids),size=1.2) +
  geom_point(aes(x=BroodYear,y=std_resids, col=std_resids),stroke=3) +
  geom_hline(yintercept = 0) +
  theme_bw(10) +
  scale_color_viridis_c(end = 0.8)

# 2 residuals more than 2 sd's away


##MCMC
autoB<-list(
  obj=objar,
  nchain=3,
  iter=iteracs,
  lowbd=lower,  # should revisit these bounds
  hibd=upper
)
posterior_ar <- posteriorsdf(autoB)

#posteriors of derived quantities -- the interesting ones
ardf <- posterior_ar$posteriors
a_ar <- (ardf$value[ardf$parameters=="a"])
alpha_ar <- exp(ardf$value[ardf$parameters=="a"])
beta_ar <- exp(ardf$value[ardf$parameters=="logbeta"])/Scale
logbeta_ar <- (ardf$value[ardf$parameters=="logbeta"])
Smax_ar <- 1/exp(ardf$value[ardf$parameters=="logbeta"])*Scale
rho_ar <- 2 * boot::inv.logit(ardf$value[ardf$parameters=="rho"])-1
sig_ar <- exp(ardf$value[ardf$parameters=="logSigObs"])*sqrt(1-rho_ar^2)


# plot of posteriors
deriv_posteriorsar<-data.frame(chains=rep(ardf$chains[ardf$parameters=="logbeta"],6),
                               parameters = rep(c("a","alpha","b","S[max]","sigma[AR]","rho"), each=length(a_ar)),
                               value = c(a_ar, alpha_ar, beta_ar, Smax_ar, sig_ar, rho_ar) )

pmar <- ggplot(deriv_posteriorsar) +
       geom_density(aes(x=value, color=chains), size=1.2) +
       facet_wrap(~parameters, scales="free",labeller = label_parsed) +
       theme_bw(10)+labs( x="Value", y="Density", color = "Chain") +
       scale_color_viridis_d(end = 0.8, option = "A") +
       theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
pmar

ggsave("../Figures/posterior_ar_model_fit.pdf", plot=pmar, width=10,height=7)


# create table with estimates

# MLE Estimates are already compiled in All_Ests
# just set up so can go straight in table
MLE_Ests_ar <- All_Ests_ar %>% filter(Param %in% c("a", "alpha", "beta", "Smax", "rhoo", "SigObs")) %>%
  arrange( factor(Param, levels = c("a", "alpha", "beta", "Smax", "rhoo", "SigObs")))

arpar <- data.frame(parameters=c("a", "alpha","b","Smax","rho", "sigma"),
                    MLE= MLE_Ests_ar$Estimate,
                    MLE_Low = MLE_Ests_ar$Lower,
                    MLE_Upper = MLE_Ests_ar$Upper,
                       median=c(median(a_ar), 
                                median(alpha_ar),
                                median(beta_ar), 
                                median(Smax_ar), 
                                median(rho_ar), 
                                median(sig_ar)),
                       "Lower 95% CI"= c(quantile(a_ar, probs=0.025),
                                         quantile(alpha_ar, probs=0.025), 
                                         quantile(beta_ar, probs=0.025), 
                                         quantile(Smax_ar, probs=0.025), 
                                         quantile(rho_ar, probs=0.025), 
                                         quantile(sig_ar, probs=0.025)),
                       "Upper 95% CI" =c(quantile(alpha_ar, probs=0.975),
                                         quantile(alpha_ar, probs=0.975), 
                                         quantile(beta_ar, probs=0.975), 
                                         quantile(Smax_ar, probs=0.975),
                                         quantile(rho_ar, probs=0.975),
                                         quantile(sig_ar, probs=0.975)))

write.csv(arpar, "../DataOut/AR_Model_Ests.csv")


#============================================================================================
# Time-Varying Ricker

mydatatimevar<-list(obs_R=SR$R/Scale, obs_S=SR$S_adj/Scale, prbeta1=3, prbeta2=3)

parameters_recursive <- list(
  a0= as.numeric(a_srm),
  logbeta = as.numeric(log(b_srm*Scale)),
  rho=.5,
  logvarphi= 0.2,
  a=rep(a_srm,length(SR$R))
)

obj_timevar<-MakeADFun(mydatatimevar, parameters_recursive, random="a", DLL="Rickerkf_ratiovar")

#newtonOption(obj_timevar, smartsearch=FALSE)

# Need a bit narrower bounds again for this MCMC to converge
# also bound beta so that Smax is between min Spawners and 2*maximum observed
Min_S_Obs <- min(SR$S_adj)/Scale
Max_S_Obs <- max(SR$S_adj)/Scale
Min_LogB <- log(1/(Max_S_Obs*2))
Max_LogB <- log(1/(Min_S_Obs))

# bound rho to be between 0 and 1
#         a0      logB   rho logvar  a's
lower <- c(0,   Min_LogB, 0, -10, rep(0,length(SR$R)))
upper <- c(10, Max_LogB, 1,  10, rep(10,length(SR$R)))

opt_timevar<-nlminb(obj_timevar$par,obj_timevar$fn,obj_timevar$gr, 
                    lower=lower, upper = upper,
                    control = list(eval.max = 1e5, iter.max = 1e5))

All_Ests_TV <- data.frame(summary(sdreport(obj_timevar)))
All_Ests_TV$Param <- row.names(All_Ests_TV)
All_Ests_TV

# Add column with upper and lower CIs
All_Ests_TV <- All_Ests_TV %>% mutate(Lower = Estimate-1.96*Std..Error, Upper = Estimate+1.96*Std..Error)
All_Ests_TV[All_Ests_TV$Param == "beta", c("Estimate", "Std..Error", "Lower", "Upper")] <- All_Ests_TV[ All_Ests_TV$Param == "beta", c("Estimate", "Std..Error", "Lower", "Upper")]/Scale
All_Ests_TV[All_Ests_TV$Param == "Smax", c("Estimate", "Std..Error", "Lower", "Upper")] <- All_Ests_TV[ All_Ests_TV$Param == "Smax", c("Estimate", "Std..Error", "Lower", "Upper")]*Scale                      


repkf<-obj_timevar$report()
SRtimevar <-SR
SRtimevar$model <- "Ricker + T.V. alpha"
SRtimevar$residuals <- repkf$residuals
SRtimevar$std_resids <- repkf$std_resids
SRtimevar$alpha <- repkf$alpha

# look at standardized residuals
ggplot(SRtimevar) +
  geom_line(aes(x=BroodYear,y=std_resids),size=1.2) +
  geom_point(aes(x=BroodYear,y=std_resids, col=std_resids),stroke=3) +
  geom_hline(yintercept = 0) +
  theme_bw(10) +
  scale_color_viridis_c(end = 0.8)

# smallest residuals, nothing close to 2

# now run as MCMC
recursiveB<-list(
  obj=obj_timevar,
  nchain=3,
  iter=iteracs,
  lowbd = lower,
  hibd = upper
)

# do outside function so can get derived values
fitmcmc1 <- tmbstan(obj_timevar, chains=recursiveB$nchain,
                    iter=recursiveB$iter, init="random",
                    lower=recursiveB$lowbd, upper=recursiveB$hibd,
                    control = list(adapt_delta = 0.98))

mc <- extract(fitmcmc1, pars=names(obj_timevar$par),
              inc_warmup=TRUE, permuted=FALSE)

recrsdf <- melt(as.array(fitmcmc1))

# need to manually pull out value for a_avg since derived quantity
All_Posts <- as.matrix(fitmcmc1)

# pull out R_Fit values -- don't know why this code doesn't work
# a_avg_med <- obj_timevar$report(All_Posts[1,-ncol(All_Posts)])$a_avg
# a_avg <- matrix(NA, nrow=nrow(All_Posts), ncol = length(a_avg_med))
# for(i in 1:nrow(All_Posts)){
#   r <- obj_timevar$report(All_Posts[i,-ncol(All_Posts)])
#   a_avg[i,] <- r$a_avg
# }

hista_rb<-(recrsdf$value[recrsdf$parameters=="a0"])
beta_rb<-exp(recrsdf$value[recrsdf$parameters=="logbeta"]) / Scale
Smax_rb<-1/exp(recrsdf$value[recrsdf$parameters=="logbeta"]) * Scale
rho_rb<-(recrsdf$value[recrsdf$parameters=="rho"])
sigtot <- 1/sqrt(exp(recrsdf$value[recrsdf$parameters=="logvarphi"])) # BD added sqrt
soa <- recrsdf[grep("a",recrsdf$parameters),]

# get median, CI's for alpha[27] - alpha[30]
a_summ <- soa %>% filter(parameters %in% c("a[27]", "a[28]", "a[29]", "a[30]")) %>%
            group_by(parameters) %>% 
            summarise(median = quantile(value, probs=c(.5)), 
                      low = quantile(value, probs=c(.025)),
                      up= quantile(value, probs = 0.975))


# for each chain, iteration, get a_avg
a_avg_post <- soa %>% filter(parameters %in% c("a[27]", "a[28]", "a[29]", "a[30]")) %>%
               group_by(iterations, chains) %>% summarise(a_avg = mean(value))

# MLE Estimates are already compiled in All_Ests
# just set up so can go straight in table
MLE_Ests_TV <- All_Ests_TV %>% 
              filter(Param %in% c("a.26", "a.27", "a.28", "a.29", "a_avg", "alpha_avg", "beta", "Smax", "rho", "theta")) %>%
              arrange( factor(Param, levels = c("a.26", "a.27", "a.28", "a.29", "a_avg", "alpha_avg", "beta", "Smax", "rho", "theta")))


timevar_par<-data.frame(parameters=c("a_[2010]","a_[2011]","a_[2012]","a_[2013]","a_[AVG]",
                                     "alpha[AVG]","b","Smax","rho", "sigma[Total]"),
                        MLE=MLE_Ests_TV$Estimate,
                        MLE_Low = MLE_Ests_TV$Lower,
                        MLE_Upper = MLE_Ests_TV$Upper,
                        median=c(a_summ$median,
                                 median(a_avg_post$a_avg),
                                 exp(median(a_avg_post$a_avg)),
                                 median(beta_rb), 
                                 median(Smax_rb),
                                 median(rho_rb),
                                 median(sigtot)),
                        low=c(a_summ$low,
                              quantile(a_avg_post$a_avg, probs=0.025), 
                              exp(quantile(a_avg_post$a_avg, probs=0.025)),
                              quantile(beta_rb, probs=0.025), 
                              quantile(Smax_rb, probs=0.025), 
                              quantile(rho_rb, probs=0.025), 
                              quantile(sigtot, probs=0.025)),
                        high=c(a_summ$up,
                               quantile(a_avg_post$a_avg, probs=0.975),
                               exp(quantile(a_avg_post$a_avg, probs=0.975)),
                               quantile(beta_rb, probs=0.975), 
                               quantile(Smax_rb, probs=0.975), 
                               quantile(rho_rb, probs=0.975), 
                               quantile(sigtot, probs=0.975)))

write.csv(timevar_par, "../DataOut/TimeVar_Ests.csv")

# collect posteriors to plot
deriv_posteriorsrb <- data.frame(chains=rep(recrsdf$chains[recrsdf$parameters=="logbeta"],4),
                               parameters = rep(c("a[avg]","b","sigma[tot]","rho"),each=length(hista_rb)),
                               value = c(a_avg_post$a_avg, beta_rb, sigtot, rho_rb) )

pmrb <- ggplot(deriv_posteriorsrb) +
       geom_density(aes(x=value, color=chains), size=1.2) +
       facet_wrap(~parameters, scales="free",labeller = label_parsed) +
       theme_bw(10)+labs(colour = "Prior", x="Value", y="Density") +
       scale_color_viridis_d(end = 0.8,option = "A")
pmrb

ggsave("../Figures/TV_Model_Posts.pdf", plot=pmrb, width=10,height=7)

#=================================================================================
# comparisons

allfits <- bind_rows(SRdiagsimple, SRdiagauto, SRtimevar)
#------------------
# residuals by BY
pp <- ggplot(allfits )  +
      geom_line(aes(x=BroodYear,y=residuals),size=1.2)  +
      geom_point(aes(x=BroodYear,y=residuals, col=residuals),stroke=3)  +
      geom_hline(yintercept = 0)  +
      theme_bw(16)  +
      scale_color_viridis_c(end = 0.8)  +
      facet_wrap(~model)  
pp 

ggsave("../Figures/Residuals_by_BroodYear.pdf", plot=pp, width=11,height=7)
#---------------------
# std. residuals by BY
pp <- ggplot(allfits )  +
  geom_line(aes(x=BroodYear,y=std_resids),size=1.2)  +
  geom_point(aes(x=BroodYear,y=std_resids, col=std_resids),stroke=3)  +
  geom_hline(yintercept = 0)  +
  theme_bw(16)  +
  scale_color_viridis_c(end = 0.8)  +
  facet_wrap(~model)  
pp 

ggsave("../Figures/Std_Residuals_by_BroodYear.pdf", plot=pp, width=11,height=7)

#-----------------------------------------------------------
#Residuals vs Spawners comparison over years with all models
psp <- ggplot(allfits )  +
        geom_point(aes(x=S_adj,y=residuals, col=residuals),stroke=3)  +
        geom_hline(yintercept = 0)  +
        theme_bw(16) + xlab("Spawners")  +
        scale_color_viridis_c(end = 0.8)  +
        facet_wrap(~model)  
psp 

ggsave("../Figures/Residuals_by_Spawners.pdf", plot=psp, width=11,height=7)

#-----------------------------------------------------------
# Std Residuals vs Spawners comparison over years with all models
psp <- ggplot(allfits )  +
  geom_point(aes(x=S_adj,y=std_resids, col=std_resids),stroke=3)  +
  geom_hline(yintercept = 0)  +
  theme_bw(16) + xlab("Spawners")  +
  scale_color_viridis_c(end = 0.8)  +
  facet_wrap(~model)  
psp 

ggsave("../Figures/Residuals_by_Spawners.pdf", plot=psp, width=11,height=7)

#------------------------------------------
# Productivity comparison with all models.
ppa <- ggplot(allfits )  +
      geom_line(aes(x=BroodYear,y=alpha, col=model),size=1.2)  +
      geom_point(aes(x=BroodYear,y=alpha, col=model))  +
      theme_bw(16)  +
      scale_color_viridis_d(end = 0.9,option = "C")  
ppa 
ggsave("../Figures/Productivity_Comparison.pdf", plot=ppa, width=10,height=7)


# AIC comparison
AICs <- c(TMBAIC(opt), TMBAIC(optar),  TMBAIC(opt_timevar))
AIC_Tab <- data.frame(model=c("Simple ricker", "Ricker + AutoCorr", "Ricker + T.V. Alpha"),
                 AIC= AICs,
                 deltaAIC= AICs - min(AICs))

write.csv(AIC_Tab, "../DataOut/AICTab.csv")
