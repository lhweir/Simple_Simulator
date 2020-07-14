

setwd("C:/Users/WeirL/Documents/GitHub/Simple_Simulator/DU2/DataOut")

library(ggplot2)
library(hrbrthemes)
library(RColorBrewer)
library(viridis)

results<-read.csv("Recovery_Results_DU2.csv")

x.labels <- c("-50","-40","-30","-20","-10","0","+10","+20","+30")
y.labels <- c("-100","-90","-80","-70","-60","-50","-40","-30","-20","-10","0","+10")


F.scenarios<-as.character(unique(results$FishScenario))
P.scenarios<-as.character(unique(results$ProdScenario))

plot.data<-expand.grid(X=P.scenarios,Y=F.scenarios)


  for (i in 1:length(plot.data$X)) {
 
  plot.data$Count.UpRT[i] <- sum(results$Upp_RT_met[results$ProdScenario == plot.data$X[i] & results$FishScenario == plot.data$Y[i]],na.rm = TRUE)/max(results$Sim)
  plot.data$Count.LowRT[i] <- sum(results$Low_RT_met[results$ProdScenario == plot.data$X[i] & results$FishScenario == plot.data$Y[i]],na.rm = TRUE)/max(results$Sim)

  if (plot.data$Y[i] == "Fish0") {
    plot.data$fishscenario[i] <- 0 }
  if (plot.data$Y[i] == "Fish+10") {
    plot.data$fishscenario[i] <- 10 }
  if (plot.data$Y[i] == "Fish10") {
    plot.data$fishscenario[i] <- -10 }
  if (plot.data$Y[i] == "Fish20") {
    plot.data$fishscenario[i] <- -20 }
  if (plot.data$Y[i] == "Fish30") {
    plot.data$fishscenario[i] <- -30 }
  if (plot.data$Y[i] == "Fish40") {
    plot.data$fishscenario[i] <- -40 }
  if (plot.data$Y[i] == "Fish50") {
    plot.data$fishscenario[i] <- -50 }
  if (plot.data$Y[i] == "Fish60") {
    plot.data$fishscenario[i] <- -60 }
  if (plot.data$Y[i] == "Fish70") {
    plot.data$fishscenario[i] <- -70 }
  if (plot.data$Y[i] == "Fish80") {
    plot.data$fishscenario[i] <- -80 }
  if (plot.data$Y[i] == "Fish90") {
    plot.data$fishscenario[i] <- -90 }
  if (plot.data$Y[i] == "Fish100") {
    plot.data$fishscenario[i] <- -100 }
    
  
  if (plot.data$X[i] == "0p_12yrs") {
    plot.data$prodscenario[i] <- 0 }
  if (plot.data$X[i] == "10p_12yrs") {
    plot.data$prodscenario[i] <- 10 }
  if (plot.data$X[i] == "20p_12yrs") {
    plot.data$prodscenario[i] <- 20 }
  if (plot.data$X[i] == "30p_12yrs") {
    plot.data$prodscenario[i] <- 30 }
  if (plot.data$X[i] == "-10p_12yrs") {
    plot.data$prodscenario[i] <- -10 }
  if (plot.data$X[i] == "-20p_12yrs") {
    plot.data$prodscenario[i] <- -20 }
  if (plot.data$X[i] == "-30p_12yrs") {
    plot.data$prodscenario[i] <- -30 }
  if (plot.data$X[i] == "-40p_12yrs") {
    plot.data$prodscenario[i] <- -40 }
  if (plot.data$X[i] == "-50p_12yrs") {
    plot.data$prodscenario[i] <- -50 }
 
  
  }


  for (i in 1:length(plot.data$X)) {
    
    if (plot.data$Count.UpRT[i] >= 0 & plot.data$Count.UpRT[i] <= 0.01) {
      plot.data$BinsUp[i] <- "Extr. unlikely" }
    if (plot.data$Count.UpRT[i] > 0.01 & plot.data$Count.UpRT[i] <= 0.1) {
      plot.data$BinsUp[i] <- "Very unlikely" }
    if (plot.data$Count.UpRT[i] > 0.1 & plot.data$Count.UpRT[i] <= 0.33) {
      plot.data$BinsUp[i] <- "Unlikely" }
    if (plot.data$Count.UpRT[i] > 0.33 & plot.data$Count.UpRT[i] <= 0.60) {
      plot.data$BinsUp[i] <- "As likely as not" }
    if (plot.data$Count.UpRT[i] > 0.60 & plot.data$Count.UpRT[i] <= 0.90) {
      plot.data$BinsUp[i] <- "Likely" }
    if (plot.data$Count.UpRT[i] > 0.90 & plot.data$Count.UpRT[i] <= 0.99) {
      plot.data$BinsUp[i] <- "Very likely" }
    if (plot.data$Count.UpRT[i] > 0.99 & plot.data$Count.UpRT[i] <= 1) {
      plot.data$BinsUp[i] <- "Certain" }
    
    if (plot.data$Count.LowRT[i] >= 0 & plot.data$Count.LowRT[i] <= 0.01) {
      plot.data$BinsLow[i] <- "Extr. unlikely" }
    if (plot.data$Count.LowRT[i] > 0.01 & plot.data$Count.LowRT[i] <= 0.1) {
      plot.data$BinsLow[i] <- "Very unlikely" }
    if (plot.data$Count.LowRT[i] > 0.1 & plot.data$Count.LowRT[i] <= 0.33) {
      plot.data$BinsLow[i] <- "Unlikely" }
    if (plot.data$Count.LowRT[i] > 0.33 & plot.data$Count.LowRT[i] <= 0.60) {
      plot.data$BinsLow[i] <- "As likely as not" }
    if (plot.data$Count.LowRT[i] > 0.60 & plot.data$Count.LowRT[i] <= 0.90) {
      plot.data$BinsLow[i] <- "Likely" }
    if (plot.data$Count.LowRT[i] > 0.90 & plot.data$Count.LowRT[i] <= 0.99) {
      plot.data$BinsLow[i] <- "Very likely" }
    if (plot.data$Count.LowRT[i] > 0.99 & plot.data$Count.LowRT[i] <= 1) {
      plot.data$BinsLow[i] <- "Certain" }
    
  }

# IPPCC categories used in Sockeye

for (i in 1:length(plot.data$X)) {
  
  if (plot.data$Count.UpRT[i] >= 0 & plot.data$Count.UpRT[i] <= 0.01) {
    plot.data$BinsUp[i] <- "Extr. unlikely" }
  if (plot.data$Count.UpRT[i] > 0.01 & plot.data$Count.UpRT[i] <= 0.1) {
    plot.data$BinsUp[i] <- "Very unlikely" }
  if (plot.data$Count.UpRT[i] > 0.1 & plot.data$Count.UpRT[i] <= 0.33) {
    plot.data$BinsUp[i] <- "Unlikely" }
  if (plot.data$Count.UpRT[i] > 0.33 & plot.data$Count.UpRT[i] <= 0.66) {
    plot.data$BinsUp[i] <- "As likely as not" }
  if (plot.data$Count.UpRT[i] > 0.66 & plot.data$Count.UpRT[i] <= 0.90) {
    plot.data$BinsUp[i] <- "Likely" }
  if (plot.data$Count.UpRT[i] > 0.90 & plot.data$Count.UpRT[i] <= 0.99) {
    plot.data$BinsUp[i] <- "Very likely" }
  if (plot.data$Count.UpRT[i] > 0.99 & plot.data$Count.UpRT[i] <= 1) {
    plot.data$BinsUp[i] <- "Certain" }
  
  if (plot.data$Count.LowRT[i] >= 0 & plot.data$Count.LowRT[i] <= 0.01) {
    plot.data$BinsLow[i] <- "Extr. unlikely" }
  if (plot.data$Count.LowRT[i] > 0.01 & plot.data$Count.LowRT[i] <= 0.1) {
    plot.data$BinsLow[i] <- "Very unlikely" }
  if (plot.data$Count.LowRT[i] > 0.1 & plot.data$Count.LowRT[i] <= 0.33) {
    plot.data$BinsLow[i] <- "Unlikely" }
  if (plot.data$Count.LowRT[i] > 0.33 & plot.data$Count.LowRT[i] <= 0.66) {
    plot.data$BinsLow[i] <- "As likely as not" }
  if (plot.data$Count.LowRT[i] > 0.66 & plot.data$Count.LowRT[i] <= 0.90) {
    plot.data$BinsLow[i] <- "Likely" }
  if (plot.data$Count.LowRT[i] > 0.90 & plot.data$Count.LowRT[i] <= 0.99) {
    plot.data$BinsLow[i] <- "Very likely" }
  if (plot.data$Count.LowRT[i] > 0.99 & plot.data$Count.LowRT[i] <= 1) {
    plot.data$BinsLow[i] <- "Certain" }
  
}



write.csv(plot.data, "PlotdataBCF.csv")

png("Heat Map - Upper Recovery Target-b&w.png", width=11, height = 8, units="in", res=800)
ggplot(plot.data, aes(prodscenario, fishscenario, fill=BinsUp)) + 
  geom_tile() +
  geom_point(aes(x=0, y=0), colour="red", size=4, show.legend=FALSE)+
  scale_fill_manual(breaks=c("Extr. unlikely","Very unlikely", "Unlikely","As likely as not","Likely"),
                    labels=c("Extr. unlikely (<1%)","Very unlikely (1-10%)", "Unlikely (10-33%)","As likely as not (33-66%)","Likely (66-90%)"),
                    values = c("black","grey15","grey43","grey70","grey96" ), name="Ability to Reach\nUpper Recovery Target")+
  scale_x_continuous(name=c("Productivity Change"),breaks=c(-50,-40,-30,-20,-10,0,10,20,30),labels=c("-50", "-40", "-30","-20","-10","0","10","20","30"))+
  scale_y_continuous(name=c("% Change in Canadian Harvest Rates"),breaks=c(10,0,-10,-20,-30,-40,-50,-60,-70,-80,-90,-100),labels=c("10", "0", "-10","-20","-30","-40","-50","-60","-70","-80","-90","-100"))+
  theme_minimal()+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=18), 
        axis.title = element_text(size=18),legend.text=element_text(size=18), legend.title=element_text(size=18))
dev.off()

# viridis(5)
# c("grey15","grey43","grey70","grey96")
#rCB - c( "#7b3294","#c2a5cf", "#f7f7f7", "#a6dba0","#008837")
#c( "#ef8a62","#fddbc7", "#f7f7f7", "#d1e5f0","#67a9cf")

png("Heat Map - Lower Recovery Target-b&w.png", width=11, height = 8, units="in", res=800)
ggplot(plot.data, aes(prodscenario, fishscenario, fill=BinsLow)) + 
  geom_tile() +
  xlab("Productivity Change") + ylab("% Change in Canadian Fishing Related Mortality")+
  scale_fill_manual(breaks=c("Extr. unlikely","Very unlikely", "Unlikely","As likely as not","Likely"),
                    labels=c("Extr. unlikely (<1%)","Very unlikely (1-10%)", "Unlikely (10-33%)","As likely as not (33-66%)","Likely (66-90%)"),
                    values = c("black","grey15","grey43","grey70","grey96"), name="Ability to Reach\nLower Recovery Target")+
  scale_x_continuous(name=c("Productivity Change"),breaks=c(-50,-40,-30,-20,-10,0,10,20,30),labels=c("-50", "-40", "-30","-20","-10","0","10","20","30"))+
  scale_y_continuous(name=c("% Change in Canadian Harvest Rates"),breaks=c(10,0,-10,-20,-30,-40,-50,-60,-70,-80,-90,-100),labels=c("10", "0", "-10","-20","-30","-40","-50","-60","-70","-80","-90","-100"))+
  theme_minimal()+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=18), 
        axis.title = element_text(size=18),legend.text=element_text(size=18), legend.title=element_text(size=18))+
  geom_point(aes(x=0, y=0), colour="red", size=4, show.legend=FALSE)
dev.off()