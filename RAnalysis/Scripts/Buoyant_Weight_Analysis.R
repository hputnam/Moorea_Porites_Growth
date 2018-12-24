#Title: Porites Buoyant Weight Data
#Author: HM Putnam
#Date Last Modified: 20181220
#See Readme file for details

rm(list=ls()) #clears workspace 

#Install Libraries
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("ggplot2" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggplot2') 
if ("sciplot" %in% rownames(installed.packages()) == 'FALSE') install.packages('sciplot') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("reshape2" %in% rownames(installed.packages()) == 'FALSE') install.packages('reshape2') 
if ("tidyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggridges') 
if ("ggridges"  %in% rownames(installed.packages()) == 'FALSE') install.packages('ggridges') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 


#Read in required libraries
##### Include Versions of libraries
library("plyr")
library("ggplot2")
library("sciplot")
library("plotrix")
library("reshape2")
library("nlme") #mixed model, repeated measures ANOVA
library("lsmeans") #mixed model posthoc  statistical comparisons
library("tidyr")
library("ggridges")
library("gridExtra")


#Analysis of juvenile massive Porites growth through time in 2018

#Timepoints
#Time0 - 20180127
#Time1 - 20180317 
#Time2 - 20180507
#Time3 - 20180827
#Time4 - 20181101

#rate1 = 20180317 - 20180127 Jan - March
#rate2 = 20180507 - 20180317 March - May
#rate3 = 20180827 - 20180507 May - August
#rate4 = 20181101 - 20180827 August - November

# Set Working Directory:
setwd("/Users/hputnam/MyProjects/Moorea_Porites_Growth/RAnalysis/Data/") #set working

# load data 
CalData<-read.csv('Bouyant_weight_calibration_curve.csv', header=T, sep=",")
#Data with Standard weights in grams in deionized water, saltwater, and the temp that they were measured at

#plot relationship between temp and Standard in fresh and salt water
plot(CalData$Temp, CalData$StandardSalt, col = 'red', ylab = 'Salt Weight (g)', xlab = 'Temp C')
par(new=TRUE)
plot(CalData$Temp, CalData$StandardFresh, col = 'blue',xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("'Fresh Weight (g)'",side=4,line=3)
legend("topleft",col=c("red","blue"),lty=1,legend=c("Salt","Fresh"), bty = "n")

#linear regression between temp and Standard
#Salt water standard measures
StandardSaltModel <- lm(CalData$StandardSalt~CalData$Temp)
summary(StandardSaltModel)
summary(StandardSaltModel)$r.squared
StandardSaltModel$coef

#Fresh water standard measures
StandardFreshModel <- lm(CalData$StandardFresh~CalData$Temp)
summary(StandardFreshModel)
summary(StandardFreshModel)$r.squared
StandardFreshModel$coef

#load coral weight data
BW.data <-read.csv('Plob_Buoyant_Weight.csv', header=T, na.strings = "NA") 
BW <- as.numeric(as.character(BW.data$Mass.g)) #assign mass data to BW for function
Temp <- BW.data$Temp.C #assign temp to temp for function

BWCalc <- function(StandardAir= 39.108, Temp, BW, CoralDensity = 2.93){
  #Parameters---------------------------------------------------------------
  # StandardAir is the weight of the Standard in air (Default set at 39.108 grams weighed on balance)
  # Temp is the temperature of the water during the coral measurement
  # BW is the buoywant weight of the coral
  # CoralDensity <- 2.03 #g cm-3, set aragonite density for Montipora from literature 
  #Montipora Arag = 2.03 g cm-3 average from table 2 in Anthony 2003 Functional Ecology 17:246-259
  # You can change the density to literatrue values for the specific coral species of interest in the function. 
  
  #Calculation------------------------------------------------------------
  #Step 1: correct the Standard weights for temperature
  # StandardFresh is the weight of the Standard in  fresh water at temp x
  # StandardSalt is the weight of the Standard in salt water at temp x
  
  # This is based on a calibration curve for Standards weighed in fresh and salt water at many temps
StandardFresh <- StandardFreshModel$coef[-1]*Temp + StandardFreshModel$coef[1] 
StandardSalt <- StandardSaltModel$coef[-1]*Temp + StandardSaltModel$coef[1] 
  
  # Step 2: Use weight in air and freshwater of the glass Standard to calculate
  # the density of the Standard (Spencer Davies eq. 4)
FreshDensity <- 1 #Fresh water has a density of 1 g/cm3
StandardDensity <- FreshDensity/(1-(StandardFresh/StandardAir)) 
  
  # Step 3: Calculate the density of seawater using the density of the Standard
  # (Spencer Davies eq. 3)
SWDensity <- StandardDensity*(1-(StandardSalt/StandardAir))
  
  # Step 4: Calculate the dry weight of the coral (Spencer Davies eq. 1)
CoralWeight <- BW/(1-(SWDensity/CoralDensity))
  
  return(CoralWeight) #returns coral dry weights in g
}
  
BW.data$Dry.Weigh.g <- BWCalc(BW=BW, Temp=Temp) #use function to calculate dry weight
hist(BW.data$Dry.Weigh.g) #examine data
boxplot(BW.data$Dry.Weigh.g ~BW.data$TimePoint) #examine data


data.0.curve <- subset(BW.data, Date==20180127)
data.1.curve <- subset(BW.data, Date==20180317)
data.2.curve <- subset(BW.data, Date==20180507)
data.3.curve <- subset(BW.data, Date==20180827)
data.4.curve <- subset(BW.data, Date==20181101)

par(mfrow=c(2,3))
boxplot(data.0.curve$Dry.Weigh.g ~data.0.curve$TimePoint, main="Jan", ylim=c(0,50)) #examine data
boxplot(data.1.curve$Dry.Weigh.g ~data.1.curve$TimePoint, main="March", ylim=c(0,50)) #examine data
boxplot(data.2.curve$Dry.Weigh.g ~data.2.curve$TimePoint, main="May", ylim=c(0,50)) #examine data
boxplot(data.3.curve$Dry.Weigh.g ~data.3.curve$TimePoint, main="August", ylim=c(0,50)) #examine data
boxplot(data.4.curve$Dry.Weigh.g ~data.4.curve$TimePoint, main="November", ylim=c(0,50)) #examine data

#Compare with single point standard
#load coral weight data
BWdata <-read.csv('Plob_Buoyant_Weight.csv', header=T, na.strings = "NA") 
BW.std.data <- read.csv("BW_standards.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
Arag <- 2.93 # density of aragonite g cm-3
BW.std.data$FW.Dens <- ((-0.000005*BW.std.data$Temp*BW.std.data$Temp)+(0.000007*BW.std.data$Temp)+1.0001) #g cm-3
BW.std.data$Ref.Dens <- ((BW.std.data$Ref.Dry.g*BW.std.data$FW.Dens)/(BW.std.data$Ref.Dry.g-BW.std.data$Ref.FW.g)) #g cm-3
BW.std.data$SW.Dens <- ((BW.std.data$Ref.Dens*(BW.std.data$Ref.Dry.g-BW.std.data$Ref.SW.g))/(BW.std.data$Ref.Dry.g))  #g cm-3

data.0 <- subset(BWdata, Date==20180127)
data.1 <- subset(BWdata, Date==20180317)
data.2 <- subset(BWdata, Date==20180507)
data.3 <- subset(BWdata, Date==20180827)
data.4 <- subset(BWdata, Date==20181101)
data.0$dry.weight <- ((as.numeric(as.character(data.0$Mass.g)))/(1-(BW.std.data[1,9]/Arag)))
data.1$dry.weight <- ((as.numeric(as.character(data.1$Mass.g)))/(1-(BW.std.data[2,9]/Arag)))
data.2$dry.weight <- ((as.numeric(as.character(data.2$Mass.g)))/(1-(BW.std.data[3,9]/Arag)))
data.3$dry.weight <- ((as.numeric(as.character(data.3$Mass.g)))/(1-(BW.std.data[4,9]/Arag)))
data.4$dry.weight <- ((as.numeric(as.character(data.4$Mass.g)))/(1-(BW.std.data[5,9]/Arag)))

par(mfrow=c(2,3))
plot(data.0.curve$Dry.Weigh.g,data.0$dry.weight)
plot(data.1.curve$Dry.Weigh.g,data.1$dry.weight)
plot(data.2.curve$Dry.Weigh.g,data.2$dry.weight)
plot(data.3.curve$Dry.Weigh.g,data.3$dry.weight)
plot(data.4.curve$Dry.Weigh.g,data.4$dry.weight)

#can use standard curve data for all time points

data <- reshape(BW.data, timevar = "TimePoint", drop = c("QC", "Mass.g","Temp.C"), idvar=c("Fragment.ID","Rack", "Depth" ), direction="wide")

#average length and width for diameter, half, and apply to equation for area of a hemisphere 2*pi*r^2

data$SA.0 <- (2* pi * (((as.numeric(as.character(data$Length.mm.Time0))+ as.numeric(as.character(data$Width.mm.Time0))))/4)^2)/100
data$SA.1 <- (2* pi * (((as.numeric(as.character(data$Length.mm.Time1))+ as.numeric(as.character(data$Width.mm.Time1))))/4)^2)/100
data$SA.2 <- (2* pi * (((as.numeric(as.character(data$Length.mm.Time2))+ as.numeric(as.character(data$Width.mm.Time2))))/4)^2)/100
data$SA.3 <- (2* pi * (((as.numeric(as.character(data$Length.mm.Time3))+ as.numeric(as.character(data$Width.mm.Time3))))/4)^2)/100
data$SA.4 <- (2* pi * (((as.numeric(as.character(data$Length.mm.Time4))+ as.numeric(as.character(data$Width.mm.Time4))))/4)^2)/100


#Growth Normalized to surface area of each timepoint mg cm-2 day-1
data$time1.growth <- (((data$Dry.Weigh.g.Time1 - data$Dry.Weigh.g.Time0)/(data$SA.1))/data$Days.Time1)*1000 #calculate growth rate 
data$time2.growth <- (((data$Dry.Weigh.g.Time2 - data$Dry.Weigh.g.Time1)/(data$SA.2))/data$Days.Time2)*1000 #calculate growth rate 
data$time3.growth <- (((data$Dry.Weigh.g.Time3 - data$Dry.Weigh.g.Time2)/(data$SA.3))/data$Days.Time3)*1000 #calculate growth rate 
data$time4.growth <- (((data$Dry.Weigh.g.Time4 - data$Dry.Weigh.g.Time3)/(data$SA.4))/data$Days.Time4)*1000 #calculate growth rate 


#removing -growth outliers
range(na.omit(data$time1.growth))
outs <- which(data$time1.growth==min((na.omit(data$time1.growth))))
data <- data[-outs,]
range(na.omit(data$time1.growth))
range(na.omit(data$time2.growth))
range(na.omit(data$time3.growth))
range(na.omit(data$time4.growth))
outs <- which(data$time4.growth==min((na.omit(data$time4.growth))))
data <- data[-outs,]
range(na.omit(data$time4.growth))
outs <- which(data$time4.growth==min((na.omit(data$time4.growth))))
data <- data[-outs,]
range(na.omit(data$time4.growth))


#Plotting Data by time
par(mfrow=c(2,2))
boxplot(data$time1.growth ~data$Date.Time1, ylim=c(0,5), main="Jan-March") 
boxplot(data$time2.growth ~data$Date.Time2, ylim=c(0,5), main="March-May") 
boxplot(data$time3.growth ~data$Date.Time3, ylim=c(0,5), main="May-August") 
boxplot(data$time4.growth ~data$Date.Time4, ylim=c(0,5), main="August-November") 

#Display Data by Time
data2 <-data[,c(1,2,3,49,50,51,52)]
data2 <- data2 %>% gather(Time, growth, time1.growth,time2.growth,time3.growth, time4.growth)
data2$Time <- gsub('time1.growth', 'Jan-March', data2$Time)
data2$Time <- gsub('time2.growth', 'March-May', data2$Time)
data2$Time <- gsub('time3.growth', 'May-Aug', data2$Time)
data2$Time <- gsub('time4.growth', 'Aug-Nov', data2$Time)
data2 <- na.omit(data2)
data2$Time <- factor(data2$Time, levels = c("Aug-Nov", "May-Aug", "March-May", "Jan-March"))


Gro.Joy <- ggplot(data2, aes(y=as.factor(Time),x=growth)) +
  geom_density_ridges(alpha=0.5) +
  scale_y_discrete(expand = c(0.01, 0)) +  
  scale_x_continuous(expand = c(0, 0))+
  theme(axis.text=element_text(size=20)) +
  xlab(" Growth mg cm-2 d-1") + 
  ylab("")
Gro.Joy

Gro.Joy.Figs <- arrangeGrob(Gro.Joy, ncol=1)
ggsave(file="../Output/Fig2_Growth_Hist.pdf", Gro.Joy, width = 4, height = 4, units = c("in"))

G.mean <- aggregate(growth ~ Time, data=data2, mean, na.rm=T)
G.se <- aggregate(growth ~ Time, data=data2, std.error, na.rm=T)
G.n <- aggregate(growth ~ Time, data=data2, length)

G <- cbind(G.mean, G.se$growth, G.n$growth)
colnames(G) <- c("Time", "mean", "se", "n")
G$Time <- factor(G$Time, levels = c("Jan-March", "March-May",   "May-Aug", "Aug-Nov"))

  
FigX <- ggplot(G, aes(x=Time, y=mean)) +
  geom_errorbar(aes(ymin=G$mean-G$se, ymax=G$mean+G$se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(), position = position_dodge(width = 0.2), size=1) + #plot points
  #scale_shape_manual(values=c(1,19)) + #sets point shape manually
  scale_x_discrete(limits=c("Jan-March","March-May", "May-Aug", "Aug-Nov")) +
  #geom_line(aes(linetype=Depth), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Time") + #Label the X Axis
  ylab("Growth mg cm-2 d-1") + #Label the Y Axis
  ylim(0, 3) + #set Y limits
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
FigX

data2$Time <- factor(data2$Time, levels = c("Jan-March", "March-May",   "May-Aug", "Aug-Nov"))
data2$Depth <- factor(data2$Depth, levels = c("Shallow", "Deep"))

Rack.Figs <- ggplot(data2, aes(x=Time, y=growth, fill=Depth)) +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values = c("gray", "white")) +
  geom_point(aes(shape=Depth), position = position_dodge(width = 0.9), size=1) +
  ylab("Growth mg cm-2 d-1") +
  theme_classic()+
  theme(legend.position=c(0.55,0.9), axis.text.x = element_text(angle=90))
Rack.Figs <- arrangeGrob(Rack.Figs, ncol=1)
ggsave(file="../Output/Fig1_Growth_by_Rack.pdf", Rack.Figs, width = 3, height = 4, units = c("in"))



#Test for differences in growth by rack
anova(aov(data$time1.growth ~data$Rack))
anova(aov(data$time2.growth ~data$Rack))
anova(aov(data$time3.growth ~data$Rack))
anova(aov(data$time4.growth ~data$Rack))


#reduced dataset of mass normalized growth
G.mass <- data[,c(1:3,49,50,51,52)] 

G.mass.1 <- data[,c(1:3,49)] #subset only T1 data
G.mass.1 <- (na.omit(G.mass.1)) 

G.mass.2 <- data[,c(1:3,50)] #subset only T2 data
G.mass.2 <- (na.omit(G.mass.2))

G.mass.3 <- data[,c(1:3,51)] #subset only T3 data
G.mass.3 <- (na.omit(G.mass.3))

G.mass.4 <- data[,c(1:3,52)] #subset only T3 data
G.mass.4 <- (na.omit(G.mass.4))

pdf("../Output/Fig3_growth_correlations.pdf", width=8.5, height=3.5)
par(mfrow=c(1,3))
plot(G.mass$time2.growth ~G.mass$time1.growth, 
     xlab = "Growth (Jan-March)", ylab = "Growth (March-May)", 
     xlim=c(0,4), ylim=c(0,3), pch = 16)
mod1.2 <- lm(G.mass$time2.growth ~G.mass$time1.growth)
clip(1,4,0,2.5)
abline(mod1.2)
mtext(paste0("Slope = ", round(mod1.2$coefficients[2], digits = 3),"    ", "R^2 = ", round(summary(mod1.2)$r.squared, digits = 2))) 

plot(G.mass$time3.growth ~G.mass$time2.growth, 
     xlab = "Growth (March-May)", ylab = "Growth (May-Aug)",
     xlim=c(0,4), ylim=c(0,3), pch = 16)
mod2.3 <- lm(G.mass$time3.growth ~G.mass$time2.growth)
clip(0,2.5,0,2.5)
abline(mod2.3)
mtext(paste0("Slope = ", round(mod2.3$coefficients[2], digits = 3),"    ", "R^2 = ", round(summary(mod2.3)$r.squared, digits = 2))) 

plot(G.mass$time4.growth ~G.mass$time3.growth, 
     xlab = "Growth (May-Aug)", ylab = "Growth (Aug-Nov)",
     xlim=c(0,4), ylim=c(0,3), pch = 16)
mod3.4 <- lm(G.mass$time4.growth ~G.mass$time3.growth)
clip(0,2.5,0,2.5)
abline(mod3.4)
mtext(paste0("Slope = ", round(mod3.4$coefficients[2], digits = 3),"    ", "R^2 = ", round(summary(mod3.4)$r.squared, digits = 2))) 
dev.off()










# #split each time point into 3rds
# #order by value large to small
# 
# G.mass.sorted <- G.mass[order(-G.mass$time1.growth),]
# G.mass.1.sorted <- G.mass.1[order(-G.mass.1$time1.growth),]
# G.mass.2.sorted <- G.mass.2[order(-G.mass.2$time2.growth),] 
# G.mass.3.sorted <- G.mass.3[order(-G.mass.3$time3.growth),] 
# G.mass.4.sorted <- G.mass.4[order(-G.mass.4$time4.growth),] 
# 
# ##### Time 1 as a function of Time 1 #####
# #empty
# ##### Time 2 as a function of Time 1 #####
# G.mass.high <- G.mass.sorted[1:c(as.numeric(round((nrow(G.mass.1.sorted)/3)))), ]
# G.mass.med <- G.mass.sorted[(c(as.numeric(round((nrow(G.mass.1.sorted)/3))))+1):(c(as.numeric(round((nrow(G.mass.1.sorted)))))-c(as.numeric(round((nrow(G.mass.1.sorted)/3))))), ]
# G.mass.low <- G.mass.sorted[((c(as.numeric(round((nrow(G.mass.1.sorted)))))-c(as.numeric(round((nrow(G.mass.1.sorted)/3)))))+1):c(as.numeric(round((nrow(G.mass.1.sorted))))), ]
# 
# pdf("../Output/FigureX_growth.pdf", width=8.5, height=11)
# par(mfrow=c(3,3), mar=c(4, 5, 3, 1))
# plot(G.mass.high$time1.growth, G.mass.high$time2.growth, col="red", xlim=c(0,4), ylim=c(0,2.5), ylab="Growth mg cm-2 d-1", xlab=" Growth mg cm-2 d-1", main="JANUARY-MARCH")
# points(G.mass.med$time1.growth, G.mass.med$time2.growth, col="blue")
# points(G.mass.low$time1.growth, G.mass.low$time2.growth, col="black")
# mtext(side=2, line=4, "MARCH-MAY", col="black", font=2, cex=0.9)
# #mtext(side=3, line=4, "JANUARY-MARCH", col="black", font=2, cex=0.9)
# abline(v=G.mass.med[1,8], col="red", lwd=1, lty=2)
# abline(v=G.mass.low[1,8], col="blue", lwd=1, lty=2)
# abline(lm(G.mass.high$time2.growth ~ G.mass.high$time1.growth), col="red")
# abline(lm(G.mass.med$time2.growth ~ G.mass.med$time1.growth), col="blue")
# abline(lm(G.mass.low$time2.growth ~ G.mass.low$time1.growth), col="black")
# 
# plot(0,type='n',axes=FALSE,ann=FALSE)
# plot(0,type='n',axes=FALSE,ann=FALSE)
# 
# ##### Time 3 as a function of Time 1 #####
# G.mass.sorted <- G.mass[order(-G.mass$time1.growth),]
# G.mass.high <- G.mass.sorted[1:c(as.numeric(round((nrow(G.mass.1.sorted)/3)))), ]
# G.mass.med <- G.mass.sorted[(c(as.numeric(round((nrow(G.mass.1.sorted)/3))))+1):(c(as.numeric(round((nrow(G.mass.1.sorted)))))-c(as.numeric(round((nrow(G.mass.1.sorted)/3))))), ]
# G.mass.low <- G.mass.sorted[((c(as.numeric(round((nrow(G.mass.1.sorted)))))-c(as.numeric(round((nrow(G.mass.1.sorted)/3)))))+1):c(as.numeric(round((nrow(G.mass.1.sorted))))), ]
# 
# 
# plot(G.mass.high$time1.growth, G.mass.high$time3.growth, col="red", xlim=c(0,4), ylim=c(0,2.5), ylab="Growth mg cm-2 d-1", xlab=" Growth mg cm-2 d-1")
# points(G.mass.med$time1.growth, G.mass.med$time3.growth, col="blue")
# points(G.mass.low$time1.growth, G.mass.low$time3.growth, col="black")
# mtext(side=2, line=4, "MAY-AUGUST", col="black", font=2, cex=0.9)
# abline(v=G.mass.med[1,8], col="red", lwd=1, lty=2)
# abline(v=G.mass.low[1,8], col="blue", lwd=1, lty=2)
# abline(lm(G.mass.high$time3.growth ~ G.mass.high$time1.growth), col="red")
# abline(lm(G.mass.med$time3.growth ~ G.mass.med$time1.growth), col="blue")
# abline(lm(G.mass.low$time3.growth ~ G.mass.low$time1.growth), col="black")
# 
# 
# ##### Time 3 as a function of Time 2 #####
# G.mass.sorted <- G.mass[order(-G.mass$time2.growth),]
# G.mass.high <- G.mass.sorted[1:c(as.numeric(round((nrow(G.mass.2.sorted)/3)))), ]
# G.mass.med <- G.mass.sorted[(c(as.numeric(round((nrow(G.mass.2.sorted)/3))))+1):(c(as.numeric(round((nrow(G.mass.2.sorted)))))-c(as.numeric(round((nrow(G.mass.2.sorted)/3))))), ]
# G.mass.low <- G.mass.sorted[((c(as.numeric(round((nrow(G.mass.2.sorted)))))-c(as.numeric(round((nrow(G.mass.2.sorted)/3)))))+1):c(as.numeric(round((nrow(G.mass.2.sorted))))), ]
# 
# plot(G.mass.high$time2.growth, G.mass.high$time3.growth, col="red", xlim=c(0,4), ylim=c(0,2.5), ylab="Growth mg cm-2 d-1",, xlab=" Growth mg cm-2 d-1", main="MARCH - MAY")
# points(G.mass.med$time2.growth, G.mass.med$time3.growth, col="blue")
# points(G.mass.low$time2.growth, G.mass.low$time3.growth, col="black")
# abline(v=G.mass.med[1,9], col="red", lwd=1, lty=2)
# abline(v=G.mass.low[1,9], col="blue", lwd=1, lty=2)
# abline(lm(G.mass.high$time3.growth ~ G.mass.high$time2.growth), col="red")
# abline(lm(G.mass.med$time3.growth ~ G.mass.med$time2.growth), col="blue")
# abline(lm(G.mass.low$time3.growth ~ G.mass.low$time2.growth), col="black")
# 
# ##### Time 3 as a function of Time 3 #####
# #empty
# plot(0,type='n',axes=FALSE,ann=FALSE)
# 
# ##### Time 4 as a function of Time 1 #####
# G.mass.sorted <- G.mass[order(-G.mass$time1.growth),]
# G.mass.high <- G.mass.sorted[1:c(as.numeric(round((nrow(G.mass.1.sorted)/3)))), ]
# G.mass.med <- G.mass.sorted[(c(as.numeric(round((nrow(G.mass.1.sorted)/3))))+1):(c(as.numeric(round((nrow(G.mass.1.sorted)))))-c(as.numeric(round((nrow(G.mass.1.sorted)/3))))), ]
# G.mass.low <- G.mass.sorted[((c(as.numeric(round((nrow(G.mass.1.sorted)))))-c(as.numeric(round((nrow(G.mass.1.sorted)/3)))))+1):c(as.numeric(round((nrow(G.mass.1.sorted))))), ]
# 
# 
# plot(G.mass.high$time1.growth, G.mass.high$time4.growth, col="red", xlim=c(0,4), ylim=c(0,2.5), ylab="Growth mg cm-2 d-1", xlab=" Growth mg cm-2 d-1")
# points(G.mass.med$time1.growth, G.mass.med$time4.growth, col="blue")
# points(G.mass.low$time1.growth, G.mass.low$time4.growth, col="black")
# mtext(side=2, line=4, "AUGUST-NOVEMBER", col="black", font=2, cex=0.9)
# abline(v=G.mass.med[1,8], col="red", lwd=1, lty=2)
# abline(v=G.mass.low[1,8], col="blue", lwd=1, lty=2)
# abline(lm(G.mass.high$time4.growth ~ G.mass.high$time1.growth), col="red")
# abline(lm(G.mass.med$time4.growth ~ G.mass.med$time1.growth), col="blue")
# abline(lm(G.mass.low$time4.growth ~ G.mass.low$time1.growth), col="black")
# 
# ##### Time 4 as a function of Time 2 #####
# G.mass.sorted <- G.mass[order(-G.mass$time2.growth),]
# G.mass.high <- G.mass.sorted[1:c(as.numeric(round((nrow(G.mass.2.sorted)/3)))), ]
# G.mass.med <- G.mass.sorted[(c(as.numeric(round((nrow(G.mass.2.sorted)/3))))+1):(c(as.numeric(round((nrow(G.mass.2.sorted)))))-c(as.numeric(round((nrow(G.mass.2.sorted)/3))))), ]
# G.mass.low <- G.mass.sorted[((c(as.numeric(round((nrow(G.mass.2.sorted)))))-c(as.numeric(round((nrow(G.mass.2.sorted)/3)))))+1):c(as.numeric(round((nrow(G.mass.2.sorted))))), ]
# 
# plot(G.mass.high$time2.growth, G.mass.high$time4.growth, col="red", xlim=c(0,4), ylim=c(0,2.5), ylab="Growth mg cm-2 d-1", xlab=" Growth mg cm-2 d-1")
# points(G.mass.med$time2.growth, G.mass.med$time4.growth, col="blue")
# points(G.mass.low$time2.growth, G.mass.low$time4.growth, col="black")
# abline(v=G.mass.med[1,9], col="red", lwd=1, lty=2)
# abline(v=G.mass.low[1,9], col="blue", lwd=1, lty=2)
# abline(lm(G.mass.high$time3.growth ~ G.mass.high$time4.growth), col="red")
# abline(lm(G.mass.med$time3.growth ~ G.mass.med$time4.growth), col="blue")
# abline(lm(G.mass.low$time3.growth ~ G.mass.low$time4.growth), col="black")
# 
# ##### Time 4 as a function of Time 3 #####
# G.mass.sorted <- G.mass[order(-G.mass$time3.growth),]
# G.mass.high <- G.mass.sorted[1:c(as.numeric(round((nrow(G.mass.3.sorted)/3)))), ]
# G.mass.med <- G.mass.sorted[(c(as.numeric(round((nrow(G.mass.3.sorted)/3))))+1):(c(as.numeric(round((nrow(G.mass.3.sorted)))))-c(as.numeric(round((nrow(G.mass.3.sorted)/3))))), ]
# G.mass.low <- G.mass.sorted[((c(as.numeric(round((nrow(G.mass.3.sorted)))))-c(as.numeric(round((nrow(G.mass.3.sorted)/3)))))+1):c(as.numeric(round((nrow(G.mass.3.sorted))))), ]
# 
# 
# plot(G.mass.high$time3.growth, G.mass.high$time4.growth, col="red", xlim=c(0,4), ylim=c(0,2.5), ylab="Growth mg cm-2 d-1", xlab=" Growth mg cm-2 d-1", main="MAY - AUGUST")
# points(G.mass.med$time3.growth, G.mass.med$time4.growth, col="blue")
# points(G.mass.low$time3.growth, G.mass.low$time4.growth, col="black")
# abline(v=G.mass.med[1,8], col="red", lwd=1, lty=2)
# abline(v=G.mass.low[1,8], col="blue", lwd=1, lty=2)
# abline(lm(G.mass.high$time4.growth ~ G.mass.high$time3.growth), col="red")
# abline(lm(G.mass.med$time4.growth ~ G.mass.med$time3.growth), col="blue")
# abline(lm(G.mass.low$time4.growth ~ G.mass.low$time3.growth), col="black")
# 
# ##### Time 4 as a function of Time 4 #####
# #empty
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #Survivorship Data
# G.mass <- data[,c(1:3,6,16,26,36,48,49,50)] #as.data.frame(cbind(data$Fragment.ID, data$Rack, data$Depth,as.character(data$Status.Time1), as.character(data$Status.Time2), as.character(data$Status.Time3), as.numeric(data$time1.growth), as.numeric(data$time2.growth), as.numeric(data$time3.growth)))
# G.mass.1 <- data[,c(1:3,16,48)]
# G.mass.2 <- data[,c(1:3,26,49)]
# G.mass.3 <- data[,c(1:3,36,50)]
# G.mass.sorted <- G.mass[order(-G.mass$time1.growth),]
# 
# #Time 1 Survivorship
# G.mass.high <- G.mass.sorted[1:c(as.numeric(round((nrow(G.mass.1.sorted)/3)))), ]
# G.mass.med <- G.mass.sorted[(c(as.numeric(round((nrow(G.mass.1.sorted)/3))))+1):(c(as.numeric(round((nrow(G.mass.1.sorted)))))-c(as.numeric(round((nrow(G.mass.1.sorted)/3))))), ]
# G.mass.low <- G.mass.sorted[((c(as.numeric(round((nrow(G.mass.1.sorted)))))-c(as.numeric(round((nrow(G.mass.1.sorted)/3)))))+1):c(as.numeric(round((nrow(G.mass.1.sorted))))), ]
# 
# par(mfrow=c(1,3))
# plot(G.mass.low$time1.growth, G.mass.low$Status.Time1, col="black",  ylab="Status, March-May", xlab=" Growth mg cm-2 d-1, Jan-March", main="Slow")
# plot(G.mass.med$time1.growth, G.mass.med$Status.Time1, col="blue",  ylab="Status, March-May", xlab=" Growth mg cm-2 d-1, Jan-March", main="Medium")
# plot(G.mass.high$time1.growth, as.factor(G.mass.high$Status.Time1), col="red",  ylab="Status, March-May", xlab=" Growth mg cm-2 d-1, Jan-March", main="Fast")
# 
# #Time 2
# G.mass.sorted <- G.mass[order(-G.mass$time2.growth),]
# G.mass.high <- G.mass.sorted[1:c(as.numeric(round((nrow(G.mass.2.sorted)/3)))), ]
# G.mass.med <- G.mass.sorted[(c(as.numeric(round((nrow(G.mass.2.sorted)/3))))+1):(c(as.numeric(round((nrow(G.mass.2.sorted)))))-c(as.numeric(round((nrow(G.mass.2.sorted)/3))))), ]
# G.mass.low <- G.mass.sorted[((c(as.numeric(round((nrow(G.mass.2.sorted)))))-c(as.numeric(round((nrow(G.mass.2.sorted)/3)))))+1):c(as.numeric(round((nrow(G.mass.2.sorted))))), ]
# 
# par(mfrow=c(3,1))
# plot(G.mass.high$time2.growth, G.mass.high$Status.Time2, col="red",  ylab="Status, May-August", xlab=" Growth mg cm-2 d-1, March-May", main="Fast")
# plot(G.mass.med$time2.growth, G.mass.med$Status.Time2, col="blue",  ylab="Status, May-August", xlab=" Growth mg cm-2 d-1, March-May", main="Medium")
# plot(G.mass.low$time2.growth, G.mass.low$Status.Time2, col="black",  ylab="Status, May-August", xlab=" Growth mg cm-2 d-1, March-May", main="Slow")
# 
# 
