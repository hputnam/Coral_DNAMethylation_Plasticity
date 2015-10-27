##Examination of DNA Methylation and Phenotypic Plasticity in Corals
#Data from Putnam NSF (NSF OCE PRF-1323822) and EPSCOR Y5 (NSF EPS-0903833)
#last modified 20151021
#See Readme file for details

rm(list=ls()) # removes all prior objects

#Read in required libraries
##### Include Versions of libraries
library("car") #levenes test
library("ggplot2") #plotting
library("plotrix") #plotting
library("reshape") #data shaping
library("plyr") #splitting, applying, and combining data
library("seacarb") #seawater carbonate chemistry
library("vegan") #calculating distance matrices
library("MetabolAnalyze") #scaling function
source('/Users/hputnam/Publications/In_Prep/Primary/Bulk_Methylation/Coral_DNAMethylation_Plasticity/R_Analysis/Scripts/opls.R') #OPLS DA analysis script by Paul Anderson used in Sogin et al 2014 (http://birg.cs.cofc.edu/index.php/O-PLS)
library("pracma") #OPLS DA analysis requirements
library("caret") #OPLS DA analysis requirements
require("gridExtra") #Arrange Plots for output

#Required Data files
#BM_Acclimation.csv
#BM_Field_Temp.csv
#BM_Tank_Temp.csv
#BM_Tank_light.csv
#BM_Daily_Measurements.csv
#BM_SWChem.csv
#BM_NBS_pH.csv
#BM_NMR_Data_0.04binned_truncated_nowater.csv
#BM_NMR_sample_info
#pH probe calibration files : Coral_DNAMethylation_Plasticity/R_Analysis/Data/pH_Calibration_Files/
#BM_Methylation.csv
#BM_Buoyant_Weight.csv


#############################################################
setwd("/Users/hputnam/Publications/In_Prep/Primary/Bulk_Methylation/Coral_DNAMethylation_Plasticity/R_Analysis/Data") #set working directory


#------------------------------------------------
##ACCLIMATION LIGHT AND TEMPERATURE ANALYSIS
#Data from period that corals were held in tank prior to experimental conditions
#load tank acclimation light and temp data
Acclim.data <- read.csv("BM_Acclimation.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

light.eq <-Acclim.data$Tank3.Light #Assign light column in dataframe
light <-(light.eq+48.339)/23.674 #Apply the cross calibration of the odyssey light to Licor cosine sensor standard 192SA cosine sensor

temp.eq <-Acclim.data$Tank3.Temp #Assign temperature column in dataframe
temp <-(temp.eq+0.4518)/1.0208 #Apply the cross calibration of temperature to standard logger #1

Acc <- data.frame(Acclim.data$Date.Time, light, temp) #combine light and temperature data

mydate <- strptime(Acc$Acclim.data.Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters
quarterhours <- format(as.POSIXct(mydate) ,format = "%H:%M") #set time as every 15 min
Acc.data <- cbind(Acc, quarterhours) #make a dataframe out of data and new times
Acc.data #View data
Acc.light.N <- sum(!is.na(Acc.data$light)) #Count sample size
Acc.light.N #View data

quarterly.light.mean <- aggregate(light ~ quarterhours, data=Acc.data, mean, na.rm=TRUE) #calculate mean of light for every 15 min interval
quarterly.light.se <- aggregate(light ~ quarterhours, data=Acc.data, std.error, na.rm=TRUE) #calculate standard error of the mean of light for every 15 min interval

light.means <- data.frame(quarterly.light.mean, quarterly.light.se$light) #combine mean and standard error results
colnames(light.means) <- c("Time", "mean", "se")  #rename columns to describe contents

Fig1 <- ggplot(light.means) + #Plot average diurnal cycle of light data
  geom_point(aes(x = Time, y = mean), colour="black") + #Plot points using time as the x axis, light as the Y axis and black points
  geom_errorbar(aes(x=Time, ymax=mean+se, ymin=mean-se), data=light.means) + #set values for standard error bars
  scale_x_discrete(breaks=c("0:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00","13:00", "14:00", "15:00", "16:00", "17:00", "18:00","19:00", "20:00", "21:00", "22:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("C") + #Label the graph with the main title
  xlab("Time") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~ 'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig1 #View figure

#------------------------------------------------
#ACCLIMATION TEMPERATURE ANALYSIS
#Tank Temperature Data
quarterly.temp.mean <- aggregate(temp ~ quarterhours, data=Acc.data, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.temp.se <- aggregate(temp ~ quarterhours, data=Acc.data, std.error, na.rm=TRUE) #calculate standard error of the mean of temperature for every 15 min interval
Acc.temp.N <- sum(!is.na(Acc.data$temp)) #Count sample size
temp.means <- data.frame(quarterly.temp.mean, quarterly.temp.se$temp) #combine mean and standard error results
colnames(temp.means) <- c("Time", "mean", "se")  #rename columns to describe contents

Fig2 <- ggplot(temp.means) + #Plot average diurnal cycle of temperature data
  geom_point(aes(x = Time, y = mean), colour="black") + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=Time, ymax=mean+se, ymin=mean-se), position=position_dodge(0.9), data=temp.means) + #set values for standard error bars and offset on the X axis for clarity
  scale_x_discrete(breaks=c("0:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00","13:00", "14:00", "15:00", "16:00", "17:00", "18:00","19:00", "20:00", "21:00", "22:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("B") + #Label the graph with the main title
  xlab("Time") + #Label the X Axis
  ylab("Temperature (°C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig2 #View figure

#------------------------------------------------
#Field Temperature Data
#load field collection/acclimation temp data
Field.data <- read.csv("BM_Field_Temp.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
mydate2 <- strptime(Field.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters
quarterhours2 <- format(as.POSIXct(mydate2) ,format = "%H:%M") #set time as every 15 min
Field.data <- cbind(Field.data, quarterhours2) #make a dataframe out of data and new times
Field.data #View data
min(Field.data$Temperature)
max(Field.data$Temperature)
quarterly.temp.mean2 <- aggregate(Temperature ~ quarterhours2, data=Field.data, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.temp.se2 <- aggregate(Temperature ~ quarterhours2, data=Field.data, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
Field.temp.N <- sum(!is.na(Field.data$Temperature)) #Count sample size
field.temp.means <- data.frame(quarterly.temp.mean2, quarterly.temp.se2$Temperature) #combine mean and standard error results
colnames(field.temp.means) <- c("Time", "mean", "se")  #rename columns to describe contents

#Plot average diurnal cycle of temp data
Fig3 <- ggplot(field.temp.means) + #Plot average diurnal cycle of temperature data
  geom_point(aes(x = Time, y = mean), colour="black") + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=Time, ymax=mean+se, ymin=mean-se), position=position_dodge(0.9), data=field.temp.means) + #set values for standard error bars and offset on the X axis for clarity
  scale_x_discrete(breaks=c("0:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00","13:00", "14:00", "15:00", "16:00", "17:00", "18:00","19:00", "20:00", "21:00", "22:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("A") + #Label the graph with the main title
  xlab("Time") + #Label the X Axis
  ylab("Temperature (°C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig3 #View figure

#------------------------------------------------
#TANK EXPERIMENTAL TEMPERATURE ANALYSIS
#Tank4 = (y4+0.5642)/1.0247 Cross Calibration to temperature standard logger #1
#Tank5 = (y5+0.2527)/1.0126 Cross Calibration to temperature standard logger #1
#Tank4=Serial Number 10489734,	Tank5=Serial Number 10489735	

#LOAD DATA
tank.data <- read.csv("BM_Tank_Temp.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
mydate.tanks <- strptime(tank.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters

Tank4 <-(tank.data$Tank4+0.5642)/1.0247 #Standardize temperature data to standard logger #1 Serial Number 
Tank5 <-(tank.data$Tank5+0.2527)/1.0126 #Standardize temperature data to standard logger #1 Serial Number 
tank.tempdata <-data.frame(mydate.tanks, Tank4, Tank5) #make a dataframe of temperature and time
colnames(tank.tempdata) <- c("Date.Time", "Tank4", "Tank5")
Tank4.temp.N <- sum(!is.na(tank.tempdata$Tank4)) #Count sample size
Tank5.temp.N <- sum(!is.na(tank.tempdata$Tank5)) #Count sample size

Fig4 <- ggplot(tank.tempdata, aes(Date.Time)) + #plot tank temperature data
  geom_line(aes(y = Tank4, colour="Ambient")) + #plot Temperature data as a line on the Y axis with date as the X axis 
  geom_line(aes(y = Tank5, colour="High")) +
  scale_colour_manual("Treatment", values = c("blue","red")) +
  xlab("Date") + #Label the X Axis
  ylab("Temperature °C") + #Label the Y Axis
  ggtitle("") + #label the main title
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig4 #View figure

tank.temp <- tank.tempdata[, 2:3] #subset the data to tank number and temperature data
tank.temp <- melt(tank.temp,  na.rm=TRUE) #rearrange the data in long format removing NA or missing data
colnames(tank.temp) <- c("Tank.Number", "Tank.Temperature") #rename the data columns

mean.tank.temp <-ddply(tank.temp, .(Tank.Number), summarize, #For each subset of a data frame, apply function then combine results into a data frame.
                  mean = (mean(Tank.Temperature)),       #take the average of the temp column summarized by Tank number
                  sem = (sd(Tank.Temperature)/sqrt(length(tank.temp)))) #calculate the SEM as the sd/sqrt of the count or data length
mean.tank.temp # display mean and sem temp levels
sem.tank.temp <-transform(mean.tank.temp, lower=mean-sem, upper=mean+sem) # add the upper and lower SEM values to the dataframe
sem.tank.temp #display mean and sem temp levels with lower and upper values for each tank

Fig5 <- ggplot(mean.tank.temp, aes(x = Tank.Number, y = mean)) + # plot mean temp by tank
  geom_point(data = mean.tank.temp, aes(y = mean), # points represent the mean
             colour = 'black', size = 5) + # plot black points of size 5
  geom_errorbar(aes(ymax=upper, ymin=lower), data=sem.tank.temp) + # plot mean and sem together
  xlab("Tanks") + #Label the X Axis
  ylab("Temperature °C") + #Label the Y Axis
  ggtitle("") + #Label the main title
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig5 #View figure

tanks.mean <- mean(tank.temp$Tank.Temperature) #Calculate the grand average of both tanks
tanks.mean #View data
tanks.se <- sd(tank.temp$Tank.Temperature)/sqrt(length(tank.temp$Tank.Temperature)) #Calculate the overal standard error of both tanks
tanks.se #View data


#Plotting diurnal cycles
tank.time <- format(as.POSIXct(mydate.tanks) ,format = "%H:%M") #Format time into only hours and minutes
tank.temperatures <- cbind(tank.tempdata, tank.time) #create a dataframe 
colnames(tank.temperatures) <- c("Date", "Tank4", "Tank5", "Time") #Rename columns to describe contents
tank.temperatures #View Data

quarterly.tank.temp.mean4 <- aggregate(Tank4 ~ Time, data=tank.temperatures, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.tank.temp.se4 <- aggregate(Tank4 ~ Time, data=tank.temperatures, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
quarterly.tank.temp.mean5 <- aggregate(Tank5 ~ Time, data=tank.temperatures, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.tank.temp.se5 <- aggregate(Tank5 ~ Time, data=tank.temperatures, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
tank.temp.means <- data.frame(quarterly.tank.temp.mean4, quarterly.tank.temp.se4$Tank4, quarterly.tank.temp.mean5$Tank5, quarterly.tank.temp.se5$Tank5) #combine mean and standard error results
colnames(tank.temp.means) <- c("Time", "Tank4.mean", "Tank4.se", "Tank5.mean", "Tank5.se")  #Rename columns to describe contents

Fig6 <- ggplot(tank.temp.means, aes(Time)) + # plot mean temp by tank
  geom_point(aes(y = Tank4.mean, colour="Ambient")) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank4.mean+Tank4.se, ymin=Tank4.mean-Tank4.se), position=position_dodge(0.9), data=tank.temp.means) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(y = Tank5.mean, colour="High")) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank5.mean+Tank5.se, ymin=Tank5.mean-Tank5.se), position=position_dodge(0.9), data=tank.temp.means) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("blue","red")) +
  scale_x_discrete(breaks=c("0:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00","13:00", "14:00", "15:00", "16:00", "17:00", "18:00","19:00", "20:00", "21:00", "22:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("A") + 
  xlab("Time") + #Label the X Axis
  ylab("Temperature (°C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig6 #View figure

#------------------------------------------------
#TANK EXPERIMENTAL LIGHT ANALYSIS

#load light data
tank.light.data <- read.csv("BM_Tank_light.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
mydate.light <- strptime(tank.light.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters

Tank4 <-(tank.light.data$Tank4-63.099)/19.652 #Apply the cross calibration of the odyssey light to Licor cosine sensor standard 192SA cosine sensor
Tank5 <-(tank.light.data$Tank5-49.729)/19.54 #Apply the cross calibration of the odyssey light to Licor cosine sensor standard 192SA cosine sensor

tank.light.data <-data.frame(mydate.light, Tank4, Tank5) #make a dataframe of light and time
colnames(tank.light.data) <- c("Date.Time", "Tank4", "Tank5") #Rename columns to describe contents

Tank4.light.N <- sum(!is.na(tank.light.data$Tank4)) #Count sample size
Tank5.light.N <- sum(!is.na(tank.light.data$Tank5)) #Count sample size

Fig7 <- ggplot(tank.light.data, aes(Date.Time)) + #plot tank light data
  geom_line(aes(y = Tank4, colour="Ambient")) + #plot light data as a line on the Y axis with date as the X axis 
  geom_line(aes(y = Tank5, colour="High")) + #plot light data as a line on the Y axis with date as the X axis for additional treatment
  scale_colour_manual("Treatment", values = c("blue","red")) + #color the treatment groups differently
  xlab("Date") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~ 'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  ggtitle("BM Tank Light") + #Label the main title
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig7 #View figure

tank.light <- tank.light.data[, 2:3] #subset the data to tank number and light data
tank.light  <- melt(tank.light ,  na.rm=TRUE) #rearrange the data in long format removing NA or missing data
colnames(tank.light ) <- c("Tank.Number", "Tank.Light") #rename the data columns

mean.tank.light <-ddply(tank.light, .(Tank.Number), summarize, #For each subset of a data frame, apply function then combine results into a data frame.
                       mean = (mean(Tank.Light)),       #take the average of the light column summarized by Tank number
                       sem = (sd(Tank.Light)/sqrt(length(tank.light)))) #calculate the SEM as the sd/sqrt of the count or data length
mean.tank.light # display mean and sem temp levels
mean.tank.light<-transform(mean.tank.light, lower=mean-sem, upper=mean+sem) # add the upper and lower SEM values to the dataframe
mean.tank.light #display mean and sem of light levels with lower and upper values for each tank

Fig8 <- ggplot(mean.tank.light, aes(x = Tank.Number, y = mean)) + # plot mean light by tank
  geom_point(data = mean.tank.light, aes(y = mean), # points represent the mean
             colour = 'black', size = 5) + # plot black points of size 5
  geom_errorbar(aes(ymax=upper, ymin=lower), data=mean.tank.light) + # plot mean and sem together
  xlab("Tanks") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~ 'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  ggtitle("BM Diurnal Tank Light") + #label the main title
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank()) #Set the plot background
Fig8 #View figure

tanks.light.mean <- mean(tank.light$Tank.Light) #Calculate the grand average of light data
tanks.light.mean #View data
tanks.light.se <- sd(tank.light$Tank.Light)/(sqrt(length(tank.temp$Tank.Light))) #Calculate the overall standard error of the light data
tanks.light.se #View data

#Plotting diurnal cycles
tank.time <- format(as.POSIXct(mydate.light) ,format = "%H:%M") #set time to hours and minutes
tank.lights <- cbind(tank.light.data, tank.time) #Combine time and light data
colnames(tank.lights) <- c("Date", "Tank4", "Tank5", "Time") #Rename columns to describe contents
tank.lights #View data

quarterly.tank.light.mean4 <- aggregate(Tank4 ~ Time, data=tank.lights, mean, na.rm=TRUE) #calculate mean of light for every 15 min interval
quarterly.tank.light.se4 <- aggregate(Tank4 ~ Time, data=tank.lights, std.error, na.rm=TRUE)  #calculate standard error of the mean of light for every 15 min interval
quarterly.tank.light.mean5 <- aggregate(Tank5 ~ Time, data=tank.lights, mean, na.rm=TRUE) #calculate mean of light for every 15 min interval
quarterly.tank.light.se5 <- aggregate(Tank5 ~ Time, data=tank.lights, std.error, na.rm=TRUE)  #calculate standard error of the mean of light for every 15 min interval
tank.light.means <- data.frame(quarterly.tank.light.mean4, quarterly.tank.light.se4$Tank4, quarterly.tank.light.mean5$Tank5, quarterly.tank.light.se5$Tank5) #combine mean and standard error results
colnames(tank.light.means) <- c("Time", "Tank4.mean", "Tank4.se", "Tank5.mean", "Tank5.se")  #Rename columns to describe contents

Fig9 <- ggplot(tank.light.means, aes(Time)) + # plot mean light by tank
  geom_point(aes(y = Tank4.mean, colour="Ambient")) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank4.mean+Tank4.se, ymin=Tank4.mean-Tank4.se), position=position_dodge(0.9), data=tank.light.means) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(y = Tank5.mean, colour="High")) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank5.mean+Tank5.se, ymin=Tank5.mean-Tank5.se), position=position_dodge(0.9), data=tank.light.means) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("blue","red")) +
  scale_x_discrete(breaks=c("0:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00","13:00", "14:00", "15:00", "16:00", "17:00", "18:00","19:00", "20:00", "21:00", "22:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("B") + #Label the graph
  xlab("Time") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~ 'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig9 #View figure

#------------------------------------------------
#SEAWATER CHEMISTRY ANALYSIS FOR CONTINUOUS MEASUREMENTS

#read in probe measurements of pH, temperature, and salinity from tanks
daily <- read.csv("BM_Daily_Measurements.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

sal <- aggregate(Salinity ~ Treatment, data=daily, mean, na.rm = TRUE) #Calculate average salinity from the tank YSI probe measurements
sal #View data
#if salinity is the same in both treatments, take grand average
#if the salinity is different by treatment, use individual salinity values manually assigned
sal <- mean(daily$Salinity) #Assign the salinity from the calculated mean value
sal #View data

# read in total alkalinity, temperature, and salinity
SW.chem <- read.csv("BM_SWChem.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

TA.mean <- aggregate(Corrected.TA ~ Treatment, data=SW.chem, mean, na.rm = TRUE) #calculate the average of total alkalinity in each tank
TA.mean #View data
TA.se <- aggregate(Corrected.TA ~ Treatment, data=SW.chem, std.error, na.rm = TRUE) #calculate the standard error of total alkalinity in each tank
TA.se #View data
TAs <- cbind(TA.mean, TA.se$Corrected.TA) #merge the mean and standard error of total alkalinity
colnames(TAs) <- c("Treatment", "mean", "se") #Rename columns to describe contents
TAs #View data

# read in NBS pH data from Aquacontrollers frequency 15min
pHs <- read.csv("BM_NBS_pH.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

mydate.pHs <- strptime(pHs$Date.Time, format="%m/%d/%y %H:%M") #Identify date format
pHs <- data.frame(mydate.pHs, pHs$High.NBS, pHs$Ambient.NBS) #combine dataframe of date and pH
colnames(pHs) <- c("Date.Time", "High.NBS", "Ambient.NBS") #Rename columns to describe contents
pH.data <- merge(pHs, tank.tempdata, by="Date.Time") #merge data sets by time and date

#Convert pH from NBS to total scale for high tank 
pH <- pH.data$High.NBS # pH data on NBS scale logged every 15 minutes
Temp.High <- pH.data$Tank5+273.15 # Temperature from Hobo loggers every 15 minutes that match pH measurement frequency converted to Kelvin
S <- sal # salinity measured daily with YSI meter calibrated to conductivity standard at 25°C

TS <- (0.14/96.062)*(S/1.80655) #concentration of SO4-2 in mol/kg-SW 
##Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966: this is .02824.*Sali./35. = .0008067.*Sali

TF <- (0.000067/18.998)*(S/1.80655) #concentration of fluoride in mol/kg-SW 
##Riley, J. P., Deep-Sea Research 12:219-220, 1965

fH <- 1.2948 - 0.002036*Temp.High + (0.0004607 -  0.000001475*Temp.High)*(S^2) # the activity coefficient of the H+ ion, which is valid for the temperatures of 20-40°C
#Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition, v. 3, 1982 (p. 80)

IonS <- (19.924*S)/(1000 - 1.005*S) # the ionic strength of Hydrogen Fluoride 
##This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4 and Zeebe Wolfgladrow Appendix A p260

lnKSO4 <- -4276.1/Temp.High + 141.328 - 23.093*log(Temp.High) + 
  (-13856/Temp.High + 324.57 - 47.986*log(Temp.High))*sqrt(IonS) +    
  (35474/Temp.High - 771.54 + 114.723*log(Temp.High))*IonS +
  (-2698/Temp.High)*sqrt(IonS)*IonS + (1776/Temp.High)*IonS^2 
KSO4 <- exp(lnKSO4)*(1-0.001005*S) #this is on the free pH scale in mol/kg-SW 
#This is from the DOE handbook 1994 ORNL/CDIAC-74

lnKF <- 1590.2/Temp.High - 12.641 + 1.525*sqrt(IonS)
KF <- exp(lnKF)*(1-0.001005*S) #the dissociation constant of HF, this is on the free pH scale in mol/kg-H2O
#Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979

#Conversion Functions
SWStoTOT.high  <- (1 + TS/KSO4)/(1 + TS/KSO4 + TF/KF) #pH scale conversion factor
NBStoTOT.high <- pH-(log(SWStoTOT.high)/log(0.1) + log(fH)/log(0.1))  #conversion for NBS to total                                  

pH.High.Out <- data.frame(pH.data$Date.Time, pH.data$Tank5, pH.data$High.NBS, NBStoTOT.high) #create dataframe of total pH, Temperature, and time
colnames(pH.High.Out) <- c("Date.Time", "High.Temp", "pH.High.NBS", "pH.High.TOTAL") #Rename columns to describe contents

#Convert pH from NBS to total scale for ambient tank 
pH<- pH.data$Ambient.NBS # pH data on NBS scale logged every 15 minutes
Temp.Amb <- pH.data$Tank4+273.15 # Temperature from Hobo loggers every 15 minutes that match pH measurement frequency

TS <- (0.14/96.062)*(S/1.80655) #concentration of SO4-2 in mol/kg-SW 
##Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966: this is .02824.*Sali./35. = .0008067.*Sali

TF <- (0.000067/18.998)*(S/1.80655) #concentration of fluoride in mol/kg-SW 
##Riley, J. P., Deep-Sea Research 12:219-220, 1965

fH <- 1.2948 - 0.002036*Temp.Amb + (0.0004607 -  0.000001475*Temp.Amb)*(S^2) # the activity coefficient of the H+ ion, which is valid for the temperatures of 20-40°C
#Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition, v. 3, 1982 (p. 80)

IonS <- (19.924*S)/(1000 - 1.005*S) # the ionic strength of Hydrogen Fluoride 
##This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4 and Zeebe Wolfgladrow Appendix A p260

lnKSO4 <- -4276.1/Temp.Amb + 141.328 - 23.093*log(Temp.Amb) + 
  (-13856/Temp.Amb + 324.57 - 47.986*log(Temp.Amb))*sqrt(IonS) +    
  (35474/Temp.Amb - 771.54 + 114.723*log(Temp.Amb))*IonS +
  (-2698/Temp.Amb)*sqrt(IonS)*IonS + (1776/Temp.Amb)*IonS^2 
KSO4 <- exp(lnKSO4)*(1-0.001005*S) #this is on the free pH scale in mol/kg-SW 
#This is from the DOE handbook 1994 ORNL/CDIAC-74

lnKF <- 1590.2/Temp.Amb - 12.641 + 1.525*sqrt(IonS)
KF <- exp(lnKF)*(1-0.001005*S) #the dissociation constant of HF, this is on the free pH scale in mol/kg-H2O
#Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979

#Conversion Functions
SWStoTOT.amb  <- (1 + TS/KSO4)/(1 + TS/KSO4 + TF/KF) #pH scale conversion factor
NBStoTOT.amb <- pH-(log(SWStoTOT.amb)/log(0.1) + log(fH)/log(0.1))  #conversion for NBS to total                                  

pH.Amb.Out <- data.frame(pH.data$Date.Time, pH.data$Tank4, pH.data$Ambient.NBS, NBStoTOT.amb) #create dataframe of pH and time
colnames(pH.Amb.Out) <- c("Date.Time", "Amb.Temp", "pH.Amb.NBS", "pH.Amb.TOTAL") #Rename columns to describe contents

pH.TOT <- cbind(pH.High.Out, pH.Amb.Out$Amb.Temp,  pH.Amb.Out$pH.Amb.NBS,	pH.Amb.Out$pH.Amb.TOTAL) #create dataframe of total pH and time
colnames(pH.TOT) <- c("Date.Time",  "High.Temp",	"pH.High.NBS",	"pH.High.TOTAL",	"Amb.Temp",	"pH.Amb.NBS",	"pH.Amb.TOTAL")

#Plot total pH for both treatments for duration of experiment
Fig10 <- ggplot(pH.TOT) + #plot pH total scale
  geom_line(aes(x = Date.Time, y = pH.High.TOTAL, col="High")) + #plot as a line
  geom_line(aes(x = Date.Time, y = pH.Amb.TOTAL, col="Ambient")) + #plot as a line
  xlab("Date") + #Label the X Axis
  ylab("pH (Total Scale)") + #Label the Y Axis
  ggtitle("All Tanks Total pH") + #Label the graph title
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig10 #View figure

pH.time <- format(as.POSIXct(pH.TOT$Date.Time) ,format = "%H:%M") #Identify time in hours and minutes
pH.TOT <- cbind(pH.TOT, pH.time) #Combine total pH and time
colnames(pH.TOT) <- c("Date.Time",  "High.Temp",  "pH.High.NBS",	"pH.High.TOTAL",	"Amb.Temp",	"pH.Amb.NBS",	"pH.Amb.TOTAL", "Time") #Rename columns to describe contents

quarterly.tank.pH.amb.mean <- aggregate(pH.Amb.TOTAL ~ Time, data=pH.TOT, mean, na.rm=TRUE) #calculate mean of pH for every 15 min interval
quarterly.tank.pH.amb.se <- aggregate(pH.Amb.TOTAL ~ Time, data=pH.TOT, std.error, na.rm=TRUE)  #calculate standard error of the mean of pH for every 15 min interval
quarterly.tank.pH.high.mean <- aggregate(pH.High.TOTAL ~ Time, data=pH.TOT, mean, na.rm=TRUE) #calculate mean of pH for every 15 min interval
quarterly.tank.pH.high.se <- aggregate(pH.High.TOTAL ~ Time, data=pH.TOT, std.error, na.rm=TRUE)  #calculate standard error of the mean of pH for every 15 min interval
tank.pH.means <- data.frame(quarterly.tank.pH.amb.mean, quarterly.tank.pH.amb.se$pH.Amb.TOTAL, quarterly.tank.pH.high.mean$pH.High.TOTAL, quarterly.tank.pH.high.se$pH.High.TOTAL) #combine mean and standard error results
colnames(tank.pH.means) <- c("Time", "pH.Amb.TOTAL.mean", "pH.Amb.TOTAL.se", "pH.High.TOTAL.mean", "pH.High.TOTAL.se")  #Rename columns to describe contents

Fig11 <- ggplot(tank.pH.means, aes(Time)) + # plot mean pH by tank
  geom_point(aes(x =Time, y = pH.Amb.TOTAL.mean, colour="Ambient")) + #plot points
  geom_errorbar(aes(x=Time, ymax=pH.Amb.TOTAL.mean+pH.Amb.TOTAL.se, ymin=pH.Amb.TOTAL.mean-pH.Amb.TOTAL.se), position=position_dodge(0.9), data=tank.pH.means) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(x = Time, y = pH.High.TOTAL.mean, colour="High")) + #plot points
  geom_errorbar(aes(x=Time, ymax=pH.High.TOTAL.mean+pH.High.TOTAL.se, ymin=pH.High.TOTAL.mean-pH.High.TOTAL.se), position=position_dodge(0.9), data=tank.pH.means) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("blue","red")) +
  scale_x_discrete(breaks=c("0:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00","13:00", "14:00", "15:00", "16:00", "17:00", "18:00","19:00", "20:00", "21:00", "22:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("A") + #Label graph
  xlab("Time") + #Label the X Axis
  ylab("pH (Total Scale)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig11 #View figure

#Generate dataframe and set parameters for seacarb calculations
For.seacarb <- pH.TOT[,c(1,2,4,5,7,8)] #create data frame of temp and total pH
For.seacarb$Salinity <- sal #set salinity from average of tanks
For.seacarb$TA.Amb <- TAs[1,2] #set TA from average of ambient tanks calculated above
For.seacarb$TA.High <- TAs[2,2] #set TA from average of high tanks calculated above
colnames(For.seacarb) <- c("Date.Time", "Temp.High", "pH.High.TOTAL",	"Temp.Amb", "pH.Amb.TOTAL", "Time", "Salinity",	"TA.Amb",	"TA.High") #Rename columns to describe contents
For.seacarb <- na.omit(For.seacarb) #omit missing values

#Calculate CO2 parameters using seacarb
carbo.high <- carb(flag=8, var1=For.seacarb$pH.High.TOTAL, var2=For.seacarb$TA.High/1000000, S= For.seacarb$Salinity, T=For.seacarb$Temp.High, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d") #calculate seawater chemistry parameters using seacarb
carbo.amb <- carb(flag=8, var1=For.seacarb$pH.Amb.TOTAL, var2=For.seacarb$TA.Amb/1000000, S= For.seacarb$Salinity, T=For.seacarb$Temp.Amb, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d") #calculate seawater chemistry parameters using seacarb

pCO2 <- data.frame(For.seacarb$Date.Time, For.seacarb$Time, carbo.amb$pCO2, carbo.high$pCO2) #make dataframe of CO2 output
colnames(pCO2) <- c("Date.Time", "Time", "pCO2.Amb", "pCO2.High") #Rename columns to describe contents

quarterly.tank.pCO2.amb.mean <- aggregate(pCO2.Amb ~ Time, data=pCO2, mean, na.rm=TRUE) #calculate mean of pCO2 for every 15 min interval
quarterly.tank.pCO2.amb.se <- aggregate(pCO2.Amb ~ Time, data=pCO2, std.error, na.rm=TRUE)  #calculate standard error of the mean of pCO2 for every 15 min interval
quarterly.tank.pCO2.high.mean <- aggregate(pCO2.High ~ Time, data=pCO2, mean, na.rm=TRUE) #calculate mean of pCO2 for every 15 min interval
quarterly.tank.pCO2.high.se <- aggregate(pCO2.High ~ Time, data=pCO2, std.error, na.rm=TRUE)  #calculate standard error of the mean of pCO2 for every 15 min interval
tank.pCO2.means <- data.frame(quarterly.tank.pCO2.amb.mean, quarterly.tank.pCO2.amb.se$pCO2.Amb, quarterly.tank.pCO2.high.mean$pCO2.High, quarterly.tank.pCO2.high.se$pCO2.High) #combine mean and standard error results
colnames(tank.pCO2.means) <- c("Time", "Amb.mean", "Amb.se", "High.mean", "High.se") #Rename columns to describe contents

#plotting averages of total pH every 15 min for a 1 day cycle from all data
Fig12 <- ggplot(tank.pCO2.means) + #plot pCO2
  geom_point(aes(x = Time, y = Amb.mean, colour="Ambient")) + #plot as line
  geom_errorbar(aes(x=Time, ymax=Amb.mean+Amb.se, ymin=Amb.mean-Amb.se), position=position_dodge(0.9), data=tank.pCO2.means) + #plot error bars
  geom_point(aes(x = Time, y = High.mean, colour="High")) + #plot as line
  geom_errorbar(aes(x=Time, ymax=High.mean+High.se, ymin=High.mean-High.se), position=position_dodge(0.9), data=tank.pCO2.means) + #plot error bars
  scale_colour_manual("Treatment", values = c("blue","red")) + #Set colors for treatments
  scale_x_discrete(breaks=c("0:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00","13:00", "14:00", "15:00", "16:00", "17:00", "18:00","19:00", "20:00", "21:00", "22:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("B") + #Label graph
  xlab("Time") + #Label the X Axis
  ylab(expression(paste('p',CO[2], ' (µatm)', sep=''))) + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig12 #View figure

#------------------------------------------------
#SEAWATER CHEMISTRY ANALYSIS FOR DISCRETE MEASUREMENTS
#Seawater chemistry table from simultaneous TA, pH, temperature and salinity measurements

#Need to load data to make conversion equations for pH from mV to total scale using tris standard

path <-("/Users/hputnam/Publications/In_Prep/Primary/Bulk_Methylation/Coral_DNA_Methylation/R_Analysis/Data/pH_Calibration_Files/")

#list all the file names in the folder to get only get the csv files
file.names<-list.files(path = path, pattern = "csv$")

pH.cals <- data.frame(matrix(NA, nrow=length(file.names), ncol=3, dimnames=list(file.names,c("Date", "Intercept", "Slope")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Calib.Data <-read.table(file.path(path,file.names[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
  model <-lm(mVTris ~ TTris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
  coe <- coef(model) #extracts the coeffecients
  pH.cals[i,2:3] <- coe #inserts them in the dataframe
  pH.cals[i,1] <- substr(file.names[i],1,8) #stores the file name in the Date column
}

colnames(pH.cals) <- c("Calib.Date",  "Intercept",	"Slope")

#merge with Seawater chemistry file
SW.chem <- merge(pH.cals, SW.chem, by="Calib.Date")

#constants for use in pH calculation 
R <- 8.31447215 #gas constant in J mol-1 K-1 
F <-96485.339924 #Faraday constant in coulombs mol-1

mvTris <- SW.chem$Temperature*SW.chem$Slope+SW.chem$Intercept #calculate the mV of the tris standard using the temperature mv relationships in the measured standard curves 
STris<-34.5 #salinity of the Tris
phTris<- (11911.08-18.2499*STris-0.039336*STris^2)*(1/(SW.chem$Temperature+273.15))-366.27059+ 0.53993607*STris+0.00016329*STris^2+(64.52243-0.084041*STris)*log(SW.chem$Temperature+273.15)-0.11149858*(SW.chem$Temperature+273.15) #calculate the pH of the tris (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
SW.chem$pH.Total<-phTris+(mvTris/1000-SW.chem$pH.MV/1000)/(R*(SW.chem$Temperature+273.15)*log(10)/F) #calculate the pH on the total scale (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)

#Calculate CO2 parameters using seacarb
carb.ouptput <- carb(flag=8, var1=SW.chem$pH.Total, var2=SW.chem$Corrected.TA/1000000, S= SW.chem$Salinity, T=SW.chem$Temperature, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d") #calculate seawater chemistry parameters using seacarb

carb.ouptput$ALK <- carb.ouptput$ALK*1000000 #convert to µmol kg-1
carb.ouptput$CO2 <- carb.ouptput$CO2*1000000 #convert to µmol kg-1
carb.ouptput$HCO3 <- carb.ouptput$HCO3*1000000 #convert to µmol kg-1
carb.ouptput$CO3 <- carb.ouptput$CO3*1000000 #convert to µmol kg-1
carb.ouptput$DIC <- carb.ouptput$DIC*1000000 #convert to µmol kg-1

carb.ouptput <- cbind(SW.chem$Measure.Date,  SW.chem$Tank,	SW.chem$Treatment, carb.ouptput) #combine the sample information with the seacarb output
colnames(carb.ouptput) <- c("Date",  "Tank",  "Treatment",	"flag",	"Salinity",	"Temperature",	"Pressure",	"pH",	"CO2",	"pCO2",	"fCO2",	"HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat", 	"Calcite.Sat") #Rename columns to describe contents

carbo.melted <- melt(carb.ouptput) #reshape the dataframe to more easily summarize all output parameters
mean.carb.output <-ddply(carbo.melted, .(Treatment, variable), summarize, #For each subset of a data frame, apply function then combine results into a data frame.
                       mean = (mean(value)),       #take the average of the parameters (variables) summarized by treatments
                       sem = (sd(value)/sqrt(length(carbo.melted)))) #calculate the SEM as the sd/sqrt of the count or data length
mean.carb.output # display mean and sem 
mean.carb.output <- mean.carb.output[with(mean.carb.output, order(variable)), ] #order the data by the variables
mean.carb.output <- mean.carb.output[-c(1:4,9,10,29,30), ] #remove non-numeric parameters extra

write.table (mean.carb.output, "Seawater_chemistry_table_Output.csv", sep=",", row.names = FALSE)

#------------------------------------------------
#METABOLOMIC ANALYSIS

#load data after truncation from 0.5-10 and removal of water peak (4.73959- 4.93955 ppm)
metabo.data <- read.csv("BM_NMR_Data_0.04binned_truncated_nowater.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

#loadsample info
metabo.info <- read.csv("BM_NMR_sample_info.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

#DATA NORMALIZATION AND SCALING
NMRData.norm<-sweep(as.matrix(metabo.data[,grep('X', colnames(metabo.data))]), 1,metabo.info$Extract.Weight.g, '/') #normalize by extractweight
NMRData.norm[NMRData.norm<0]<-0 #Assign negative values to 0 value

#Scale the data
NMRData.scale <- scaling(NMRData.norm, type="pareto") #scale the data using pareto scaling to give more weight to medium features without inflating the baseline noise, the calculate is done by dividing each variable by the square root of its SD

#Mean center the data
NMRData.scale <- scale(NMRData.scale, center=TRUE, scale=FALSE) #scale the data about the mean. This calculation is done by subtracting the grand mean from each individual data point

#calculate similarity matrix of euclidean distance
dist.euc <- vegdist(NMRData.scale, method = "euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) #calculates euclidean distance on the transformed data matrix

#Run PCA and look at summary to identify outliers
pca.res <- prcomp(dist.euc, retx=TRUE)
summary(pca.res)
biplot(pca.res)
Y <- metabo.info$Treatment #SELECTS FOR GROUPING
ordiellipse(pca.res, group=rep(c(1), c(length(Y))), conf=0.99, kind='sd', draw='lines')  

#remove outliers identified outside 99%cl
cleaned.all.data <- NMRData.norm[-c(3,4,7,11,15,16,19,41,42,44,48,51,53),] #remove outliers identified in PCA outside 99%cl
cleaned.all.info <- metabo.info[-c(3,4,7,11,15,16,19,41,42,44,48,51,53),] #remove sample info of outliers identified in PCA outside 99%cl

cleaned.data <- NMRData.scale[-c(3,4,7,11,15,16,19,41,42,44,48,51,53),] #remove outliers identified in PCA outside 99%cl
cleaned.info <- metabo.info[-c(3,4,7,11,15,16,19,41,42,44,48,51,53),] #remove sample info of outliers identified in PCA outside 99%cl

metabo.Species <- cleaned.info$Species #identify Species
metabo.Treatment <- cleaned.info$Treatment #Identify Treatments

sample.counts <- aggregate(cleaned.info["Coral.ID"], by=cleaned.info[c("Species","Treatment")], FUN=length) #calculate sample size 
sample.counts #View data

#REcalculate similarity matrix of euclidean distance of only cleaned data
dist.euc <- vegdist(cleaned.data, method = "euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) #calculates euclidean distance on the transformed data matrix

#Run PCA and look at summary to identify outliers
pca.res <- prcomp(dist.euc, retx=TRUE) #run PCA analysis
summary(pca.res) #summarize PCA output
biplot(pca.res) #plot PCs and loadings
Treats <- cleaned.info$Treatment #identifies the grouping information
ordiellipse(pca.res,group=rep(c(1), c(length(Treats))), conf=0.99, kind='sd', draw='lines')  

scores <- as.data.frame(pca.res$x) # extract all PCs

PCA.plot <- ggplot(scores, aes(x = PC1, y = PC2, shape=metabo.Species, col=metabo.Treatment)) + # plot of PC1 and PC2
  geom_point(size=4) + # use points + 
  scale_shape_manual(values=c(16,17, 18, 3)) + #create your own scale and values for scale
  theme(panel.background = element_rect(fill='white', colour='black')) + #removes gray background
  theme(panel.grid.major = element_line(colour="white"), panel.grid.minor = element_line(colour="white")) + #Set the panel gridlines
  theme(legend.key = element_rect(fill='white')) + #Set plot legend key
  theme(panel.border = element_rect(fill=NA, size = 1, colour = "black")) #Set the border
PCA.plot

#OPLS Analysis of All Data
ALL.X<-as.matrix(cleaned.data) #Identify Data
class(ALL.X)<-'numeric' #changes class data to numeric 
ALL.Y<-as.matrix(as.numeric(as.factor(as.vector(cleaned.info$Species)))) #Select the grouping for model testing
ALL.ulabels <- unique(as.vector(cleaned.info$Species)) #Identify the labels for the groupings
set.seed(10) #sets the value on the uniform random-number generator so that the analysis can be reproduced every time
ALL.resultsALL<-n.group.opls(ALL.X,ALL.Y, num_permutations=100, CV= 10, nIterations=100, min_num_OPLS_fact=1) #run the OPLS DA model (Anderson)
ALL.Q2<-ALL.resultsALL$Q2 #GIVES Q2 VALUE
ALL.pval<-ALL.resultsALL$helper.results$pvalue #GIVES P-VALUE FOR MODEL 
ALL.modelALL <- ALL.resultsALL$helper.results$model
ALL.OPLSResults<-data.frame(model=c('Species'), Q2=c(ALL.Q2), pval=c(ALL.pval))
ALL.OPLSResults

#visualize
t.opls <- as.numeric(as.vector(ALL.modelALL$t)) #assign numeric t values from OPLS results
t.ortho.opls <- as.numeric(as.vector(ALL.modelALL$t_ortho)) #assign numeric t-orthogonal values from OPLS results
All.metabo.data <- data.frame(metabo.Species, t.opls, t.ortho.opls) #create a dataframe
colnames(All.metabo.data) <- c("Species", "t", "t.ortho") #Rename columns

Fig13 <- ggplot(All.metabo.data, aes(t, t.ortho)) + 
  geom_point(aes(x=t, y=t.ortho, colour=Species, shape=Species)) +
  scale_colour_manual("Species", values = c("black","gray")) + #Set colors for Species
  scale_shape_manual(values = c(17,19)) + #Set shape for Species
  xlim(-1600, 1600) +
  ylim(-2500, 2500) +
  ggtitle("A All Species") + #Label graph
  xlab("t") + #Label the X Axis
  ylab("t-Orthogonal") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        legend.position="top", #Set position of legend in graphic
        plot.title=element_text(hjust=0, face="bold")) #Justify the title to the top left and use Bold italics font
Fig13  #View figure


#MONTIPORA OUTLIER CHECK
MC.data <-subset(cleaned.all.data, cleaned.all.info$Species=='Montipora capitata') #Select only the Montipora data from the dataset from which outliers have been removed.
MC.info<-subset(cleaned.all.info, cleaned.all.info$Species=='Montipora capitata') #Select only the Montipora info
which(colSums(MC.data) ==0) #check to see if any bins (columns) are present in the whole analysis, but not in Montipora only
MC.data <- MC.data [, colSums(MC.data  != 0) > 0]
which(colSums(MC.data) ==0)

#assign treatment from info
MC.Treatment <- MC.info$Treatment

#Scale the data
MC.NMRData.scale <-scaling(MC.data, type="pareto")

#Mean center the data
MC.NMRData.scale <- scale(MC.NMRData.scale, center=TRUE, scale=FALSE)

#calculate similarity matrix of euclidean distance
MC.dist.euc <- vegdist(MC.NMRData.scale, method = "euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = TRUE) #calculates euclidean distance on the transformed data matrix
MC.pca.res <- prcomp(MC.dist.euc, retx=TRUE)
summary(MC.pca.res)

#check for outliers
biplot(MC.pca.res) #plot PC1 and PC2
ordiellipse(MC.pca.res,group=rep(c(1), c(length(MC.Treatment))), conf=0.99, kind='sd', draw='lines')  #plot 99% confidence interval ellipse

#remove outliers and assess sample size
cleaned.MC.data <- MC.NMRData.scale[-c(5,37,22,9,13),] #remove outliers identified in PCA
cleaned.MC.info <- MC.info[-c(5,37,22,9,13),] #remove sample info of outliers identified in PCA
cleaned.MC.Treat <- cleaned.MC.info$Treatment #Set the treatment information for the filtered data
MC.counts <- aggregate(cleaned.MC.info["Coral.ID"], by=cleaned.MC.info[c("Species","Treatment")], FUN=length) #check the sample size for each Treatment
MC.counts #view the sample size for each Treatment

#calculate similarity matrix of euclidean distance
MC.dist.euc <- vegdist(cleaned.MC.data, method = "euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = TRUE) #calculates euclidean distance on the transformed data matrix
MC.pca.res <- prcomp(MC.dist.euc, retx=TRUE)
summary(MC.pca.res)

#Double check for outliers
biplot(MC.pca.res)
ordiellipse(MC.pca.res,group=rep(c(1), c(length(cleaned.MC.Treat))), conf=0.99, kind='sd', draw='lines')  

# extract all PCs
MC.scores <- as.data.frame(MC.pca.res$x)
MC.PCA.plot <- ggplot(MC.scores, aes(x = PC1, y = PC2, shape=cleaned.MC.Treat , col=cleaned.MC.Treat )) +
  geom_point(size=4) + # use points + 
  scale_shape_manual(values=c(16,16)) + #create your own scale and values for scale
  theme(panel.background = element_rect(fill='white', colour='black')) + #removes gray background
  theme(panel.grid.major = element_line(colour="white"), panel.grid.minor = element_line(colour="white")) + #Set the panel gridlines
  theme(legend.key = element_rect(fill='white')) + #Set plot legend key
  theme(panel.border = element_rect(fill=NA, size = 1, colour = "black")) + #Set the border
  ggtitle("M. capitata")
MC.PCA.plot

#Montipora OPLS Data
MC.X<-as.matrix(cleaned.MC.data) #SELECT DATA -- USER INPUT
class(MC.X)<-'numeric' #changes class data to numeric 
MC.Y<-as.matrix(as.numeric(as.factor(as.vector(cleaned.MC.info$Treatment)))) #SELECT GROUPING -- USER INPUT BY CHANGING COLUMN USED
MC.ulabels <- unique(as.vector(cleaned.MC.info$Treatment)) #PICK LABELS 
set.seed(10) #sets the value on the uniform random-number generator so that the analysis can be reproduced every time
MC.resultsALL<-n.group.opls(MC.X,MC.Y, num_permutations=100, CV= 10, nIterations=100, min_num_OPLS_fact=1) #RUNS MODEL
MC.Q2<-MC.resultsALL$Q2 #GIVES Q2 VALUE
MC.pval<-MC.resultsALL$helper.results$pvalue #GIVES P-VALUE FOR MODEL 
MC.modelALL <- MC.resultsALL$helper.results$model #lists results
MC.OPLSResults<-data.frame(model=c('M. capitata'), Q2=c(MC.Q2), pval=c(MC.pval)) #lists Q2 and p
MC.OPLSResults

#visualize
MC.t.opls <- as.numeric(as.vector(MC.modelALL$t)) #assign numeric t values from OPLS results
MC.t.ortho.opls <- as.numeric(as.vector(MC.modelALL$t_ortho)) #assign numeric t-orthogonal values from OPLS results
MC.metabo.data <- data.frame(cleaned.MC.Treat, MC.t.opls, MC.t.ortho.opls) #create a dataframe
colnames(MC.metabo.data) <- c("Treatment", "t", "t.ortho") #Rename columns

Fig14 <- ggplot(MC.metabo.data, aes(t, t.ortho)) + 
  geom_point(aes(x=t, y=t.ortho, colour=Treatment)) +
  scale_colour_manual("Treatment", values = c("blue","red")) + #Set colors for Species
  xlim(-1600, 1600) +
  ylim(-2500, 2500) +
  ggtitle("B  Montipora capitata") + #Label graph
  xlab("t") + #Label the X Axis
  ylab("t-Orthogonal") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        legend.position="top", #Set position of legend in graphic
        plot.title=element_text(hjust=0, face="bold.italic")) #Justify the title to the top left and use Bold italics font
Fig14  #View figure

#POCILLIPORA OUTLIER CHECK
PD.data <-subset(cleaned.all.data, cleaned.all.info$Species=='Pocillopora damicornis') #Select only the Pocillopora data from the dataset from which outliers have been removed.
PD.info<-subset(cleaned.all.info, cleaned.all.info$Species=='Pocillopora damicornis') #Select only the Pocillopora info
which(colSums(PD.data) ==0) #check to see if any bins (columns) are present in the whole analysis, but not in Pocillopora only
PD.data <- PD.data [, colSums(PD.data  != 0) > 0] #remove the coloumns that have zeros
which(colSums(PD.data) ==0) #check again for any zero columns

#assign treatment from info
PD.Treatment <- PD.info$Treatment

#Scale the data
PD.NMRData.scale <- scaling(PD.data, type="pareto")

#Mean center the data
PD.NMRData.scale <- scale(PD.NMRData.scale, center=TRUE, scale=FALSE)

#calculate similarity matrix of euclidean distance
PD.dist.euc <- vegdist(PD.NMRData.scale, method = "euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) #calculates euclidean distance on the transformed data matrix
PD.pca.res <- prcomp(PD.dist.euc, retx=TRUE)
summary(PD.pca.res)

biplot(PD.pca.res)
ordiellipse(PD.pca.res,group=rep(c(1), c(length(PD.Treatment))), conf=0.99, kind='sd', draw='lines')  

#remove outliers and assess sample size
cleaned.PD.data <- PD.NMRData.scale[-c(17),] #remove outliers identified in PCA
cleaned.PD.info <- PD.info[-c(17),] #remove sample info of outliers identified in PCA
cleaned.PD.Treat <- cleaned.PD.info$Treatment #Set the treatment information for the filtered data
PD.counts <- aggregate(cleaned.PD.info["Coral.ID"], by=cleaned.PD.info[c("Species","Treatment")], FUN=length) #check the sample size for each Treatment
PD.counts #view the sample size for each Treatment

#calculate similarity matrix of euclidean distance
PD.dist.euc <- vegdist(cleaned.PD.data, method = "euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = TRUE) #calculates euclidean distance on the transformed data matrix
PD.pca.res <- prcomp(PD.dist.euc, retx=TRUE)
summary(PD.pca.res)

#Double check for outliers
biplot(PD.pca.res)
ordiellipse(PD.pca.res,group=rep(c(1), c(length(cleaned.PD.Treat))), conf=0.99, kind='sd', draw='lines')  

# extract all PCs
PD.scores <- as.data.frame(PD.pca.res$x)
PD.PCA.plot <- ggplot(PD.scores, aes(x = PC1, y = PC2, shape=cleaned.PD.Treat, col=cleaned.PD.Treat)) +
  geom_point(size=4) + # use points + 
  scale_shape_manual(values=c(17,17)) + #create your own scale and values for scale
  theme(panel.background = element_rect(fill='white', colour='black')) + #removes gray background
  theme(panel.grid.major = element_line(colour="white"), panel.grid.minor = element_line(colour="white")) + #Set the panel gridlines
  theme(legend.key = element_rect(fill='white')) + #Set plot legend key
  theme(panel.border = element_rect(fill=NA, size = 1, colour = "black")) + #Set the border
  ggtitle("P. damicornis")
PD.PCA.plot

PD.X<-as.matrix(cleaned.PD.data) #SELECT DATA -- USER INPUT
class(PD.X)<-'numeric' #changes class data to numeric 
PD.Y<-as.matrix(as.numeric(as.factor(as.vector(cleaned.PD.info$Treatment)))) #SELECT GROUPING -- USER INPUT BY CHANGING COLUMN USED
PD.ulabels <- unique(as.vector(cleaned.PD.info$Treatment)) #PICK LABELS 
set.seed(10)
PD.resultsALL<-n.group.opls(PD.X,PD.Y, num_permutations=100, CV= 10, nIterations=100, min_num_OPLS_fact=1) #RUNS MODEL
PD.Q2<-PD.resultsALL$Q2 #GIVES Q2 VALUE
PD.pval<-PD.resultsALL$helper.results$pvalue #GIVES P-VALUE FOR MODEL 
PD.modelALL <- PD.resultsALL$helper.results$model
PD.OPLSResults<-data.frame(model=c('P. damicornis'), Q2=c(PD.Q2), pval=c(PD.pval))
PD.OPLSResults #View Results

#visualize
PD.t.opls <- as.numeric(as.vector(PD.modelALL$t)) #assign numeric t values from OPLS results
PD.t.ortho.opls <- as.numeric(as.vector(PD.modelALL$t_ortho)) #assign numeric t-orthogonal values from OPLS results
PD.metabo.data <- data.frame(cleaned.PD.Treat, PD.t.opls, PD.t.ortho.opls) #create a dataframe
colnames(PD.metabo.data) <- c("Treatment", "t", "t.ortho") #Rename columns

Fig15 <- ggplot(PD.metabo.data, aes(t, t.ortho)) + 
  geom_point(aes(x=t, y=t.ortho, colour=Treatment)) +
  scale_colour_manual("Treatment", values = c("blue","red")) + #Set colors for Species
  xlim(-1600, 1600) +
  ylim(-2500, 2500) +
  ggtitle("C  Pocillopora damicornis") + #Label graph
  xlab("t") + #Label the X Axis
  ylab("t-Orthogonal") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        legend.position="top", #Set position of legend in graphic
        plot.title=element_text(hjust=0, face="bold.italic")) #Justify the title to the top left and use Bold italics font
Fig15  #View figure

#------------------------------------------------
#Growth Analysis

#load weight data
weight <- read.csv("BM_Buoyant_Weight.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

PD.Arag <- 2.78 #g cm-3, set aragonite density for Pocillopora from literature
MC.Arag <- 2.03 #g cm-3, set aragonite density for Montipora from literature
#Pocillopora Arag = 2.78 g cm-3 Spencer-Davies 1989 Mar Bio 101:389-395, Al-Sofyani and Floos 2013 Oecologia 55:917-935 
#Montipora Arag = 2.03 g cm-3 average from table 2 in Anthony 2003 Functional Ecology 17:246-259

#TimeInitial date=20140501
Temp0 <- 23.92 #Temperature of the measurement
RefDry0 <- 39.110 #weight of solid glass dry reference in air
RefSW0 <- 21.129 #weight of solid glass reference in seawater
RefFW0 <- 21.582 #weight of solid glass reference in freshwater
DenFW0 <- ((-0.000005*Temp0*Temp0)+(0.000007*Temp0)+1.0001) #g cm-3 Calculates density of the freshwater as a function of temperature
DenRef0 <- ((RefDry0*DenFW0)/(RefDry0-RefFW0)) #g cm-3 Equation 4 from Spencer Davies 1989 calculates the density of the reference weight
DenSW0 <- ((DenRef0*(RefDry0-RefSW0))/(RefDry0))  #g cm-3 Equation 3 from Spencer Davies 1989 using the density of the calcualted reference weight

#TimeFinal date=20140612
Temp5 <- 26.24 #Temperature of the measurement
RefDry5 <- 39.110 #weight of solid glass dry reference in air
RefSW5 <- 21.128 #weight of solid glass reference in seawater
RefFW5 <- 21.572 #weight of solid glass reference in freshwater
DenFW5 <- ((-0.000005*Temp5*Temp5)+(0.000007*Temp5)+1.0001) #g cm-3 Calculates density of the freshwater as a function of temperature
DenRef5 <- ((RefDry5*DenFW5)/(RefDry5-RefFW5)) #g cm-3 Equation 4 from Spencer Davies 1989 calculates the density of the reference weight
DenSW5 <- ((DenRef5*(RefDry5-RefSW5))/(RefDry5))  #g cm-3 Equation 3 from Spencer Davies 1989 using the density of the calcualted reference weight

days <- 43

#Montipora Dry Weight
#Subset weights by species
MC.weight <-subset(weight, weight$Species=='Montipora capitata') 
PD.weight <-subset(weight, weight$Species=='Pocillopora damicornis') 

MC.dryweight.initial <- ((MC.weight$Initial_01May)/(1-(DenSW0/MC.Arag))) #Calculates the dry weight in grams using Equation 1 from Spencer Davies 1989 calculates the density of the reference weight
MC.dryweight.final <- ((MC.weight$Initial_12June)/(1-(DenSW5/MC.Arag))) #Calculates the dry weight in grams using Equation 1 from Spencer Davies 1989 calculates the density of the reference weight

PD.dryweight.initial <- ((PD.weight$Initial_01May)/(1-(DenSW0/PD.Arag))) #Calculates the dry weight in grams using Equation 1 from Spencer Davies 1989 calculates the density of the reference weight
PD.dryweight.final <- ((PD.weight$Initial_12June)/(1-(DenSW5/PD.Arag))) #Calculates the dry weight in grams using Equation 1 from Spencer Davies 1989 calculates the density of the reference weight

MC.growth.rate <- ((MC.dryweight.final-MC.dryweight.initial)/MC.dryweight.initial)*100/days #calculates growth in % per day (grams changed/grams initial/day)
PD.growth.rate <- ((PD.dryweight.final-PD.dryweight.initial)/PD.dryweight.initial)*100/days #calculates growth in % per day (grams changed/grams initial/day)


MC.growth <- data.frame(MC.weight,MC.growth.rate)
colnames(MC.growth) <- c("Species",	"Coral.ID",	"Treatment",	"Initial.weight",	"Final.weight",	"Growth.Rate")

PD.growth <- data.frame(PD.weight,PD.growth.rate)
colnames(PD.growth) <- c("Species",  "Coral.ID",	"Treatment",	"Initial.weight",	"Final.weight",	"Growth.Rate")

growth.rate <-rbind(MC.growth,PD.growth)

growth.counts <- aggregate(growth.rate["Growth.Rate"], by=growth.rate[c("Species","Treatment")], FUN=length)
growth.means <- aggregate(Growth.Rate ~ Species + Treatment, data=growth.rate, mean, na.rm = TRUE)
growth.se <- aggregate(Growth.Rate ~ Species + Treatment, data=growth.rate, std.error, na.rm = TRUE)
growth.means <- cbind(growth.means, growth.se$Growth.Rate) 
colnames(growth.means) <- c("Species", "Treatment", "Growth.mean", "Growth.se")
growth.means[, "TS"] <- c("MC Amb", "PD Amb", "MC High", "PD High")
growth.means

Fig16 <- ggplot(growth.means, aes(x=Species, y=Growth.mean, fill=Treatment)) + 
  geom_bar() +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3, # set line weight
           show_guide=FALSE) + #remove slash from legend
  scale_fill_manual(values=c("blue", "red")) +  #Set bar color
  geom_errorbar(aes(ymin=Growth.mean-Growth.se, ymax=Growth.mean+Growth.se),
                size=.3,    # set line weight
                width=.2,   # set error bar width
                position=position_dodge(.9)) + # offset the lines on the x axis so they match the bars and don't overlap
  xlab("Species") + #Label the X Axis
  ylab("Growth (% d-1)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig16 #View figure


###Testing ANOVA Assumptions
###Homogeneity of variance
#get unstandardized predicted and residual values
growth.results <-aov(log10(Growth.Rate) ~ Treatment + Species + Treatment*Species, data=growth.rate)
growth.stats <- anova(growth.results)
unstandardizedPredicted <- predict(growth.results)
unstandardizedResiduals <- resid(growth.results)

#get standardized values and plot to determine normality
standardizedPredicted <- (unstandardizedPredicted - mean(unstandardizedPredicted)) / sd(unstandardizedPredicted)
standardizedResiduals <- (unstandardizedResiduals - mean(unstandardizedResiduals)) / sd(unstandardizedResiduals)
#create standardized residuals plot
plot(standardizedPredicted, standardizedResiduals, main = "Standardized Residuals Plot", xlab = "Standardized Predicted Values", ylab = "Standardized Residuals")
#add horizontal line
abline(0,0)
leveneTest(log10(Growth.Rate) ~Treatment * Species, data=growth.rate)

###Data are normally distributed
results.stdres <- rstandard(growth.results)
qqnorm(results.stdres) # normal quantile plot
qqline(results.stdres) # adding a qline of comparison
shapiro.test(results.stdres) #runs a normality test using shapiro-wilk test on the standardized residuals

growth.HSD <- TukeyHSD(growth.results, ordered = FALSE, conf.level = 0.95)
growth.HSD

#calculate backtransformed means and SE
growth.rate$rate.trans <- log10(growth.rate$Growth.Rate)
sp.growth.trans.avg <- aggregate(rate.trans ~ Species, data=growth.rate, mean, na.rm = TRUE) #calculate the averages by Species 
sp.growthbacktrans.avg <- 10^(sp.growth.trans.avg$rate.trans) #backtransform the average

growth.trans.avg <- aggregate(rate.trans  ~ Species + Treatment, data=growth.rate, mean, na.rm = TRUE) #calculate the averages by Species and Treatment
growth.trans.se <- aggregate(rate.trans  ~ Species + Treatment, data=growth.rate, std.error, na.rm = TRUE) #calculate the standard errors by Species and Treatment
growth.trans.std <- aggregate(rate.trans  ~ Species + Treatment, data=growth.rate, sd, na.rm = TRUE) #calculate the standard deviation by Species and Treatment
growth.trans <- growth.trans.avg #generate dataframe of averages
colnames(growth.trans) <- c("Species", "Treatment", "avg")
growth.trans$Upper <- growth.trans$avg+growth.trans.se$rate.trans
growth.trans$Lower <- growth.trans$avg-growth.trans.se$rate.trans
back.trans.growth <- growth.trans
back.trans.growth$avg <- 10^(growth.trans$avg)
back.trans.growth$Upper <- 10^(growth.trans$Upper)
back.trans.growth$Lower <- 10^(growth.trans$Lower)
back.trans.growth

#Plot backtransformed means and SE
Fig17 <- ggplot(data=back.trans.growth, aes(x=factor(Species), y=(avg), fill=Treatment)) + #plot methylation data averages
  geom_bar()+
  geom_bar(stat="identity", position=position_dodge(), #use values and offset position so bars don't overlap
           colour="black", # set line color,
           size=.3, #set line weight
           show_guide=FALSE) + #remove slash from legend
  scale_fill_manual(values=c("blue", "red")) + #set fill color for bars
  geom_errorbar(aes(ymin=Lower, ymax=Upper), #plot standard errors
                width=.2,                    # set width of the error bars
                position=position_dodge(.9)) + #offset position so bars don't overlap
  ggtitle("A") + #set plot title
  xlab("Species") + #Label the X Axis
  ylab("Growth (% d-1)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig17 #View figure

#------------------------------------------------
#DNA METHYLATION ANALYSIS

#load DNA Methylation data
data <- read.csv("BM_Methylation.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
data <- na.omit(data) #omit data rows that include NA
counts <- aggregate(data["Methylation"], by=data[c("Species","Treatment")], FUN=length) #count the sample sizes
avg <- aggregate(Methylation ~ Species + Treatment, data=data, mean, na.rm = TRUE) #calculate the averages by Species and Treatment
se <- aggregate(Methylation ~ Species + Treatment, data=data, std.error, na.rm = TRUE) #calculate the standard errors by Species and Treatment
methyl <- cbind(avg, se$Methylation) #generate dataframe of averages and standard errors
colnames(methyl) <- c("Species", "Treatment", "avg", "se" ) #rename columns of dataframe

Fig18 <- ggplot(data=methyl, aes(x=factor(Species), y=(avg), fill=Treatment)) + #plot methylation data averages
  geom_bar()+
  geom_bar(stat="identity", position=position_dodge(), #use values and offset position so bars don't overlap
           colour="black", # set line color,
           size=.3, #set line weight
           show_guide=FALSE) + #remove slash from legend
  scale_fill_manual(values=c("blue", "red")) + #set fill color for bars
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), #plot standard errors
                width=.2,                    # set width of the error bars
                position=position_dodge(.9)) + #offset position so bars don't overlap
  ggtitle("B") + #set plot title
  xlab("Species") + #Label the X Axis
  ylab("% DNA Methylation")+ #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig18 #View figure

###Testing ANOVA Assumptions
methylation.results <- aov((Methylation^0.25) ~Treatment * Species, data=data)
methylation.stats <- anova(methylation.results)

###Testing ANOVA Assumptions
###Homogeneity of variance
#get unstandardized predicted and residual values
unstandardizedPredicted <- predict(methylation.results)
unstandardizedResiduals <- resid(methylation.results)

#get standardized values and plot to determine homogeneity
standardizedPredicted <- (unstandardizedPredicted - mean(unstandardizedPredicted)) / sd(unstandardizedPredicted)
standardizedResiduals <- (unstandardizedResiduals - mean(unstandardizedResiduals)) / sd(unstandardizedResiduals)
#create standardized residuals plot
plot(standardizedPredicted, standardizedResiduals, main = "Standardized Residuals Plot", xlab = "Standardized Predicted Values", ylab = "Standardized Residuals")
#add horizontal line
abline(0,0)
leveneTest((Methylation^0.33) ~Treatment * Species, data=data)

###Data are normally distributed
results.stdres <- rstandard(methylation.results)
qqnorm(results.stdres) # normal quantile plot
qqline(results.stdres) # adding a qline of comparison
shapiro.test(results.stdres) #runs a normality test using shapiro-wilk test on the standardized residuals

methylation.HSD <- TukeyHSD(methylation.results, ordered = FALSE, conf.level = 0.95) #calculate pairwise differences between groups using Tukey's HSD Test
methylation.HSD #View results

MCPDavg <- aggregate(Methylation ~ Species, data=data, mean, na.rm = TRUE) #Calculate methylation average for each species
MCPDavg #View results

#calculate backtransformed means and SE
data$Methylation.trans <- (data$Methylation)^0.25
sp.trans.avg <- aggregate(Methylation.trans ~ Species, data=data, mean, na.rm = TRUE) #calculate the averages by Species 
sp.backtrans.avg <- (sp.trans.avg$Methylation.trans)^4 #backtransform the average

trans.avg <- aggregate(Methylation.trans ~ Species + Treatment, data=data, mean, na.rm = TRUE) #calculate the averages by Species and Treatment
trans.se <- aggregate(Methylation.trans ~ Species + Treatment, data=data, std.error, na.rm = TRUE) #calculate the standard errors by Species and Treatment
trans.std <- aggregate(Methylation.trans ~ Species + Treatment, data=data, sd, na.rm = TRUE) #calculate the standard deviation by Species and Treatment
trans.methyl <- trans.avg #generate dataframe of averages
colnames(trans.methyl) <- c("Species", "Treatment", "avg")
trans.methyl$Upper <- trans.methyl$avg+trans.se$Methylation.trans
trans.methyl$Lower <- trans.methyl$avg-trans.se$Methylation.trans
back.trans.methyl <- trans.methyl
back.trans.methyl$avg <- (trans.methyl$avg)^4
back.trans.methyl$Upper <- (trans.methyl$Upper)^4
back.trans.methyl$Lower <- (trans.methyl$Lower)^4
back.trans.methyl

#Plot backtransformed means and SE
Fig19 <- ggplot(data=back.trans.methyl, aes(x=factor(Species), y=(avg), fill=Treatment)) + #plot methylation data averages
  geom_bar()+
  geom_bar(stat="identity", position=position_dodge(), #use values and offset position so bars don't overlap
           colour="black", # set line color,
           size=.3, #set line weight
           show_guide=FALSE) + #remove slash from legend
  scale_fill_manual(values=c("blue", "red")) + #set fill color for bars
  geom_errorbar(aes(ymin=Lower, ymax=Upper), #plot standard errors
                width=.2,                    # set width of the error bars
                position=position_dodge(.9)) + #offset position so bars don't overlap
  ggtitle("B") + #set plot title
  xlab("Species") + #Label the X Axis
  ylab("% DNA Methylation")+ #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig19 #View figure


#------------------------------------------------
#CAPTURE ALL STATISTICAL OUTPUT TO A FILE

capture.output(MC.OPLSResults, PD.OPLSResults, ALL.OPLSResults, growth.results, growth.stats, growth.HSD, methylation.results, methylation.stats, methylation.HSD, file="Bulk_Methylation_Statistical_Results.txt")

Figure.S1 <- arrangeGrob(Fig3, Fig2, Fig1, ncol=3)
ggsave(file="FigureS1_Manuscript.pdf", Figure.S1, width = 11, height = 4, units = c("in"))

Figure.S2 <- arrangeGrob(Fig6, Fig9, ncol=2)
ggsave(file="FigureS2_Manuscript.pdf", Figure.S2, width = 11, height = 4, units = c("in"))

Figure.1 <- arrangeGrob(Fig11, Fig12, ncol=2)
ggsave(file="Figure1_Manuscript.pdf", Figure.1, width = 11, height = 4, units = c("in"))

Figure.2 <- arrangeGrob(Fig13, Fig14, Fig15, ncol=1)
ggsave(file="Figure2_Manuscript.pdf", Figure.2, width = 5, height = 11, units = c("in"))

Figure.3 <- arrangeGrob(Fig17, Fig19, ncol=2)
ggsave(file="Figure3_Manuscript.pdf", Figure.3, width = 11, height = 4, units = c("in"))


