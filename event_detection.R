library(tidyverse)
library(stringr)
library(imputeTS)
library(plyr)
library(lubridate)
library(Rcpp)

impute_data <- function(x){
  salinity <- read.csv(sprintf("%s/Salinity_Hourly.csv", x))
  i_skip = which(salinity$X...BEGIN.HEADER == "#Time UTC (yyyy-mm-ddThh:mm:ss.fffZ)")
  salinity <- read.csv(sprintf("%s/Salinity_Hourly.csv", x), skip = i_skip, col.names = c("Time", "PSU", "QC Flag", "Count"))[-1,]
  salinity$Time <- as.POSIXct(salinity$Time, format = "%Y-%m-%dT%H:%M:%S", tz = "GMT")
  salinity <- filter(salinity, hour(Time) %in% c(00, 12))
  salinity <- salinity[!(month(salinity$Time) == 2 & day(salinity$Time) == 29),]
  
  temperature <- read.csv(sprintf("%s/Temperature_Hourly.csv", x))
  i_skip = which(temperature$X...BEGIN.HEADER == "#Time UTC (yyyy-mm-ddThh:mm:ss.fffZ)")
  temperature <- read.csv(sprintf("%s/Temperature_Hourly.csv", x), skip = i_skip, col.names = c("Time", "Celsius", "QC Flag", "Count"))[-1,]
  temperature$Time <- as.POSIXct(temperature$Time, format = "%Y-%m-%dT%H:%M:%S", tz = "GMT")
  temperature <- filter(temperature, hour(temperature$Time) %in% c(00, 12))
  temperature <- temperature[!(month(temperature$Time) == 2 & day(temperature$Time) == 29),]
  
  oxygen <- read.csv(sprintf("%s/Oxygen_Hourly.csv", x))
  i_skip = which(oxygen$X...BEGIN.HEADER == "#Time UTC (yyyy-mm-ddThh:mm:ss.fffZ)")
  oxygen <- read.csv(sprintf("%s/Oxygen_Hourly.csv", x), skip = i_skip, col.names = c("Time", "ML_L", "QC Flag", "Count"))[-1,]
  oxygen$Time <- as.POSIXct(oxygen$Time, format = "%Y-%m-%dT%H:%M:%S", tz = "GMT")
  oxygen <- filter(oxygen, hour(oxygen$Time) %in% c(00, 12))
  oxygen <- oxygen[!(month(oxygen$Time) == 2 & day(oxygen$Time) == 29),]
  
  
  impute_salinity <- na_interpolation(salinity$PSU)
  impute_oxygen <- na_interpolation(oxygen$ML_L)
  impute_temp <- na_interpolation(temperature$Celsius)
  
  out <- data.frame(Celsius = impute_temp, ML_L = impute_oxygen, PSU = impute_salinity, Time = salinity$Time)
  outm <- matrix(c(out$Celsius, out$ML_L, out$PSU), ncol = 3)
  
  return(outm)
}

find_cprob <- function(x){
  bchange <- bcp(x)
  #plot(bchange)
  
  return(bchange$posterior.prob)
}

yrs <- 2009:2019
len <- length(yrs)
output <- matrix(ncol=len, nrow = 728)

for(i in 1:len){
  temp <- impute_data(yrs[i]) %>% find_cprob()
  output[,i] <- temp
}

start_date = as.POSIXct('01-01 00:30:00', format = "%m-%d %H:%M:%S")
end_date = as.POSIXct('12-30 23:30:00', format =  "%m-%d %H:%M:%S")

#returns the posterior probabilities of an event (change-point) at each timestep
Time <- seq.POSIXt(from = start_date, to = end_date, by = 'hours')
Time <- Time[hour(Time) %in% c(00, 12)]
output <- data.frame(output, Time)
names(output) <- c(yrs, "Time")
output <- na.omit(output)

write.csv(output, "dwr_bcp.csv")

#trains best arima model
library(forecast)
fit <- auto.arima(c(output$`2009`, output$`2010`, output$`2011`))
prediction <- forecast(fit, h = 365)
plot(prediction)

fit2 <- arima(c(output$`2009`, output$`2010`, output$`2011`), order = c(2,0,0))
prediction2 <- forecast(fit2, h=365)
plot(prediction2)

library(sarima)
library(astsa)
fit3 <- sarima.for(c(output$`2009`, output$`2010`, output$`2011`), n.ahead = 365, 4,0,0,1,0,0,12)
prediction3 <- forecast(fit3, h=365)


