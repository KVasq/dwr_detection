library(tidyverse)
library(stringr)
library(ggplot2)
library(Rcpp)
library(lubridate)
library(plyr)
library(bcp)
#library(mtsdi)
#library(mice)
library(imputeTS)

#Turbidity

turbidity <- read.csv("2009/Turbidity_Hourly.csv", skip = 57)[-1,]
#na.omit(turbidity)

t_ <- data.frame(matrix(NA, nrow = 654, ncol = 4))
names(t_) <- c("Time", "NTU", "QC Flag", "Count")
t_$Time <- seq.POSIXt(as.POSIXct('2009-12-03 18:30:00'), as.POSIXct('2009-12-30 23:30:00'), by = 'hours')
t_$`QC Flag`<- 9
t_$Count <- 0

names(turbidity) <- c("Time", "NTU", "QC Flag", "Count")

turbidity$Time <- as.POSIXct(turbidity$Time, format = "%Y-%m-%dT%H:%M:%S")

turbidity <- rbind(turbidity, t_)




turbidity <- filter(turbidity, hour(Time) %in% c(00, 12))

ggplot(turbidity, aes(x = Time,y = NTU, group = 1)) +
  geom_line()

#Temperature

temperature <- read.csv("2009/Temp_Hourly.csv", skip = 57)[-1,]
#temperature <- na.omit(temperature)

names(temperature) <- c("Time", "Celsius", "QC Flag", "Count")

temperature$Time <- as.POSIXct(temperature$Time, format = "%Y-%m-%dT%H:%M:%S", tz = "GMT")

temperature <- filter(temperature, hour(Time) %in% c(00, 12))

ggplot(temperature, aes(x = Time,y = Celsius, group = 1)) +
  geom_line()

#temp_change <- e.divisive(matrix(temperature$Celsius))

#ggplot(temperature, aes(x = as.numeric(row.names(temperature)), y = Celsius, group = 1)) +
#  geom_line() +
#  geom_vline(xintercept=temp_change$estimates)

#Oxygen

oxygen <- read.csv("2009/Oxygen_Hourly.csv", skip = 57)[-1,]
#oxygen <- na.omit(oxygen)

names(oxygen) <- c("Time", "ML_L", "QC Flag", "Count")

oxygen$Time <- as.POSIXct(oxygen$Time, format = "%Y-%m-%dT%H:%M:%S", tz = "GMT")

oxygen <- filter(oxygen, hour(Time) %in% c(00, 12))

ggplot(oxygen, aes(x = Time,y = ML_L, group = 1)) +
  geom_line()

#oxy_change <- e.divisive(matrix(oxygen$ML_L))

#ggplot(oxygen, aes(x = as.numeric(row.names(oxygen)), y = ML_L, group = 1)) +
#  geom_line() +
#  geom_vline(xintercept=oxy_change$estimates)

#Pressure

pressure <- read.csv("2009/Pressure_Hourly.csv", skip = 57)[-1,]
#na.omit()

names(pressure) <- c("Time", "Decibar", "QC Flag", "Count")

pressure$Time <- as.POSIXct(pressure$Time, format = "%Y-%m-%dT%H:%M:%S", tz = "GMT")

pressure <- filter(pressure, hour(Time) %in% c(00, 12))

ggplot(pressure, aes(x = Time,y = Decibar, group = 1)) +
  geom_line()

#pressure_change <- e.divisive(matrix(pressure$Decibar), min.size = 9)

#ggplot(pressure, aes(x = as.numeric(row.names(pressure)), y = Decibar, group = 1)) +
#  geom_line() +
#  geom_vline(xintercept=pressure_change$estimates)

#Salinity

salinity <- read.csv("2009/Salinity_Hourly.csv", skip = 57)[-1,]
#na.omit()

names(salinity) <- c("Time", "PSU", "QC Flag", "Count")

salinity$Time <- as.POSIXct(salinity$Time, format = "%Y-%m-%dT%H:%M:%S", tz = "GMT")

salinity <- filter(salinity, hour(Time) %in% c(00, 12))

ggplot(salinity, aes(x = Time,y = PSU, group = 1)) +
  geom_line()

#sal_change <- e.divisive(matrix(salinity$PSU))

#ggplot(salinity, aes(x = as.numeric(row.names(salinity)), y = PSU, group = 1)) +
#  geom_line() +
#  geom_vline(xintercept=sal_change$estimates)

#Multivariate Change-Point Estimation

rbind(salinity$PSU, temperature$Celsius, oxygen$ML_L)

change_test <- kcpa(rbind.fill.matrix(mat1,mat2,mat3), 25, -1)

ggplot(salinity, aes(x = as.numeric(row.names(salinity)), y = PSU, group = 1)) +
  geom_line() +
  geom_vline(xintercept=change_test$estimates)

ggplot(temperature, aes(x = as.numeric(row.names(temperature)), y = Celsius, group = 1)) +
  geom_line() +
  geom_vline(xintercept=change_test$estimates)

ggplot(oxygen, aes(x = as.numeric(row.names(oxygen)), y = ML_L, group = 1)) +
  geom_line() +
  geom_vline(xintercept=change_test$estimates)

mat1 <- matrix(salinity$PSU, nrow = 1)
mat2 <- matrix(oxygen$ML_L, nrow = 1)
mat3 <- matrix(temperature$Celsius, nrow = 1)

rbind.fill.matrix(mat1,mat2,mat3)

rbind(mat1,mat2,mat3)

#unvariate bcp test

bchange_test <- bcp(salinity$PSU)

bchange_test <- bcp(oxygen$ML_L)

bchange_test <- bcp(temperature$Celsius)

bchange_test <- bcp(turbidity$NTU)

bchange_test <- bcp(pressure$Decibar)

plot(bchange_test)

bcp_df <- data.frame(temperature$Celsius, oxygen$ML_L, salinity$PSU, turbidity$NTU, pressure$Decibar)
names(bcp_df) <- c("Celsius", "ML_L", "PSU", "NTU", "Decibar")

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

bcp_df[is.nan(bcp_df)] <- NA


prep_bcp_df <- edaprep(bcp_df)

#multivariate time series imputation test
#mtsdi

f <- ~bcp_df$ML_L + bcp_df$PSU + bcp_df$NTU + bcp_df$Celsius + bcp_df$Decibar
i <- mnimput(f, bcp_df, method = "arima", ar.control = list(order=cbind(c(2,0,1), c(2,0,1), c(2,0,1), c(2,0,1), c(2,0,1))))
out <- predict(i)
names(out) <- c("ML_L", "PSU", "NTU", "Celsius", "Decibar")
out <- mutate(out, Time = temperature$Time)
ggplot(out, aes(x = Time, y = Decibar, group = 1)) +
  geom_line()


#imputeTS
impute_temp <- na_interpolation(temperature$Celsius, option = "spline")
ggplot_na_distribution(impute_temp)
out$Celsius <- impute_temp
temperature$Celsius <- impute_temp

ggplot(out, aes(x = Time, y = impute_temp, group = 1)) +
  geom_line()


impute_pressure <- na_interpolation(pressure$Decibar)
out$Decibar <- impute_pressure

pressure$Decibar <- pressure[!pressure %in% boxplot.stats(out$Decibar)$out,]$Decibar

ggplot(out, aes(x = Time, y = NTU, group = 1)) +
  geom_line()

impute_salinity <- na_interpolation(salinity$PSU)
impute_oxygen <- na_interpolation(oxygen$ML_L)
impute_turbidity <- na_interpolation(turbidity$NTU)

out <- data.frame(impute_oxygen, impute_salinity, impute_turbidity, impute_temp, impute_pressure, pressure$Time)
names(out) <- c("ML_L", "PSU", "NTU", "Celsius", "Decibar", "Time")

ggplot(out, aes(x = Time, y = PSU, group = 1)) +
  geom_line()

#mice
mice(bcp_df, m = 1)


turbidity$NTU <- out$`bcp_df$NTU`
temperature$Celsius <- out$`bcp_df$Celsius`
salinity$PSU <- out$`bcp_df$PSU`
oxygen$ML_L <- out$`bcp_df$ML_L`
pressure$Decibar <- out$`bcp_df$Decibar`



#multivariate bcp test

bchange_test <- bcp(out$PSU)

bchange_test <- bcp(oxygen$ML_L)

bchange_test <- bcp(temperature$Celsius)

bchange_test <- bcp(turbidity$NTU)

bchange_test <- bcp(pressure$Decibar)

plot(bchange_test)

outm <- matrix(c(out$Celsius, out$ML_L, out$PSU), ncol = 3)

bcp(out$PSU)
multi_bchange <- bcp(outm)
plot(multi_bchange)


#correlation matrix
library(Hmisc)

cor(out[1:4])

outm <- matrix(c(out$NTU, out$ML_L, out$PSU, out$Celsius), ncol = 4)
colnames(outm) <- c('turbidity', 'oxygen', 'salinity', 'temperature')
cormatrix <- cor(outm)

library(corrplot)
col1 <- colorRampPalette(c('#7F0000', 'red', '#FF7F00', '#ffbfbf', 'white',
                           'cyan', '#007FFF', 'blue', '#00007F'))
corrplot(cormatrix, col = col1(50))

pal <- colorRampPalette(c('green', 'white', 'red')) (20)
heatmap(cormatrix, col = pal, symm = TRUE)

#Prophet test

library(prophet)

outp <- out[c(2,6)]

names(outp) <- c('y', 'ds')

p <- prophet(outp)

future <- make_future_dataframe(p, period = 7)

plot(p, predict(p, future)) + add_changepoints_to_plot(p)
