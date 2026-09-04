library(readr)
# data <- read_csv("Descargas/Ali_amb_temp_sard.vs.eez.csv")
View(data)

data <- read.csv("C:/Users/sebas/Desktop/TG-Seminario/adaptive-density-weak-dep/scripts/real_data_study/Ali_amb_temp_sard.vs.eez.csv")

# cantidad de observaciones por cada variable
length(data$temp.media) # 257
length(data$temp.max) # 257
length(data$temp.min) # 257

# creando el df con las tres variables que nos interesan
data_sardinas <- data.frame(data[,12],data[,13],data[,14])
data_sardinas
colnames(data_sardinas) <- c("tmin.sard", "tmax.sard", "tmean.sard")
par(mfrow = c(1, 3))
acf(data_sardinas$tmin.sard); acf(data_sardinas$tmax.sard); acf(data_sardinas$tmean.sard)
ndiffs(data_sardinas$tmin.sard) # 0
ndiffs(data_sardinas$tmax.sard) # 1
ndiffs(data_sardinas$tmean.sard) # 0

par(mfrow = c(1, 3))
acf(data_sardinas$tmin.sard); acf(diff(data_sardinas$tmax.sard, diff =1)); acf(diff(data_sardinas$tmean.sard, diff=1))

mmodelo1 <- auto.arima(x = data_sardinas$tmin.sard)
mmodelo1 # ARIMA(2,1,2) with non-zero mean

mmodelo2 <- auto.arima(x = data_sardinas$tmax.sard)
mmodelo2 # ARIMA(0,1,2) 

mmodelo3 <- auto.arima(x = data_sardinas$tmax.sard)
mmodelo3 # ARIMA(0,1,2) 

decompose(data_sardinas$tmin.sard)
# tapply(data_sardinas$temp.min, data_sardinas$temp.media, data_sardinas$temp.max, mean)
