library(readr)
data <- read_csv("Descargas/Ali_amb_temp_sard.vs.eez.csv")
View(data)

# cantidad de observaciones por cada variable
length(data$temp.media) # 257
length(data$temp.max) # 257
length(data$temp.min) # 257

# creando el df con las tres variables que nos interesan
data_sardinas <- data.frame(data[,3],data[,4],data[,5])
data_sardinas
# tapply(data_sardinas$temp.min, data_sardinas$temp.media, data_sardinas$temp.max, mean)
