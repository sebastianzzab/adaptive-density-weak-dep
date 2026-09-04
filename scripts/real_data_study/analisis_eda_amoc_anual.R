# ============================================================
# Analisis exploratorio de la serie AMOC Proxy Subpolar Gyre
# HadISST
#
# La serie original es mensual. Este script:
# 1. Lee los datos mensuales.
# 2. Calcula el promedio anual.
# 3. Obtiene estadisticas descriptivas.
# 4. Grafica la serie observada anual.
# 5. Realiza una descomposicion aditiva de tendencia y residuo.
# 6. Culmina con un histograma de la serie anual.
#
# Nota metodologica
# Al agregar los datos mensuales a frecuencia anual se elimina la
# estacionalidad intra-anual. Por ello, la componente estacional
# de la serie anual se representa como cero. La descomposicion
# conserva las componentes observada, tendencia e irregular.
# ============================================================

# -----------------------------
# 1. Configuracion
# -----------------------------
library(readxl)#importar hojas de Excel a R#
library(aTSA) #análisis detallados de series temporales
library(astsa)
library(car)
library(lubridate)
library(vars)
library(tseries)#permitirá manipular datos correspondiente a series temporale#
library(foreign)#sirve para importar datos en diveros formatos#
library(quantmod)#unciones para la descarga, manipulación, visualización y construcción de modelos con datos financieros#
library(forecast)# Contiene el modelo ARIMA
library(tseries) #Para series de tiempo
library(TSA)     #Para series de tiempo
library(urca)    #Para hacer el Test de Raiz Unitaria (detectar hay o no estacionariedad)
library(ggplot2) #Para hacer gráficos
library(dplyr)   #Para la manipulación de datos (filtrar, seleccionar, agregar, transformar)#
library(stats)   #Se usa para diversas pruebas estadísticas (medias,varianza, arima,etc)#
library(lmtest)
library(mFilter)
library(dynlm)
library(nlme)
library(NTS)
library(broom)
library(rugarch)


archivo_entrada <- "AMOC_Proxy_SubpolarGyre_HadISST.csv"
directorio_salida <- "salidas_eda_amoc"

if (!dir.exists(directorio_salida)) {
  dir.create(directorio_salida, recursive = TRUE)
}

# -----------------------------
# 2. Lectura y preparacion
# -----------------------------

datos_mensuales <- read.csv(
  archivo_entrada,
  header = TRUE,
  stringsAsFactors = FALSE,
  fileEncoding = "ASCII"
)

datos_mensuales$date <- as.Date(datos_mensuales$date)
datos_mensuales$year <- as.integer(datos_mensuales$year)
datos_mensuales$amoc_proxy_sst <- as.numeric(datos_mensuales$amoc_proxy_sst)

if (anyNA(datos_mensuales$date) ||
    anyNA(datos_mensuales$year) ||
    anyNA(datos_mensuales$amoc_proxy_sst)) {
  stop("Hay valores faltantes o no validos en las columnas necesarias.")
}

# Promedio de los meses disponibles para cada año
serie_anual <- aggregate(
  amoc_proxy_sst ~ year,
  data = datos_mensuales,
  FUN = mean,
  na.rm = TRUE
)

names(serie_anual)[names(serie_anual) == "amoc_proxy_sst"] <- "promedio_anual"
serie_anual <- serie_anual[order(serie_anual$year), ]
rownames(serie_anual) <- NULL

# write.csv(
#   serie_anual,
#   file.path(directorio_salida, "serie_AMOC_promedio_anual.csv"),
#   row.names = FALSE
# )

# Serie temporal anual regular
serie_ts <- ts(
  serie_anual$promedio_anual,
  start = min(serie_anual$year),
  frequency = 1
)

# -----------------------------
# 3. Estadisticas descriptivas
# -----------------------------

estadisticas_descriptivas <- data.frame(
  estadistica = c(
    "Numero de observaciones",
    "Anio inicial",
    "Anio final",
    "Media",
    "Mediana",
    "Desviacion estandar",
    "Varianza",
    "Minimo",
    "Maximo",
    "Rango",
    "Primer cuartil",
    "Tercer cuartil",
    "Rango intercuartil",
    "Coeficiente de variacion_pct"
  ),
  valor = c(
    length(serie_anual$promedio_anual),
    min(serie_anual$year),
    max(serie_anual$year),
    mean(serie_anual$promedio_anual),
    median(serie_anual$promedio_anual),
    sd(serie_anual$promedio_anual),
    var(serie_anual$promedio_anual),
    min(serie_anual$promedio_anual),
    max(serie_anual$promedio_anual),
    diff(range(serie_anual$promedio_anual)),
    unname(quantile(serie_anual$promedio_anual, 0.25)),
    unname(quantile(serie_anual$promedio_anual, 0.75)),
    IQR(serie_anual$promedio_anual),
    100 * sd(serie_anual$promedio_anual) /
      mean(serie_anual$promedio_anual)
  )
)

print(estadisticas_descriptivas)

# write.csv(
#   estadisticas_descriptivas,
#   file.path(directorio_salida, "estadisticas_descriptivas_AMOC_anual.csv"),
#   row.names = FALSE
# )

# -----------------------------
# 4. Grafico de la serie observada
# -----------------------------

# png(
#   filename = file.path(directorio_salida, "01_serie_observada_anual.png"),
#   width = 1200,
#   height = 700,
#   res = 120
# )
# plot(
#   serie_ts,
#   type = "o",
#   pch = 16,
#   cex = 0.55,
#   col = "steelblue4",
#   # main = "Serie observada - promedio anual",
#   xlab = "Año",
#   ylab = "AMOC proxy SST"
# )
# grid()
# dev.off()

plot(serie_ts, type = "o", pch = 16, cex = 0.55,
  col = "steelblue4",# main = "Serie observada - promedio anual",
  xlab = "Año", ylab = "AMOC proxy SST")
acf(serie_ts)
pacf(serie_ts)
ndiffs(x = serie_ts) # 1
par(mfrow = c(1, 2))
plot(serie_ts, type = "o", pch = 16, cex = 0.55,
     col = "steelblue4",# main = "Serie observada - promedio anual",
     xlab = "Año", ylab = "AMOC proxy SST")
plot(diff(serie_ts, 1), type = "o", pch = 16, cex = 0.55, 
     col = "steelblue4", xlab = "Año", ylab = "AMOC proxy SST diferenciada")
par(mfrow = c(1, 1))
# decompose(serie_ts)

mmodelo <- auto.arima(x = serie_ts)
mmodelo # ARIMA(2,1,2)

# Realizar prueba de raíz unitaria#
#H0: Serie No estacionaria: Hay raiz unitaria H1: Serie Estacionaria: No hay raiz unitaria#
adf.test(serie_ts) # p-value = 0.3228
adf.test(diff(serie_ts,1)) # p-value = 0.01
#Se aplica la regla del P-valor y se concluye que la serie es estacionaria en media#
diff_serie_ts <- diff(serie_ts, diff=1)

acf(diff_serie_ts, lag.max = 60)
pacf(diff_serie_ts, lag.max = 60)

# -----------------------------
# 5. Descomposicion aditiva
# -----------------------------

# Como la serie es anual, no existe una frecuencia intra-anual
# que permita estimar una estacionalidad mensual. Se usa una
# tendencia suavizada de 11 años y se calcula el residuo.
ventana_tendencia <- 11

tendencia <- as.numeric(
  stats::filter(
    serie_anual$promedio_anual,
    filter = rep(1 / ventana_tendencia, ventana_tendencia),
    sides = 2
  )
)

componente_estacional <- rep(0, nrow(serie_anual))
componente_irregular <- serie_anual$promedio_anual -
  tendencia -
  componente_estacional

descomposicion <- data.frame(
  year = serie_anual$year,
  observado = serie_anual$promedio_anual,
  tendencia = tendencia,
  estacional = componente_estacional,
  irregular = componente_irregular
)

# write.csv(
#   descomposicion,
#   file.path(directorio_salida, "descomposicion_AMOC_anual.csv"),
#   row.names = FALSE
# )
# 
# png(
#   filename = file.path(directorio_salida, "02_descomposicion_AMOC_anual.png"),
#   width = 1200,
#   height = 1000,
#   res = 120
# )

par(mfrow = c(4, 1), mar = c(3, 4, 2, 1))

plot(
  descomposicion$year,
  descomposicion$observado,
  type = "l",
  col = "black",
  main = "Observada",
  xlab = "",
  ylab = "Valor"
)
#grid()

plot(
  descomposicion$year,
  descomposicion$tendencia,
  type = "l",
  col = "firebrick",
  lwd = 2,
  main = paste0("Tendencia suavizada (ventana de ", ventana_tendencia, " años)"),
  xlab = "",
  ylab = "Valor"
)
#grid()

plot(
  descomposicion$year,
  descomposicion$estacional,
  type = "l",
  col = "darkgreen",
  main = "Componente estacional",
  xlab = "",
  ylab = "Valor"
)
abline(h = 0, lty = 2)
#grid()

plot(
  descomposicion$year,
  descomposicion$irregular,
  type = "h",
  col = "purple4",
  main = "Componente irregular / residuo",
  xlab = "Año",
  ylab = "Residuo"
)
abline(h = 0, lty = 2)
#grid()

par(mfrow = c(1, 1))
#dev.off()

# -----------------------------
# 6. Histograma de la serie anual
# -----------------------------

# png(
#   filename = file.path(directorio_salida, "03_histograma_AMOC_anual.png"),
#   width = 1000,
#   height = 700,
#   res = 120
# )

hist(
  serie_anual$promedio_anual,
  breaks = "FD",
  # col = "skyblue2",
  border = "white",
  # main = "Histograma del promedio anual",
  main = "",
  xlab = "AMOC proxy SST",
  ylab = "Frecuencia"
)
abline(
  v = mean(serie_anual$promedio_anual),
  col = "red3",
  lwd = 2,
  lty = 2
)
legend(
  "topright",
  legend = "Media",
  col = "red3",
  lwd = 2,
  lty = 2,
  bty = "n"
)
#grid()
#dev.off()

# -----------------------------
# 7. Salida resumida en consola
# -----------------------------

cat("\nAnalisis completado.\n")
cat("Archivo de entrada:", archivo_entrada, "\n")
cat("Periodo anual:", min(serie_anual$year), "-", max(serie_anual$year), "\n")
cat("Numero de años:", nrow(serie_anual), "\n")
cat("Resultados guardados en:", directorio_salida, "\n")
