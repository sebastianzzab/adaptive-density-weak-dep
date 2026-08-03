# Instalar el paquete si no está disponible
# install.packages("nor1mix")
library(nor1mix)
# Si se desea usar computación en paralelo para acelerar uniroot en simulaciones
library(parallel) 

# =====================================================================
# 1. Definición de Parámetros y Proceso AR(1)
# =====================================================================
n <- 2200
phi <- 0.75

# Generación del AR(1) subyacente
# La varianza teórica de un AR(1) es sigma^2 / (1 - phi^2). Aquí sigma = 1.
Z <- arima.sim(list(ar = phi), n = n)
sigma_Z <- sqrt(1 / (1 - phi^2))

# Transformación integral de probabilidad para obtener dependientes U(0,1)
U <- pnorm(Z, mean = 0, sd = sigma_Z)

# =====================================================================
# 2. Configuración de las Densidades de Marron & Wand
# =====================================================================
# El paquete nor1mix provee objetos predefinidos MW.nm1 hasta MW.nm15.
# Almacenamos estas funciones en una lista para facilitar el acceso por índice.
mw_list <- list(
  MW.nm1,  MW.nm2,  MW.nm3,  MW.nm4,  MW.nm5,
  MW.nm6,  MW.nm7,  MW.nm8,  MW.nm9,  MW.nm10,
  MW.nm11, MW.nm12, MW.nm13, MW.nm14, MW.nm15
)

# =====================================================================
# 3. Función Generadora mediante Transformada Inversa Numérica
# =====================================================================
# Dado que la distribución empírica en nor1mix no está truncada,
# ampliamos el intervalo de búsqueda a valores seguros (ej. -15 a 15).
generar_dependiente_mw <- function(U, mw_obj, limite = c(-15, 15)) {
  
  # Usamos mclapply (o sapply) para optimizar el cálculo de raíces en vectores largos.
  # El método uniroot invierte la función pnorMix (CDF de la mezcla).
  Y_list <- mclapply(U, function(u) {
    res <- uniroot(function(x) pnorMix(x, mw_obj) - u, interval = limite)
    return(res$root)
  }, mc.cores = detectCores() - 1) 
  
  return(unlist(Y_list))
}

# =====================================================================
# 4. Prueba del Esquema (Ejemplo: Densidad #10 - Claw Density)
# =====================================================================
indice_mw <- 3
densidad_objetivo <- mw_list[[indice_mw]]

# Generación de la muestra Y fuertemente estructurada y dependiente
Y_mw <- generar_dependiente_mw(U, densidad_objetivo)

# =====================================================================
# 5. Visualización y Diagnóstico
# =====================================================================
par(mfrow = c(1, 2))

# Diagnóstico de dependencia débil
acf(Y_mw, lag.max = 100, main = paste("ACF - Densidad MW #", indice_mw))

# Comparación de densidades empírica y teórica
xgrid <- seq(-4, 4, length = 500)
# dnorMix calcula la PDF exacta de la mezcla
g_teorica <- dnorMix(xgrid, densidad_objetivo)
ylim_max <- max(g_teorica) * 1.2

hist(Y_mw, prob = TRUE, breaks = 50, col = "lightgray", border = "white",
     ylim = c(0, ylim_max), main = "Densidad Empírica vs Teórica",
     xlab = "Y", ylab = "Densidad")

lines(density(Y_mw), col = "blue", lwd = 2)
lines(xgrid, g_teorica, col = "red", lwd = 2)

legend("topright", legend = c("Estimación por Núcleo", "Teórica (nor1mix)"),
       col = c("blue", "red"), lwd = 2, cex = 0.8)