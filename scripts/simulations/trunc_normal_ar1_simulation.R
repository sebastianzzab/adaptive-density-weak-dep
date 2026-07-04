# ==============================================================================
# TÍTULO: Simulación de Procesos AR(1) con Distribución Normal Truncada y 
#         Comparación de Estimadores de Densidad Kernel (nrd0 vs. ucv)
#
# DESCRIPCIÓN: Este script genera series temporales con dependencia lineal AR(1)
#              transformadas a una distribución Normal Truncada mediante el 
#              método de la transformada inversa. Posteriormente, evalúa el
#              impacto del tamaño muestral (n) y el parámetro de persistencia (phi)
#              en la estimación de la densidad, comparando dos métodos de 
#              ancho de banda (h) frente a la densidad teórica real.
#
# AUTOR: Sebastian Zabala
# FECHA: Julio, 2026
# ==============================================================================

# ==============================================================================
# 1. PAQUETES Y CONFIGURACIÓN
# ==============================================================================
# Instala el paquete si no lo tienes: install.packages("truncnorm")
library(truncnorm)

# ==============================================================================
# 2. FUNCIÓN DE SIMULACIÓN (Normal Truncada Dependiente AR(1))
# ==============================================================================
gnormtrunc_dep <- function(nsim, phi, a, b, m, de) {
  # Generacion de AR(1)
  if (phi == 0) {
    Z <- rnorm(nsim) # Si phi es 0, es simplemente ruido blanco estandar
  } else {
    Z <- arima.sim(list(ar=phi), n=nsim) 
  }
  
  # Calculamos la varianza teórica del proceso AR(1)
  sd_Z <- if (phi == 0) 1 else sqrt(1 / (1 - phi^2))
  
  # Generar variables uniformes
  U <- pnorm(Z, mean=0, sd=sd_Z)
  
  # Soporte para la normal truncada, dependiente.
  F_a <- pnorm(a, mean = m, sd = de)
  F_b <- pnorm(b, mean = m, sd = de)
  p <- F_b - F_a
  
  # Inversa de la función de distribución (qnorm) ajustada a la escala
  X = qnorm(p * U + F_a, mean=m, sd=de)
  
  par(mfrow=c(1,3))
  # Gráficos internos solicitados de la serie
  plot(seq(1, nsim, 1), X, xlab="Índice", ylab="X", main=paste("Dispersión (phi =", phi, ", n =", nsim, ")"))
  plot(X, type = "l", xlab="Índice", ylab="X", main="Gráfico de Líneas")
  acf(X, main=paste("ACF para phi =", phi))
  par(mfrow=c(1,1))
  
  return(X)
}

# ==============================================================================
# 3. FUNCIÓN ASISTENTE PARA SOBREPONER DENSIDADES
# ==============================================================================
graficar_densidades_trunc <- function(X, bw_metodo, a, b, m, de, titulo = "") {
  # 1. Calcular estimación de densidad kernel
  est_dens <- density(X, bw = bw_metodo)
  
  # 2. Extraer el valor numérico de h (ancho de banda) redondeado a 4 decimales
  h_estimado <- round(est_dens$bw, 4)
  
  # 3. Graficar la estimación del Kernel incluyendo el valor de h en el título
  plot(est_dens, 
       main = paste("Densidad", titulo, "[", bw_metodo, "]"), 
       sub = paste("Ancho de banda estimado (h) =", h_estimado),
       col = "blue", lwd = 2, ylim = c(0, max(est_dens$y) * 1.2),
       xlab = "Valores de X", ylab = "Densidad")
  
  # 4. Sobreponer la teórica usando el paquete truncnorm
  curve(dtruncnorm(x, a = a, b = b, mean = m, sd = de), 
        col = "darkred", lwd = 2, lty = 2, add = TRUE)
  
  # 5. Leyenda explicativa enriquecida con el valor de h
  legend("topright", 
         legend = c(paste("Kernel (h =", h_estimado, ")"), "Teórica (truncnorm)"), 
         col = c("blue", "darkred"), lty = c(1, 2), lwd = 2, bty = "n")
}

# ==============================================================================
# 4. NUEVA FUNCIÓN: DIAGNÓSTICO DE ERROR LOCAL (SESGO Y ERROR CUADRÁTICO)
# ==============================================================================
diagnostico_error_local <- function(X, a, b, m, de, titulo = "") {
  # Configurar layout para 2 gráficos apilados verticalmente
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  
  # 1. Definir un grid común para alinear todas las curvas matemáticamente
  # Se extiende un poco más allá de a y b para ver la fuga de densidad (efecto borde)
  x_grid <- seq(a - 0.5, b + 0.5, length.out = 512)
  y_teo <- dtruncnorm(x_grid, a = a, b = b, mean = m, sd = de)
  
  # 2. Estimar con ambos métodos
  d_nrd0 <- density(X, bw = "nrd0")
  # Suprimimos los warnings de UCV para no saturar la consola en esta fase
  d_ucv <- suppressWarnings(density(X, bw = "ucv"))
  
  # 3. Interpolar las estimaciones al grid común
  y_nrd0 <- approx(d_nrd0$x, d_nrd0$y, xout = x_grid)$y
  y_ucv  <- approx(d_ucv$x, d_ucv$y, xout = x_grid)$y
  
  # Limpiar NAs por interpolación fuera del rango original de estimación
  y_nrd0[is.na(y_nrd0)] <- 0
  y_ucv[is.na(y_ucv)] <- 0
  
  # 4. GRÁFICO 1: Residuos (Sesgo Local)
  sesgo_nrd0 <- y_nrd0 - y_teo
  sesgo_ucv  <- y_ucv - y_teo
  y_lim_sesgo <- max(abs(c(sesgo_nrd0, sesgo_ucv)))
  
  plot(x_grid, sesgo_nrd0, type = "l", col = "blue", lwd = 2,
       ylim = c(-y_lim_sesgo, y_lim_sesgo),
       ylab = "Residuo (Est - Teo)", xlab = "x",
       main = paste("Sesgo Local (Análisis de Bordes):", titulo))
  lines(x_grid, sesgo_ucv, col = "darkorange", lwd = 2, lty = 2)
  abline(h = 0, col = "darkgray", lty = 1) # Línea de perfección
  abline(v = c(a, b), col = "red", lty = 3) # Límites de truncamiento
  
  legend("topright", legend = c("Sesgo nrd0", "Sesgo ucv", "Límites (a,b)"),
         col = c("blue", "darkorange", "red"), lty = c(1, 2, 3), lwd = 2, bty = "n")
  
  # 5. GRÁFICO 2: Error Cuadrático Local
  se_nrd0 <- sesgo_nrd0^2
  se_ucv  <- sesgo_ucv^2
  
  plot(x_grid, se_nrd0, type = "l", col = "blue", lwd = 2,
       ylab = "Error Cuadrático", xlab = "x",
       main = paste("Penalización de Error Cuadrático (SE):", titulo))
  lines(x_grid, se_ucv, col = "darkorange", lwd = 2, lty = 2)
  abline(v = c(a, b), col = "red", lty = 3)
  
  # Restaurar layout gráfico a la normalidad
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}


# ==============================================================================
# 5. EJECUCIÓN DE SIMULACIONES Y GRÁFICOS DE COMPARACIÓN
# ==============================================================================

# [Tus bloques de simulación originales se mantienen intactos. 
# A modo de ejemplo, te muestro cómo aplicar la nueva función al final del Bloque 3.]

# ------------------------------------------------------------------------------
# Bloque 3: Simulaciones con phi = 0.9 y normal truncada (0,1) en (-2,2)
# ------------------------------------------------------------------------------
X31 <- gnormtrunc_dep(nsim = 200, phi= 0.9, a = -2, b = 2, m = 0, de = 1)
par(mfrow=c(1,2))
graficar_densidades_trunc(X31, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X31")
graficar_densidades_trunc(X31, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X31")

X32 <- gnormtrunc_dep(nsim = 500, phi= 0.9, a = -2, b = 2, m = 0, de = 1)
par(mfrow=c(1,2))
graficar_densidades_trunc(X32, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X32")
graficar_densidades_trunc(X32, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X32")

X33 <- gnormtrunc_dep(nsim = 1000, phi= 0.9, a = -2, b = 2, m = 0, de = 1)
par(mfrow=c(1,2))
graficar_densidades_trunc(X33, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X33")
graficar_densidades_trunc(X33, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X33")

X34 <- gnormtrunc_dep(nsim = 2000, phi= 0.9, a = -2, b = 2, m = 0, de = 1)
par(mfrow=c(1,2))
graficar_densidades_trunc(X34, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X34")
graficar_densidades_trunc(X34, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X34")

# >>> Aplicación de la herramienta de diagnóstico visual a la muestra X34 <<<
diagnostico_error_local(X31, a = -2, b = 2, m = 0, de = 1, titulo = "X34 (n=2000, phi=0.9)")

# # Simulaciones con phi = 0 y normal truncada (0,1) en (-2,2)

X11 <- gnormtrunc_dep(nsim = 200, phi= 0, a = -2, b = 2, m = 0, de = 1)

# Estimaciones
dsilX11 <- density(X11, bw = "nrd0"); plot(dsilX11);
dcvX11 <- density(X11, bw = "ucv"); plot(dcvX11)

X12 <- gnormtrunc_dep(nsim = 500, phi= 0, a = -2, b = 2, m = 0, de = 1)

# Estimaciones
dsilX12 <- density(X12, bw = "nrd0"); plot(dsilX12)
dcvX12 <- density(X12, bw = "ucv"); plot(dcvX12)

X13 <- gnormtrunc_dep(nsim = 1000, phi= 0, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX13 <- density(X13, bw = "nrd0"); plot(dsilX13)
dcvX13 <- density(X13, bw = "ucv"); plot(dcvX13)

X14 <- gnormtrunc_dep(nsim = 2000, phi= 0, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX14 <- density(X14, bw = "nrd0"); plot(dsilX14)
dcvX14 <- density(X14, bw = "ucv"); plot(dcvX14)

# Simulaciones con phi = 0.5 y normal truncada (0,1) en (-2,2)
X21 <- gnormtrunc_dep(nsim = 200, phi= 0.5, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX21 <- density(X21, bw = "nrd0"); plot(dsilX21)
dcvX21 <- density(X21, bw = "ucv"); plot(dcvX21)

X22 <- gnormtrunc_dep(nsim = 500, phi= 0.5, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX22 <- density(X22, bw = "nrd0"); plot(dsilX22)
dcvX22 <- density(X22, bw = "ucv"); plot(dcvX22)

X23 <- gnormtrunc_dep(nsim = 1000, phi= 0.5, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX23 <- density(X23, bw = "nrd0"); plot(dsilX23)
dcvX23 <- density(X23, bw = "ucv"); plot(dcvX23)

X24 <- gnormtrunc_dep(nsim = 2000, phi= 0.5, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX24 <- density(X24, bw = "nrd0"); plot(dsilX24)
dcvX24 <- density(X24, bw = "ucv"); plot(dcvX24)

# Simulaciones con phi = 0.9 y normal truncada (0,1) en (-2,2)

X31 <- gnormtrunc_dep(nsim = 200, phi= 0.9, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX31 <- density(X31, bw = "nrd0"); plot(dsilX31)
dcvX31 <- density(X31, bw = "ucv"); plot(dcvX31)

X32 <- gnormtrunc_dep(nsim = 500, phi= 0.9, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX32 <- density(X32, bw = "nrd0"); plot(dsilX32)
dcvX32 <- density(X32, bw = "ucv"); plot(dcvX32)

X33 <- gnormtrunc_dep(nsim = 1000, phi= 0.9, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX33 <- density(X33, bw = "nrd0"); plot(dsilX33)
dcvX33 <- density(X33, bw = "ucv"); plot(dcvX33)

X34 <- gnormtrunc_dep(nsim = 2000, phi= 0.9, a = -2, b = 2, m = 0, de = 1)
# Estimaciones
dsilX34 <- density(X34, bw = "nrd0"); plot(dsilX34)
dcvX34 <- density(X34, bw = "ucv"); plot(dcvX34)
