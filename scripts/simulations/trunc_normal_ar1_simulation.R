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

# Esteblciendo semilla inicial
set.seed(1234)

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
# 3. FUNCIÓN PARA ESTIMAR LA DENSIDAD
# ==============================================================================

est_dens_func <- function(X, bw_metodo, lim_inf, lim_sup) {
  
  # 1. Calcular estimación de densidad kernel
  # Nota: Es buena práctica fijar 'n' (ej. n = 512) para asegurar el mismo tamaño siempre
  est_dens <- density(X, bw = bw_metodo, from = lim_inf, to = lim_sup, n = 512)
  
  # 2. Extraer el valor numérico de h (ancho de banda) redondeado a 4 decimales
  h_estimado <- round(est_dens$bw, 4)
  
  # 3. Extraer la malla con la que se construyó la estimación
  malla <- est_dens$x
  
  # 4. Extraer las estimaciones de la densidad de la malla
  densidades <- est_dens$y
  
  # 5. Retornar las variables en una lista con nombres descriptivos
  return(list(
    h_estimado = h_estimado,
    malla      = malla,
    densidades = densidades
  ))
}

# ==============================================================================
# 3. FUNCIÓN ASISTENTE PARA SOBREPONER DENSIDADES (CORREGIDA Y AUTOMATIZADA)
# ==============================================================================
graficar_densidades_trunc <- function(densidad_est, bw_metodo, a, b, m, de, titulo = "") {
  
  # 1. Extraer los límites dinámicamente desde la malla ya calculada
  limite_inferior <- min(densidad_est$malla)
  limite_superior <- max(densidad_est$malla)
  
  # 2. Graficar la estimación del Kernel usando los vectores x e y
  plot(x = densidad_est$malla, 
       y = densidad_est$densidades, 
       type = "l", 
       main = paste("Densidad", titulo, "[", bw_metodo, "]"), 
       sub = paste("Ancho de banda estimado (h) =", round(densidad_est$h_estimado, 4)),
       col = "blue", lwd = 2, 
       xlim = c(limite_inferior, limite_superior), # Límites inyectados automáticamente
       ylim = c(0, max(densidad_est$densidades) * 1.2),
       xlab = "Valores de X", ylab = "Densidad")
  
  # 3. Sobreponer la teórica usando el paquete truncnorm
  curve(dtruncnorm(x, a = a, b = b, mean = m, sd = de), 
        col = "darkred", lwd = 3, lty = 3, add = TRUE)
  
  # 4. Leyenda explicativa enriquecida con el valor de h
  legend("topright", 
         legend = c(paste("Kernel (h =", round(densidad_est$h_estimado, 4), ")"), "Teórica (truncnorm)"), 
         col = c("blue", "darkred"), lty = c(1, 2), lwd = 2, bty = "n")
}



# ==============================================================================
# 4. DIAGNÓSTICO DE ERROR LOCAL (SESGO Y ERROR CUADRÁTICO) (ACTUALIZADO)
# ==============================================================================
diagnostico_error_local <- function(X, a, b, m, de, lim_inf, lim_sup, titulo = "") {
  # Configurar layout para 2 gráficos apilados verticalmente
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  
  # 1. Estimar con ambos métodos forzando la grilla exacta (sin interpolación)
  # Usamos from, to y n=512 para "congelar" la malla
  d_nrd0 <- density(X, bw = "nrd0", from = lim_inf, to = lim_sup, n = 512)
  
  # Suprimimos los warnings de UCV para no saturar la consola en esta fase
  d_ucv  <- suppressWarnings(density(X, bw = "ucv", from = lim_inf, to = lim_sup, n = 512))
  
  # 2. Extraer la malla exacta generada por R y evaluar la densidad teórica ahí
  x_grid <- d_nrd0$x
  y_teo <- dtruncnorm(x_grid, a = a, b = b, mean = m, sd = de)
  
  # 3. GRÁFICO 1: Residuos (Sesgo Local)
  # Calculamos el error directamente sin usar approx()
  sesgo_nrd0 <- d_nrd0$y - y_teo
  sesgo_ucv  <- d_ucv$y - y_teo
  
  # Establecer límites simétricos para el eje Y del sesgo
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
  
  # 4. GRÁFICO 2: Error Cuadrático Local
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

# ------------------------------------------------------------------------------
# Bloque 1: Simulaciones con phi = 0 y normal truncada (0,1) en (-2,2)
# ------------------------------------------------------------------------------

# --- X11 (n = 200) ---
X11 <- gnormtrunc_dep(nsim = 200, phi = 0, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X11 <- est_dens_func(X11, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X11  <- suppressWarnings(est_dens_func(X11, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X11, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X11")
graficar_densidades_trunc(dens_ucv_X11, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X11")
diagnostico_error_local(X11, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X11 (n=200, phi=0)")

# --- X12 (n = 500) ---
X12 <- gnormtrunc_dep(nsim = 500, phi = 0, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X12 <- est_dens_func(X12, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X12  <- suppressWarnings(est_dens_func(X12, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X12, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X12")
graficar_densidades_trunc(dens_ucv_X12, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X12")
diagnostico_error_local(X12, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X12 (n=500, phi=0)")

# --- X13 (n = 1000) ---
X13 <- gnormtrunc_dep(nsim = 1000, phi = 0, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X13 <- est_dens_func(X13, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X13  <- suppressWarnings(est_dens_func(X13, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X13, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X13")
graficar_densidades_trunc(dens_ucv_X13, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X13")
diagnostico_error_local(X13, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X13 (n=1000, phi=0)")

# --- X14 (n = 2000) ---
X14 <- gnormtrunc_dep(nsim = 2000, phi = 0, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X14 <- est_dens_func(X14, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X14  <- suppressWarnings(est_dens_func(X14, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X14, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X14")
graficar_densidades_trunc(dens_ucv_X14, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X14")
diagnostico_error_local(X14, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X14 (n=2000, phi=0)")


# ------------------------------------------------------------------------------
# Bloque 2: Simulaciones con phi = 0.5 y normal truncada (0,1) en (-2,2)
# ------------------------------------------------------------------------------

# --- X21 (n = 200) ---
X21 <- gnormtrunc_dep(nsim = 200, phi = 0.5, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X21 <- est_dens_func(X21, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X21  <- suppressWarnings(est_dens_func(X21, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X21, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X21")
graficar_densidades_trunc(dens_ucv_X21, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X21")
diagnostico_error_local(X21, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X21 (n=200, phi=0.5)")

# --- X22 (n = 500) ---
X22 <- gnormtrunc_dep(nsim = 500, phi = 0.5, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X22 <- est_dens_func(X22, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X22  <- suppressWarnings(est_dens_func(X22, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X22, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X22")
graficar_densidades_trunc(dens_ucv_X22, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X22")
diagnostico_error_local(X22, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X22 (n=500, phi=0.5)")

# --- X23 (n = 1000) ---
X23 <- gnormtrunc_dep(nsim = 1000, phi = 0.5, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X23 <- est_dens_func(X23, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X23  <- suppressWarnings(est_dens_func(X23, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X23, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X23")
graficar_densidades_trunc(dens_ucv_X23, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X23")
diagnostico_error_local(X23, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X23 (n=1000, phi=0.5)")

# --- X24 (n = 2000) ---
X24 <- gnormtrunc_dep(nsim = 2000, phi = 0.5, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X24 <- est_dens_func(X24, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X24  <- suppressWarnings(est_dens_func(X24, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X24, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X24")
graficar_densidades_trunc(dens_ucv_X24, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X24")
diagnostico_error_local(X24, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X24 (n=2000, phi=0.5)")


# ------------------------------------------------------------------------------
# Bloque 3: Simulaciones con phi = 0.9 y normal truncada (0,1) en (-2,2)
# ------------------------------------------------------------------------------

# --- X31 (n = 200) ---
X31 <- gnormtrunc_dep(nsim = 200, phi = 0.9, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X31 <- est_dens_func(X31, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X31  <- suppressWarnings(est_dens_func(X31, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X31, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X31")
graficar_densidades_trunc(dens_ucv_X31, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X31")
diagnostico_error_local(X31, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X31 (n=200, phi=0.9)")

# --- X32 (n = 500) ---
X32 <- gnormtrunc_dep(nsim = 500, phi = 0.9, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X32 <- est_dens_func(X32, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X32  <- suppressWarnings(est_dens_func(X32, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X32, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X32")
graficar_densidades_trunc(dens_ucv_X32, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X32")
diagnostico_error_local(X32, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X32 (n=500, phi=0.9)")

# --- X33 (n = 1000) ---
X33 <- gnormtrunc_dep(nsim = 1000, phi = 0.9, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X33 <- est_dens_func(X33, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X33  <- suppressWarnings(est_dens_func(X33, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X33, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X33")
graficar_densidades_trunc(dens_ucv_X33, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X33")
diagnostico_error_local(X33, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X33 (n=1000, phi=0.9)")

# --- X34 (n = 2000) ---
X34 <- gnormtrunc_dep(nsim = 2000, phi = 0.9, a = -2, b = 2, m = 0, de = 1)
dens_nrd0_X34 <- est_dens_func(X34, "nrd0", lim_inf = -3, lim_sup = 3)
dens_ucv_X34  <- suppressWarnings(est_dens_func(X34, "ucv", lim_inf = -3, lim_sup = 3))

par(mfrow=c(1,2))
graficar_densidades_trunc(dens_nrd0_X34, "nrd0", a = -2, b = 2, m = 0, de = 1, titulo = "X34")
graficar_densidades_trunc(dens_ucv_X34, "ucv",  a = -2, b = 2, m = 0, de = 1, titulo = "X34")
diagnostico_error_local(X34, a = -2, b = 2, m = 0, de = 1, lim_inf = -3, lim_sup = 3, titulo = "X34 (n=2000, phi=0.9)")