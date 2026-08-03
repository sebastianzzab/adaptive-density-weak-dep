# ==============================================================================
# FUNCIÓN DE SIMULACIÓN (Normal Truncada Dependiente AR(1))
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
  
  return(X)
}

# =====================================================================
# Función Generadora mediante Transformada Inversa Numérica
# =====================================================================
# Dado que la distribución empírica en nor1mix no está truncada,
# ampliamos el intervalo de búsqueda a valores seguros (ej. -15 a 15).
generar_dependiente_mw <- function(nsim, phi, mw_obj, limite = c(-15,15)) {
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
  
  # Usamos mclapply (o sapply) para optimizar el cálculo de raíces en vectores largos.
  # El método uniroot invierte la función pnorMix (CDF de la mezcla).
  Y_list <- lapply(U, function(u) {
    res <- uniroot(function(x) pnorMix(x, mw_obj) - u, interval = limite)
    return(res$root)
  })
  
  return(unlist(Y_list))
}

# ==============================================================================
# FUNCIÓN ASISTENTE PARA SOBREPONER DENSIDADES
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
# FUNCIÓN ASISTENTE PARA SOBREPONER DENSIDADES DE MARRON & WAND
# ==============================================================================
graficar_densidades_mw <- function(densidad_est, densidad_objetivo, bw_metodo, titulo = "") {
  
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
  curve(dnorMix(x, densidad_objetivo),
        col = "darkred", lwd = 3, lty = 3, add = TRUE)
  # curve(dtruncnorm(x, a = a, b = b, mean = m, sd = de), 
  #       col = "darkred", lwd = 3, lty = 3, add = TRUE)
  
  # 4. Leyenda explicativa enriquecida con el valor de h
  legend("topright", 
         legend = c(paste("Kernel (h =", round(densidad_est$h_estimado, 4), ")"), "Teórica (nor1mix)"), 
         col = c("blue", "darkred"), lty = c(1, 2), lwd = 2, bty = "n")
}

# ==============================================================================
# FUNCIÓN DE DIAGNÓSTICO DE ERROR LOCAL (ERROR CUADRÁTICO)
# ==============================================================================
diagnostico_error_local <- function(X, X_nrd0, X_ucv, a, b, m, de, lim_inf, lim_sup, titulo = "") {
  # Configurar layout para 2 gráficos apilados verticalmente
  # par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  
  # # 1. Estimar con ambos métodos forzando la grilla exacta (sin interpolación)
  # # Usamos from, to y n=512 para "congelar" la malla
  # d_nrd0 <- density(X, bw = "nrd0", from = lim_inf, to = lim_sup, n = 512)
  # 
  # # Suprimimos los warnings de UCV para no saturar la consola en esta fase
  # d_ucv  <- suppressWarnings(density(X, bw = "ucv", from = lim_inf, to = lim_sup, n = 512))
  
  # 2. Extraer la malla exacta generada por R y evaluar la densidad teórica ahí
  x_grid <- X_nrd0$malla
  y_teo <- dtruncnorm(x_grid, a = a, b = b, mean = m, sd = de)
  
  # 3. GRÁFICO 1: Residuos (Sesgo Local)
  # Calculamos el error directamente sin usar approx()
  sesgo_nrd0 <- X_nrd0$densidades - y_teo
  sesgo_ucv <- X_ucv$densidades - y_teo
  
  # sesgo_nrd0 <- d_nrd0$y - y_teo
  # sesgo_ucv  <- d_ucv$y - y_teo
  
  # Establecer límites simétricos para el eje Y del sesgo
  y_lim_sesgo <- max(abs(c(sesgo_nrd0, sesgo_ucv)))
  
  # plot(x_grid, sesgo_nrd0, type = "l", col = "blue", lwd = 2,
  #      ylim = c(-y_lim_sesgo, y_lim_sesgo),
  #      ylab = "Residuo (Est - Teo)", xlab = "x",
  #      main = paste("Sesgo Local (Análisis de Bordes):", titulo))
  # lines(x_grid, sesgo_ucv, col = "darkorange", lwd = 2, lty = 2)
  # abline(h = 0, col = "darkgray", lty = 1) # Línea de perfección
  # abline(v = c(a, b), col = "red", lty = 3) # Límites de truncamiento
  # 
  # legend("topright", legend = c("Sesgo nrd0", "Sesgo ucv", "Límites (a,b)"),
  #        col = c("blue", "darkorange", "red"), lty = c(1, 2, 3), lwd = 2, bty = "n")
  
  # 4. GRÁFICO 2: Error Cuadrático Local
  se_nrd0 <- sesgo_nrd0^2
  se_ucv  <- sesgo_ucv^2
  
  plot(x_grid, se_nrd0, type = "l", col = "blue", lwd = 2,
       ylab = "Error Cuadrático", xlab = "x",
       main = paste("Penalización de Error Cuadrático (SE):", titulo))
  lines(x_grid, se_ucv, col = "darkorange", lwd = 2, lty = 2)
  abline(v = c(a, b), col = "red", lty = 3)
  legend("topright", legend = c("nrd0", "ucv", "Límites (a,b)"),
         col = c("blue", "darkorange", "red"), lty = c(1, 2, 3), lwd = 2, bty = "n")
  
  # Restaurar layout gráfico a la normalidad
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

# ==============================================================================
# FUNCIÓN DE DIAGNÓSTICO DE ERROR LOCAL (ERROR CUADRÁTICO)
# ==============================================================================
diagnostico_error_localmw <- function(X, X_nrd0, X_ucv, densidad_objetivo,lim_inf, lim_sup, titulo = "") {
  # Configurar layout para 2 gráficos apilados verticalmente
  # par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  
  # # 1. Estimar con ambos métodos forzando la grilla exacta (sin interpolación)
  # # Usamos from, to y n=512 para "congelar" la malla
  # d_nrd0 <- density(X, bw = "nrd0", from = lim_inf, to = lim_sup, n = 512)
  # 
  # # Suprimimos los warnings de UCV para no saturar la consola en esta fase
  # d_ucv  <- suppressWarnings(density(X, bw = "ucv", from = lim_inf, to = lim_sup, n = 512))
  
  # 2. Extraer la malla exacta generada por R y evaluar la densidad teórica ahí
  x_grid <- X_nrd0$malla
  y_teo <- dnorMix(x_grid, densidad_objetivo)
  # y_teo <- dtruncnorm(x_grid, a = a, b = b, mean = m, sd = de)
  
  # 3. GRÁFICO 1: Residuos (Sesgo Local)
  # Calculamos el error directamente sin usar approx()
  sesgo_nrd0 <- X_nrd0$densidades - y_teo
  sesgo_ucv <- X_ucv$densidades - y_teo
  
  # sesgo_nrd0 <- d_nrd0$y - y_teo
  # sesgo_ucv  <- d_ucv$y - y_teo
  
  # Establecer límites simétricos para el eje Y del sesgo
  y_lim_sesgo <- max(abs(c(sesgo_nrd0, sesgo_ucv)))
  
  # plot(x_grid, sesgo_nrd0, type = "l", col = "blue", lwd = 2,
  #      ylim = c(-y_lim_sesgo, y_lim_sesgo),
  #      ylab = "Residuo (Est - Teo)", xlab = "x",
  #      main = paste("Sesgo Local (Análisis de Bordes):", titulo))
  # lines(x_grid, sesgo_ucv, col = "darkorange", lwd = 2, lty = 2)
  # abline(h = 0, col = "darkgray", lty = 1) # Línea de perfección
  # abline(v = c(a, b), col = "red", lty = 3) # Límites de truncamiento
  # 
  # legend("topright", legend = c("Sesgo nrd0", "Sesgo ucv", "Límites (a,b)"),
  #        col = c("blue", "darkorange", "red"), lty = c(1, 2, 3), lwd = 2, bty = "n")
  
  # 4. GRÁFICO 2: Error Cuadrático Local
  se_nrd0 <- sesgo_nrd0^2
  se_ucv  <- sesgo_ucv^2
  
  plot(x_grid, se_nrd0, type = "l", col = "blue", lwd = 2,
       ylab = "Error Cuadrático", xlab = "x",
       main = paste("Penalización de Error Cuadrático (SE):", titulo))
  lines(x_grid, se_ucv, col = "darkorange", lwd = 2, lty = 2)
  abline(v = c(lim_inf, lim_sup), col = "red", lty = 3)
  legend("topright", legend = c("nrd0", "ucv", "Límites (a,b)"),
         col = c("blue", "darkorange", "red"), lty = c(1, 2, 3), lwd = 2, bty = "n")
  
  # Restaurar layout gráfico a la normalidad
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}