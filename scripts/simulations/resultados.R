###############################################################################
# Script Modificado: Estimacion_Densidades_y_MonteCarlo.R
# Descripción: Integración de replicación individual, simulación Monte Carlo, 
# método GL matricial y muestreo eficiente de Normal Truncada (Transformada Inversa).
###############################################################################

# Carga de paquetes necesarios
if (!require(matrixStats)) install.packages("matrixStats")
if (!require(nor1mix)) install.packages("nor1mix")
if (!require(doParallel)) install.packages("doParallel")
if (!require(foreach)) install.packages("foreach")

library(matrixStats)
library(nor1mix)
library(doParallel)
library(foreach)
library(stats)

# =====================================================================
# 1. DEFINICIÓN DE DENSIDADES Y GENERADORES EFICIENTES
# =====================================================================

# Lista de las 15 densidades de Marron & Wand predefinidas
mw_list <- list(MW.nm1, MW.nm2, MW.nm3, MW.nm4, MW.nm5,
                MW.nm6, MW.nm7, MW.nm8, MW.nm9, MW.nm10,
                MW.nm11, MW.nm12, MW.nm13, MW.nm14, MW.nm15)

# Funciones de la Mezcla Bimodal
GMix <- function(x) { 0.5 * pnorm(x, mean = -2, sd = 1) + 0.5 * pnorm(x, mean = 2, sd = 1) }
grid.mix.x <- seq(-8, 8, length.out = 10000)
grid.mix.u <- GMix(grid.mix.x)
QMix <- function(u) { approx(x = grid.mix.u, y = grid.mix.x, xout = u, rule = 2)$y }

# Generación de variables uniformes bajo dependencia AR(1)
generar_uniformes_ar1 <- function(nsim, phi) {
  if (phi == 0) {
    Z <- rnorm(nsim)
  } else {
    Z <- arima.sim(list(ar = phi), n = nsim) 
  }
  sd_Z <- if (phi == 0) 1 else sqrt(1 / (1 - phi^2))
  return(pnorm(Z, mean = 0, sd = sd_Z))
}

# ---------------------------------------------------------------------
# NUEVO: Generador eficiente para Normal Truncada (Transformada Inversa)
# ---------------------------------------------------------------------
generar_truncnorm_dep <- function(n, phi, a, b, m = 0, de = 1) {
  U <- generar_uniformes_ar1(nsim = n, phi = phi)
  
  # Calculamos las probabilidades acumuladas en los extremos de truncamiento
  F_a <- pnorm(a, mean = m, sd = de)
  F_b <- pnorm(b, mean = m, sd = de)
  p <- F_b - F_a
  
  # Transformada Inversa adaptada al dominio truncado
  return(qnorm(p * U + F_a, mean = m, sd = de))
}

# Generador universal unificado
generar_muestra <- function(n, phi, tipo, mw_idx = 1, a = -2, b = 2) {
  U <- generar_uniformes_ar1(nsim = n, phi = phi)
  
  if (tipo == "truncnorm") {
    X <- generar_truncnorm_dep(n, phi, a, b, 0, 1)
    c_norm <- pnorm(b) - pnorm(a)
    gsup <- dnorm(0) / c_norm # Asumiendo m=0 y a < 0 < b para el supremo
    
    y_teo_func <- function(x) {
      y <- dnorm(x) / c_norm
      y[x < a | x > b] <- 0 # Densidad 0 fuera de los límites
      return(y)
    }
  } else if (tipo == "normal") {
    X <- qnorm(U)
    gsup <- 1 / sqrt(2 * pi)
    y_teo_func <- dnorm
  } else if (tipo == "lognormal") {
    X <- qlnorm(U, meanlog = 0, sdlog = 0.5)
    gsup <- exp(0.125) / (0.5 * sqrt(2 * pi))
    y_teo_func <- function(x) dlnorm(x, meanlog = 0, sdlog = 0.5)
  } else if (tipo == "mezcla") {
    X <- QMix(U)
    gsup <- 0.5 * dnorm(-2, -2, 1) + 0.5 * dnorm(-2, 2, 1)
    y_teo_func <- function(x) 0.5 * dnorm(x, -2, 1) + 0.5 * dnorm(x, 2, 1)
  } else if (tipo == "mw") {
    mw_obj <- mw_list[[mw_idx]]
    X <- qnorMix(U, mw_obj)
    opt <- optimize(function(x) dnorMix(x, mw_obj), interval = c(-5, 5), maximum = TRUE)
    gsup <- opt$objective
    y_teo_func <- function(x) dnorMix(x, mw_obj)
  } else {
    stop("Tipo de densidad no válido.")
  }
  return(list(X = X, gsup = gsup, y_teo_func = y_teo_func))
}


# =====================================================================
# 2. IMPLEMENTACIÓN MATRICIAL EFICIENTE DEL MÉTODO GL (Kernel Normal)
# =====================================================================

Vh_vectorizado <- function(H_vec, n, gsup, gamma) {
  delta.n <- sqrt(log(n))
  K2 <- sqrt(1 / (2 * sqrt(pi)))
  termino_ppal <- sqrt(2 * gamma * gsup) * K2 * (2) * (1 + delta.n)
  tasa <- (log(n)^(-1/2)) / sqrt(n * H_vec)
  return(termino_ppal * tasa)
}

GL_matrix_estimator <- function(X, gridcal, H, V_H) {
  n <- length(X)
  dist_base <- outer(gridcal, X, "-")
  
  Matriz_GH <- matrix(0, nrow = length(gridcal), ncol = length(H))
  for(j in seq_along(H)) {
    Matriz_GH[, j] <- rowSums(dnorm(dist_base / H[j])) / (n * H[j])
  }
  
  Criterio_A_Matriz <- matrix(0, nrow = length(gridcal), ncol = length(H))
  for(i in seq_along(H)) {
    h_actual <- H[i]
    GHH_temp <- matrix(0, nrow = length(gridcal), ncol = length(H))
    
    for(j in seq_along(H)) {
      hp_actual <- H[j]
      hs <- sqrt(h_actual^2 + hp_actual^2)
      GHH_temp[, j] <- rowSums(dnorm(dist_base / hs)) / (n * hs)
    }
    
    Diferencia <- abs(GHH_temp - Matriz_GH)
    Penalizado <- sweep(Diferencia, MARGIN = 2, STATS = V_H, FUN = "-")
    Penalizado[Penalizado < 0] <- 0 
    Criterio_A_Matriz[, i] <- rowMaxs(Penalizado)
  }
  
  Criterio_Final <- sweep(Criterio_A_Matriz, MARGIN = 2, STATS = V_H, FUN = "+")
  indices_minimos <- max.col(-Criterio_Final, ties.method = "first")
  
  return(Matriz_GH[cbind(1:length(gridcal), indices_minimos)])
}


# =====================================================================
# 3. FUNCIONES DE DIAGNÓSTICO Y REPLICACIÓN INDIVIDUAL
# =====================================================================

plot_composite <- function(gridcal, y_teo, dens_nrd0, dens_ucv, dens_sj, dens_gl, 
                           err_nrd0, err_ucv, err_sj, err_gl, titulo) {
  
  layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
  par(mar = c(4, 4, 3, 1) + 0.1)
  
  # --- Panel Superior: Densidades ---
  y_max <- max(c(y_teo, dens_nrd0, dens_ucv, dens_sj, dens_gl)) * 1.1
  plot(gridcal, y_teo, type = "l", col = "black", lwd = 3, lty = 2,
       ylim = c(0, y_max), xlab = "x", ylab = "Densidad", 
       main = paste("Comparación de Estimadores -", titulo))
  
  lines(gridcal, dens_nrd0, col = "#0072B2", lwd = 2) 
  lines(gridcal, dens_ucv, col = "#D55E00", lwd = 2)  
  lines(gridcal, dens_sj, col = "#009E73", lwd = 2)   
  lines(gridcal, dens_gl, col = "#CC79A7", lwd = 2)   
  
  legend("topright", legend = c("Teórica", "Silverman", "UCV", "S-J (Plug-in)", "GL"),
         col = c("black", "#0072B2", "#D55E00", "#009E73", "#CC79A7"),
         lty = c(2, 1, 1, 1, 1), lwd = c(3, 2, 2, 2, 2), bty = "n")
  
  # --- Panel Inferior: Error Cuadrático Local ---
  err_max <- max(c(err_nrd0, err_ucv, err_sj, err_gl))
  plot(gridcal, err_gl, type = "l", col = "#CC79A7", lwd = 2,
       ylim = c(0, err_max), xlab = "x", ylab = "Error Cuadrático",
       main = "Error Cuadrático Local")
  
  lines(gridcal, err_nrd0, col = "#0072B2", lwd = 2)
  lines(gridcal, err_ucv, col = "#D55E00", lwd = 2)
  lines(gridcal, err_sj, col = "#009E73", lwd = 2)
  
  layout(1) 
}

# ---------------------------------------------------------------------
# NUEVO: Función modular para una única replicación y gráfico automático
# ---------------------------------------------------------------------
analizar_y_graficar_replicacion <- function(n, phi, tipo_densidad, mw_idx = 1, 
                                            gamma_gl = 0.05, a = -2, b = 2, 
                                            generar_plot = TRUE) {
  
  # 1. Generación de muestra universal
  info_muestra <- generar_muestra(n = n, phi = phi, tipo = tipo_densidad, mw_idx = mw_idx, a = a, b = b)
  X <- info_muestra$X
  gsup <- info_muestra$gsup
  
  # 2. Configurar malla (Truncada utiliza malla fija [-4, 4] para consistencia general)
  lim_inf <- if(tipo_densidad == "truncnorm") -4 else min(X) - 1.5 * sd(X)
  lim_sup <- if(tipo_densidad == "truncnorm") 4 else max(X) + 1.5 * sd(X)
  n_grid <- 512
  
  # 3. Estimaciones Clásicas (Misma grilla forzada)
  d_nrd0 <- density(X, bw = "nrd0", from = lim_inf, to = lim_sup, n = n_grid)
  d_sj   <- density(X, bw = "SJ",   from = lim_inf, to = lim_sup, n = n_grid)
  d_ucv  <- tryCatch({
    suppressWarnings(density(X, bw = "ucv", from = lim_inf, to = lim_sup, n = n_grid))
  }, error = function(e) return(d_nrd0)) # Fallback a Silverman si falla la convergencia
  
  # 4. Estimación GL
  gridcal <- d_nrd0$x
  u <- seq(0, floor(log(n)) * (2/3), by = 0.1)
  H <- exp(-u)
  V_H <- Vh_vectorizado(H, n, gsup, gamma_gl)
  dens_gl <- GL_matrix_estimator(X, gridcal, H, V_H)
  
  # 5. Errores Cuadráticos e ISE
  y_teo <- info_muestra$y_teo_func(gridcal)
  delta_x <- gridcal[2] - gridcal[1]
  
  err_nrd0 <- (d_nrd0$y - y_teo)^2
  err_ucv  <- (d_ucv$y - y_teo)^2
  err_sj   <- (d_sj$y - y_teo)^2
  err_gl   <- (dens_gl - y_teo)^2
  
  ise_nrd0 <- sum(err_nrd0) * delta_x
  ise_ucv  <- sum(err_ucv) * delta_x
  ise_sj   <- sum(err_sj) * delta_x
  ise_gl   <- sum(err_gl) * delta_x
  
  # 6. Gráfico Compuesto
  if (generar_plot) {
    titulo <- paste(toupper(tipo_densidad), "(n =", n, ", phi =", phi, ")")
    plot_composite(gridcal, y_teo, d_nrd0$y, d_ucv$y, d_sj$y, dens_gl, 
                   err_nrd0, err_ucv, err_sj, err_gl, titulo)
  }
  
  # 7. Empaquetar y retornar resultados de forma transparente
  resultados <- list(
    Datos = X,
    Malla = gridcal,
    Densidad_Teorica = y_teo,
    Estimaciones = list(Silverman = d_nrd0$y, UCV = d_ucv$y, SJ = d_sj$y, GL = dens_gl),
    ISE = c(Silverman = ise_nrd0, UCV = ise_ucv, SJ = ise_sj, GL = ise_gl)
  )
  
  return(invisible(resultados))
}


# =====================================================================
# 4. ESTUDIO DE SIMULACIÓN MONTE CARLO (B = 1000)
# =====================================================================

ejecutar_montecarlo <- function(B = 1000, n = 200, phi = 0, tipo_densidad = "truncnorm", mw_idx = 1, gamma_gl = 0.05) {
  
  # Prevenir sobrecalentamiento dejando núcleos libres durante el proceso intensivo
  num_cores <- max(1, parallel::detectCores() - 2)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  cat(sprintf("Iniciando Monte Carlo (%d iteraciones). Densidad: %s. Núcleos: %d\n", B, tipo_densidad, num_cores))
  
  # Integración de la función modular dentro de `foreach`
  errores_mc <- foreach(m = 1:B, .combine = 'rbind', 
                        .packages = c("stats", "matrixStats", "nor1mix"),
                        .export = c("analizar_y_graficar_replicacion", "generar_muestra", "generar_truncnorm_dep",
                                    "generar_uniformes_ar1", "GL_matrix_estimator", "Vh_vectorizado", 
                                    "GMix", "QMix", "grid.mix.x", "grid.mix.u", "mw_list")) %dopar% {
                                      
                                      # Se ejecuta la función en modo silente (sin gráfico) para recolectar el ISE de cada método
                                      res_iter <- analizar_y_graficar_replicacion(n = n, phi = phi, tipo_densidad = tipo_densidad, 
                                                                                  mw_idx = mw_idx, gamma_gl = gamma_gl, 
                                                                                  generar_plot = FALSE)
                                      return(res_iter$ISE)
                                    }
  
  stopCluster(cl)
  
  mise_res <- colMeans(errores_mc)
  return(list(MISE = mise_res, Errores_Crudos = errores_mc))
}

# =====================================================================
# EJEMPLOS DE USO 
# =====================================================================

# 1. Ejecutar una única replicación y graficar para Densidad Normal Truncada
replicacion_unica <- analizar_y_graficar_replicacion(n = 60, phi = 0, tipo_densidad = "mezcla")
print(replicacion_unica$ISE)

# 2. Ejecutar Estudio Monte Carlo Completo usando las funciones modulares
# res_mc_trunc <- ejecutar_montecarlo(B = 1000, n = 200, phi = 0.5, tipo_densidad = "truncnorm")
# print(res_mc_trunc$MISE)