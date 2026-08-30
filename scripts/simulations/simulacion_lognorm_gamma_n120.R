###############################################################################
# Script: Estudio de simulación: Caso Lognormal y Gamma con n=120
# Realizado por: Sebastian Zabala
# Fecha: Septiembre 2026
###############################################################################

setwd("~/Desktop/tesis_sebastian/adaptive-density-weak-dep/scripts/calibration_study/PRP_Lognormal")

# Instalación de paquetes necesarios
if (!require(matrixStats)) install.packages("matrixStats")
if (!require(nor1mix)) install.packages("nor1mix")
if (!require(doParallel)) install.packages("doParallel")
if (!require(foreach)) install.packages("foreach")
if (!require(doRNG)) install.packages("doRNG")
if (!require(xtable)) install.packages("xtable")

# Carga de paquetes necesarios
library(xtable)
library(doRNG)
library(matrixStats)
library(nor1mix)
library(doParallel)
library(foreach)
library(stats)
library(tidyr)

# =====================================================================
# DEFINICIÓN DE DENSIDADES Y GENERADORES EFICIENTES
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
generar_uniformes_ar1 <- function(nsim, phi, me = 0, de = 1) {
  if (phi == 0) {
    Z <- rnorm(nsim)
  } else {
    Z <- arima.sim(list(ar = phi), n = nsim) 
  }
  sd_Z <- if (phi == 0) 1 else sqrt(1 / (1 - phi^2))
  return(pnorm(Z, mean = 0, sd = sd_Z))
}

# Generador para Normal Truncada mediante Transformada Inversa
generar_truncnorm_dep <- function(n, phi, a, b, m = 0, de = 1) {
  U <- generar_uniformes_ar1(nsim = n, phi = phi)
  
  # Calculamos las probabilidades acumuladas en los extremos de truncamiento
  F_a <- pnorm(a, mean = m, sd = de)
  F_b <- pnorm(b, mean = m, sd = de)
  p <- F_b - F_a
  
  # Transformada Inversa adaptada al dominio truncado
  return(qnorm(p * U + F_a, mean = m, sd = de))
}

# Generador de muestras
generar_muestra <- function(n, phi, tipo, mw_idx = 1, a = -2, b = 2, me = 0, de = 1) {
  U <- generar_uniformes_ar1(nsim = n, phi = phi, me, de)
  
  if (tipo == "truncnorm") {
    X <- generar_truncnorm_dep(n, phi, a, b, 0, 1)
    c_norm <- pnorm(b) - pnorm(a)
    gsup <- dnorm(0) / c_norm 
    
    y_teo_func <- function(x) {
      y <- dnorm(x) / c_norm
      y[x < a | x > b] <- 0 
      return(y)
    }
  } else if (tipo == "normal-est") {
    X <- qnorm(U)
    gsup <- 1 / sqrt(2 * pi)
    y_teo_func <- dnorm
    lim_inf <- -8
    lim_sup <- 8
  } else if (tipo == "normal") {
    X <- qnorm(U, mean = me, sd = de)                      
    gsup <- 1 / (sqrt(2 * pi) * de)
    y_teo_func <- function(x) dnorm(x, mean = me, sd = de)
    lim_inf <- -8
    lim_sup <- 8
  } else if (tipo == "lognormal") {
    X <- qlnorm(U, meanlog = 0, sdlog = 0.5)
    gsup <- exp(0.125) / (0.5 * sqrt(2 * pi))
    y_teo_func <- function(x) dlnorm(x, meanlog = 0, sdlog = 0.5)
    lim_inf <- 0
    lim_sup <- qlnorm(0.999999, 0, 0.5)
  } else if (tipo == "mezcla") {
    X <- QMix(U)
    gsup <- 0.5 * dnorm(-2, -2, 1) + 0.5 * dnorm(-2, 2, 1)
    y_teo_func <- function(x) 0.5 * dnorm(x, -2, 1) + 0.5 * dnorm(x, 2, 1)
    lim_inf <- -8
    lim_sup <- 8
  } else if (tipo == "mw") {
    mw_obj <- mw_list[[mw_idx]]
    X <- qnorMix(U, mw_obj)
    opt <- optimize(function(x) dnorMix(x, mw_obj), interval = c(-5, 5), maximum = TRUE)
    gsup <- opt$objective
    y_teo_func <- function(x) dnorMix(x, mw_obj)
    lim_inf <- -8
    lim_sup <- 8
  } else if (tipo == "logistica") {
    s <- sqrt(3) / pi
    # Cambiamos location a 3
    X <- qlogis(U, location = 3, scale = s)
    
    # El valor máximo (gsup) ocurre en la media
    gsup <- dlogis(3, location = 3, scale = s) 
    y_teo_func <- function(x) dlogis(x, location = 3, scale = s)
    lim_inf <- -8
    lim_sup <- 8
  } else if (tipo == "laplace") {
    b_param <- 1 / sqrt(2)
    U_centrada <- U - 0.5 
    
    # Sumamos 3 para desplazar la media
    X <- 3 - b_param * sign(U_centrada) * log(1 - 2 * abs(U_centrada))
    
    # El valor máximo se mantiene igual, pero actualizamos la función teórica
    gsup <- 1 / (2 * b_param)
    y_teo_func <- function(x) (1 / (2 * b_param)) * exp(-abs(x - 3) / b_param)
    lim_inf <- -8
    lim_sup <- 8
  } else if (tipo == "gamma") {
    # 1. Emparejamos momentos con la Lognormal(0, 0.5)
    media_ln <- exp(0.125)
    var_ln <- (exp(0.25) - 1) * exp(0.25)
    
    # 2. Despejamos los parámetros shape y rate de la Gamma
    shape_g <- (media_ln^2) / var_ln
    rate_g <- media_ln / var_ln
    
    X <- qgamma(U, shape = shape_g, rate = rate_g)
    
    # El supremo de la densidad Gamma ocurre en su moda: (shape - 1)/rate
    moda_g <- (shape_g - 1) / rate_g
    gsup <- dgamma(moda_g, shape = shape_g, rate = rate_g)
    y_teo_func <- function(x) dgamma(x, shape = shape_g, rate = rate_g)
    lim_inf <- -8
    lim_sup <- 8
  } else if (tipo == "weibull") {
    # 1. Fijamos una forma asimétrica positiva fuerte (shape = 1.5)
    shape_w <- 1.5
    media_ln <- exp(0.125)
    
    # 2. Ajustamos la escala para que coincida exactamente con la media de la Lognormal
    scale_w <- media_ln / gamma(1 + 1/shape_w)
    
    X <- qweibull(U, shape = shape_w, scale = scale_w)
    
    # El supremo de la Weibull ocurre en su moda
    moda_w <- scale_w * ((shape_w - 1)/shape_w)^(1/shape_w)
    gsup <- dweibull(moda_w, shape = shape_w, scale = scale_w)
    y_teo_func <- function(x) dweibull(x, shape = shape_w, scale = scale_w)
    lim_inf <- -8
    lim_sup <- 8
  }else {
    stop("Tipo de densidad no válido.")
  }
  return(list(X = X, gsup = gsup, y_teo_func = y_teo_func, lim_inf = lim_inf, lim_sup = lim_sup))
}

# =====================================================================
# OBTENCIÓN DE PUNTOS CRÍTICOS TEÓRICOS
# =====================================================================

obtener_puntos_evaluacion <- function(tipo, mw_idx = 1, me = 0, de = 1) {
  nombres <- c("Moda", "P05", "Q1", "Mediana", "Q3", "P95")
  pts <- numeric(6)
  
  if (tipo == "normal") {
    pts[1] <- me
    pts[2:6] <- qnorm(c(0.05, 0.25, 0.50, 0.75, 0.95), mean = me, sd = de)
    
  } else if (tipo == "lognormal") {
    pts[1] <- exp(0 - 0.5^2) # Moda de LN(0, 0.5)
    pts[2:6] <- qlnorm(c(0.05, 0.25, 0.50, 0.75, 0.95), meanlog = 0, sdlog = 0.5)
    
  } else if (tipo == "logistica") {
    s <- sqrt(3) / pi
    pts[1] <- 3
    pts[2:6] <- qlogis(c(0.05, 0.25, 0.50, 0.75, 0.95), location = 3, scale = s)
    
  } else if (tipo == "laplace") {
    b_param <- 1 / sqrt(2)
    pts[1] <- 3
    # Función cuantil inversa para Laplace
    q_laplace <- function(p, m, b) {
      ifelse(p < 0.5, m + b * log(2 * p), m - b * log(2 * (1 - p)))
    }
    pts[2:6] <- q_laplace(c(0.05, 0.25, 0.50, 0.75, 0.95), m = 3, b = b_param)
    
  } else if (tipo == "gamma") {
    media_ln <- exp(0.125); var_ln <- (exp(0.25) - 1) * exp(0.25)
    shape_g <- (media_ln^2) / var_ln; rate_g <- media_ln / var_ln
    pts[1] <- (shape_g - 1) / rate_g
    pts[2:6] <- qgamma(c(0.05, 0.25, 0.50, 0.75, 0.95), shape = shape_g, rate = rate_g)
    
  } else if (tipo == "weibull") {
    shape_w <- 1.5; media_ln <- exp(0.125)
    scale_w <- media_ln / gamma(1 + 1/shape_w)
    pts[1] <- scale_w * ((shape_w - 1)/shape_w)^(1/shape_w)
    pts[2:6] <- qweibull(c(0.05, 0.25, 0.50, 0.75, 0.95), shape = shape_w, scale = scale_w)
    
  } else if (tipo == "mezcla") {
    pts[1] <- 2 # Moda del componente derecho
    pts[2:6] <- QMix(c(0.05, 0.25, 0.50, 0.75, 0.95))
    
  } else if (tipo == "mw") {
    # Cálculo de puntos críticos para familia Marron-Wand (Garra)
    mw_obj <- mw_list[[mw_idx]]
    opt <- optimize(function(x) dnorMix(x, mw_obj), interval = c(-5, 5), maximum = TRUE)
    pts[1] <- opt$maximum # La moda (Pico más alto)
    pts[2:6] <- qnorMix(c(0.05, 0.25, 0.50, 0.75, 0.95), mw_obj)
  }
  
  names(pts) <- nombres
  return(pts)
}


# =====================================================================
# IMPLEMENTACIÓN MATRICIAL EFICIENTE DEL MÉTODO GL (Kernel Normal)
# =====================================================================

Vh_vectorizado <- function(H_vec, n, gsup, gamma, phi) {
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
  
  densidad_estimada <- Matriz_GH[cbind(1:length(gridcal), indices_minimos)]
  h_seleccionadas <- H[indices_minimos] # La ventana usada en cada punto
  
  return(list(densidad = densidad_estimada, h_local = h_seleccionadas))
}


# =====================================================================
# FUNCIONES DE DIAGNÓSTICO Y REPLICACIÓN INDIVIDUAL
# =====================================================================

plot_composite <- function(X, gridcal, y_teo, dens_nrd0, h_nrd0, dens_ucv, h_ucv, dens_sj, h_sj, dens_gl, h_gl, 
                           err_nrd0, err_ucv, err_sj, err_gl, titulo) {
  
  layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
  par(mar = c(4, 4, 3, 1) + 0.1)
  
  # --- Panel Superior: Densidades con Histograma ---
  y_max <- max(c(y_teo, dens_nrd0, dens_ucv, dens_sj, dens_gl)) * 1.1
  
  # 1. Dibujamos el histograma de fondo primero (prob = TRUE es crucial para que coincida con las densidades)
  hist(X, prob = TRUE, breaks = 30, col = "gray90", border = "white",
       xlim = range(gridcal), ylim = c(0, y_max), 
       main = "", xlab = "x", ylab = "Densidad")
  
  # 2. Agregamos la densidad teórica (ahora con lines() en vez de plot())
  lines(gridcal, y_teo, col = "black", lwd = 3, lty = 2)
  
  # 3. Agregamos el resto de las estimaciones
  lines(gridcal, dens_nrd0, col = "#0072B2", lwd = 2) 
  lines(gridcal, dens_ucv, col = "#D55E00", lwd = 2)  
  lines(gridcal, dens_sj, col = "#009E73", lwd = 2)   
  lines(gridcal, dens_gl, col = "blue", lwd = 2)   
  
  legend("topright", 
         legend = as.expression(c(
           bquote("Teórica"),
           bquote("ROT (h = " * .(round(h_nrd0, 2)) * ")"),
           bquote("UCV (h = " * .(round(h_ucv, 2)) * ")"),
           bquote("SJ (h = " * .(round(h_sj, 2)) * ")"),
           bquote("GL (" * hat(h) * " = " * .(round(h_gl, 4)) * ")")
         )),
         col = c("black", "#0072B2", "#D55E00", "#009E73", "blue"),
         lty = c(2, 1, 1, 1, 1), 
         lwd = c(3, 2, 2, 2, 2), 
         bty = "n")
  
  # --- Panel Inferior: Error Cuadrático Local ---
  err_max <- max(c(err_nrd0, err_ucv, err_sj, err_gl))
  plot(gridcal, err_gl, type = "l", col = "blue", lwd = 2,
       ylim = c(0, err_max), xlab = "x", ylab = "Error Cuadrático Local")
  
  lines(gridcal, err_nrd0, col = "#0072B2", lwd = 2)
  lines(gridcal, err_ucv, col = "#D55E00", lwd = 2)
  lines(gridcal, err_sj, col = "#009E73", lwd = 2)
  
  # === Cálculo del ISE para la leyenda ===
  delta_x <- gridcal[2] - gridcal[1]
  ise_nrd0 <- sum(err_nrd0) * delta_x
  ise_ucv  <- sum(err_ucv) * delta_x
  ise_sj   <- sum(err_sj) * delta_x
  ise_gl   <- sum(err_gl) * delta_x
  
  legend("topright", 
         legend = c( 
           paste0("ROT (ISE = ", round(ise_nrd0, 4), ")"), 
           paste0("UCV (ISE = ", round(ise_ucv, 4), ")"), 
           paste0("SJ (ISE = ", round(ise_sj, 4), ")"), 
           paste0("GL (ISE = ", round(ise_gl, 4), ")")
         ),
         col = c( "#0072B2", "#D55E00", "#009E73", "blue"),
         lty = c(1, 1, 1, 1), 
         lwd = c(2, 2, 2, 2), 
         bty = "n")
  
  layout(1) 
}

# ---------------------------------------------------------------------
# Función modular para una única replicación y gráfico automático
# ---------------------------------------------------------------------
evaluar_muestra_completa <- function(n, phi, tipo_densidad, mw_idx = 1, gamma_gl = 0.05, me = 0, de = 1) {
  
  # 1. GENERACIÓN DE LA MUESTRA ÚNICA
  info_muestra <- generar_muestra(n = n, phi = phi, tipo = tipo_densidad, mw_idx = mw_idx, me = me, de = de)
  X <- info_muestra$X
  gsup <- info_muestra$gsup
  y_teo_func <- info_muestra$y_teo_func
  
  # 2. DEFINICIÓN DE LAS DOS MALLAS DE EVALUACIÓN
  # A) Malla Global
  lim_inf <- 0 ; lim_sup <- qlnorm(0.999999, 0, 0.5); n_grid <- 200
  grid_global <- seq(lim_inf, lim_sup, length.out = n_grid)
  delta_x <- grid_global[2] - grid_global[1]
  y_teo_global <- y_teo_func(grid_global)
  
  # B) Malla Puntual (6 puntos críticos)
  grid_puntual <- obtener_puntos_evaluacion(tipo = tipo_densidad, mw_idx = mw_idx, me = me, de = de)
  y_teo_puntual <- y_teo_func(grid_puntual)
  
  # 3. ESTIMACIÓN CLÁSICA (Malla Global nativa)
  d_nrd0 <- density(X, bw = "nrd0", from = lim_inf, to = lim_sup, n = n_grid)
  d_sj   <- density(X, bw = "SJ",   from = lim_inf, to = lim_sup, n = n_grid)
  d_ucv  <- tryCatch({ suppressWarnings(density(X, bw = "ucv", from = lim_inf, to = lim_sup, n = n_grid))
  }, error = function(e) return(d_nrd0))
  
  # 4. EXTRACCIÓN PUNTUAL CLÁSICA
  # Extraemos las ventanas (h) calculadas por las reglas clásicas
  h_nrd0 <- d_nrd0$bw
  h_ucv  <- d_ucv$bw
  h_sj   <- d_sj$bw
  
  # =====================================================================
  # ESTIMADOR POR NÚCLEO (Kernel Gaussiano)
  # =====================================================================
  estimador_nucleo_clasico <- function(X, puntos_eval, h) {
    # Para cada punto de evaluación, calcula el promedio del kernel sobre toda la muestra X
    sapply(puntos_eval, function(x_eval) {
      mean(dnorm(x_eval, mean = X, sd = h))
    })
  }
  
  # Calculamos la densidad analítica exacta en los 6 puntos
  dens_nrd0_pt <- estimador_nucleo_clasico(X, grid_puntual, h_nrd0)
  dens_ucv_pt  <- estimador_nucleo_clasico(X, grid_puntual, h_ucv)
  dens_sj_pt   <- estimador_nucleo_clasico(X, grid_puntual, h_sj)
  
  # 5. ESTIMACIÓN GL (Se aprovecha la matriz H para ambos)
  u <- seq(0, floor(log(n)) * (2/3), by = 0.1)
  H <- exp(-u)
  V_H <- Vh_vectorizado(H, n, gsup, gamma_gl, phi)
  
  gl_global <- GL_matrix_estimator(X, gridcal = grid_global, H, V_H)$densidad
  gl_puntual <- GL_matrix_estimator(X, gridcal = grid_puntual, H, V_H)$densidad
  
  # 6. CÁLCULO DE ERRORES 
  # A) Error Integrado (ISE)
  ise_nrd0 <- sum((d_nrd0$y - y_teo_global)^2) * delta_x
  ise_ucv  <- sum((d_ucv$y  - y_teo_global)^2) * delta_x
  ise_sj   <- sum((d_sj$y   - y_teo_global)^2) * delta_x
  ise_gl   <- sum((gl_global - y_teo_global)^2) * delta_x
  
  ISE_vector <- c(Silverman = ise_nrd0, UCV = ise_ucv, SJ = ise_sj, GL = ise_gl)
  
  # B) Error Puntual Cuadrático (SE)
  err_nrd0_pt <- (dens_nrd0_pt - y_teo_puntual)^2
  err_ucv_pt  <- (dens_ucv_pt  - y_teo_puntual)^2
  err_sj_pt   <- (dens_sj_pt   - y_teo_puntual)^2
  err_gl_pt   <- (gl_puntual   - y_teo_puntual)^2
  
  SE_matriz <- cbind(ROT = err_nrd0_pt, UCV = err_ucv_pt, SJ = err_sj_pt, GL = err_gl_pt)
  rownames(SE_matriz) <- names(grid_puntual)
  
  # Retorna TODO para esta muestra
  return(list(ISE = ISE_vector, SE_Puntual = SE_matriz, Muestra = X))
}


# =====================================================================
# ESTUDIO DE SIMULACIÓN MONTE CARLO (B = 1000)
# =====================================================================

ejecutar_montecarlo <- function(B = 1000, n = 120, phi = 0, tipo_densidad = "lognormal", mw_idx = 1, gamma_gl = 0.05, me = 0, de = 1) {
  
  # Crear el clúster dejando núcleos libres para el sistema operativo
  num_cores <- max(1, parallel::detectCores() - 2)
  cl <- makeCluster(num_cores)
  on.exit({ stopCluster(cl); closeAllConnections() }, add = TRUE)
  registerDoParallel(cl)
  
  cat(sprintf("Iniciando Monte Carlo  (%d iter). Densidad: %s\n", B, tipo_densidad))
  
  # Ejecutamos las B réplicas (.combine = 'list' para guardar resultados complejos)
  resultados_mc <- foreach(m = 1:B, .packages = c("stats", "matrixStats", "nor1mix"),
                           .export = c("evaluar_muestra_completa", "generar_muestra", 
                                       "obtener_puntos_evaluacion", "generar_uniformes_ar1", 
                                       "GL_matrix_estimator", "Vh_vectorizado", "QMix")) %dorng% {
                                         
                                         evaluar_muestra_completa(n = n, phi = phi, tipo_densidad = tipo_densidad, 
                                                                  mw_idx = mw_idx, gamma_gl = gamma_gl, me = me, de = de)
                                       }
  
  # 1. Extraer y promediar los ISE para obtener el MISE Global
  lista_ISE <- lapply(resultados_mc, function(res) res$ISE)
  matriz_ISE <- do.call(rbind, lista_ISE)
  MISE_final <- colMeans(matriz_ISE)
  
  # 2. Extraer y promediar las matrices SE para obtener el MSE Puntual
  lista_SE_puntual <- lapply(resultados_mc, function(res) res$SE_Puntual)
  # Usamos Reduce para sumar todas las matrices elemento a elemento, luego dividimos por B
  MSE_Puntual_final <- Reduce("+", lista_SE_puntual) / B
  
  # 3. Extraer y apilar las Muestras (Matriz B x n)
  lista_muestras <- lapply(resultados_mc, function(res) res$Muestra)
  matriz_muestras <- do.call(rbind, lista_muestras)
  
  # Retorna el panorama completo
  return(list(
    MISE_Global = MISE_final, 
    MSE_Puntual = MSE_Puntual_final,
    Muestras = matriz_muestras,
    Matriz_ISE_Cruda = matriz_ISE 
  ))
}

#############################################################################
# INICIO DE SIMULACION
############################################################################

# 1. Cargar los gammas precalculados ANTES de ejecutar los escenarios
load("Resultados_PRP_Lognormal_n120.RData")

# 2. Definir el diseño experimental
niveles_phi <- c(0, 0.5, 0.9)
distribuciones <- c("lognormal", "gamma")

# Función auxiliar para mapear el phi con el gamma cargado desde el .RData
obtener_gamma <- function(phi) {
  phi_str <- as.character(round(phi, 2))
  
  if(phi_str == "0")    return(Resultados_PRP_Lognormal_n120$gam.phi0)
  if(phi_str == "0.5")  return(Resultados_PRP_Lognormal_n120$gam.phi05)
  if(phi_str == "0.9")  return(Resultados_PRP_Lognormal_n120$gam.phi09)
  
  stop(paste("Error crítico: No se encontró un gamma calibrado para phi =", phi))
}

# 3. Ejecutar el bucle sobre la grilla de parámetros
resultados_lista <- list()     # LISTA PARA EL MISE
resultados_mse_lista <- list() # LISTA PARA EL MSE
resultados_muestras_lista <- list() # LISTA PARA GUARDAR LAS MATRICES DE DATOS
resultados_anova_lista <- list() # LISTA PARA EL ANOVA
fila <- 1

# Fijamos el tamaño de muestra
n_actual <- 120

for (phi_actual in niveles_phi) {
  gamma_actual <- obtener_gamma(phi_actual)
  
  for (dist_actual in distribuciones) {
    
    # Ejecutamos el Monte Carlo para la combinación actual
    res <- ejecutar_montecarlo(B = 1000, n = 120, phi = phi_actual, 
                               tipo_densidad = dist_actual, 
                               me = 3, de = 1, 
                               gamma_gl = gamma_actual)
    
    # A) EXTRACCIÓN Y GUARDADO DEL MISE GLOBAL
    mise_rot <- res$MISE_Global["Silverman"]
    mise_ucv <- res$MISE_Global["UCV"]
    mise_sj  <- res$MISE_Global["SJ"]
    mise_gl  <- res$MISE_Global["GL"]
    
    resultados_lista[[fila]] <- data.frame(
      Phi = phi_actual,
      Distribucion = tools::toTitleCase(dist_actual),
      ROT = mise_rot,
      UCV = mise_ucv,
      SJ = mise_sj,
      GL = mise_gl
    )
    
    # B) EXTRACCIÓN Y GUARDADO DEL MSE PUNTUAL
    # Extraemos la matriz 6x4 de esta iteración exacta
    matriz_mse <- res$MSE_Puntual
    
    # La convertimos en data.frame y le agregamos sus etiquetas
    df_mse_temp <- as.data.frame(matriz_mse)
    df_mse_temp$Punto <- rownames(matriz_mse)
    df_mse_temp$Phi <- phi_actual
    df_mse_temp$Distribucion <- tools::toTitleCase(dist_actual)
    
    # Reordenamos las columnas (Etiquetas primero, luego los métodos)
    df_mse_temp <- df_mse_temp[, c("Phi", "Distribucion", "Punto", "ROT", "UCV", "SJ", "GL")]
    
    # Guardamos en la nueva lista
    resultados_mse_lista[[fila]] <- df_mse_temp
    
    # C) GUARDADO DE LAS MUESTRAS
    # Guardamos la matriz de 1000 x n junto con sus metadatos
    resultados_muestras_lista[[fila]] <- list(
      n = n_actual,
      Phi = phi_actual,
      Distribucion = nombre_dist,
      Matriz_Datos = res$Muestras
    )
    
    # D) GUARDAR LA MATRIZ PARA EL ANOVA
    resultados_anova_lista[[fila]] <- list(
      Phi = phi_actual,
      Distribucion = nombre_dist,
      Matriz_Datos = res$Matriz_ISE_Cruda # Capturamos la matriz 1000x4
    )
    
    fila <- fila + 1
  }
}

# 4. Consolidar resultados en Data Frames finales
# Tabla 1: MISE
df_final_mise <- do.call(rbind, resultados_lista)
df_final_mise <- df_final_mise[order(df_final_mise$Distribucion, df_final_mise$Phi), ]
rownames(df_final_mise) <- NULL # Limpiar nombres de fila

# Tabla 2: MSE Puntual
df_final_mse <- do.call(rbind, resultados_mse_lista)
rownames(df_final_mse) <- NULL # Limpiar nombres de fila

# # =====================================================================
# # EXPORTACIÓN A HOJA DE CÁLCULO
# # =====================================================================
if (!require(writexl)) install.packages("writexl")
library(writexl)

setwd("~/Desktop/tesis_sebastian/adaptive-density-weak-dep/scripts/simulations/Resultados")

# 1. Guardar tablas en Excel
ruta_xlsx <- "mise_simulacion_lognorm_gamma_n120.xlsx"
write_xlsx(df_final_mise, path = ruta_xlsx)
cat("Archivo XLSX guardado en:", file.path(getwd(), ruta_xlsx), "\n")

ruta_xlsx2 <- "mse_simulacion_lognorm_gamma_n120.xlsx"
write_xlsx(df_final_mse, path = ruta_xlsx2)
cat("Archivo XLSX guardado en:", file.path(getwd(), ruta_xlsx2), "\n")

# 2. Guardar las muestras generadas en RData
ruta_muestras <- "muestras_generadas_lognorm_gamma_n120.RData"
save(resultados_muestras_lista, file = ruta_muestras)

cat("Archivos de resultados (XLSX) y Muestras (RData) guardados exitosamente en:\n", getwd(), "\n")

# =====================================================================
# GUARDADO DE MATRICES PARA ANOVA
# =====================================================================
ruta_anova <- "matrices_ise_anova_lognorm_gamma_n120.RData"
save(resultados_anova_lista, file = ruta_anova)
cat("Matrices para ANOVA guardadas en:", file.path(getwd(), ruta_anova), "\n")

# # =====================================================================
# # 6. GENERACIÓN DEL CÓDIGO LATEX (FORMATO APA 7)
# # =====================================================================
# 
# # Formateamos los números a 4 decimales
# df_final[, 3:6] <- lapply(df_final[, 3:6], function(x) sprintf("%.4f", x))
# 
# # Renombrar columnas para la tabla
# colnames(df_final) <- c("$\\phi$", "Distribución", "ROT", "UCV", "SJ", "GL")
# 
# # Crear el objeto xtable con las directrices APA 7
# tabla_latex <- xtable(
#   df_final, 
#   caption = "Comparación del MISE entre selectores de ancho de ventana según distribución y nivel de dependencia AR(1)",
#   label = "tab:mise_comparacion",
#   align = c("l", "c", "l", "c", "c", "c", "c") # Alineación de columnas
# )
# 
# # Imprimir en consola el código LaTeX listo para copiar y pegar
# print(tabla_latex, 
#       include.rownames = FALSE, 
#       sanitize.colnames.function = identity, # Permite que LaTeX procese el símbolo phi
#       booktabs = TRUE,                       # Usa toprule, midrule y bottomrule (Exigencia APA 7)
#       caption.placement = "top",             # APA 7 exige el título arriba
#       table.placement = "htbp")
