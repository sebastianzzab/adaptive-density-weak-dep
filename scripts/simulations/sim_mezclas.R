###############################################################################
# Script: Estudio de Simulación Monte Carlo - Grupo 2 (Normal y Laplace)
# Autor: Sebastian Zabala
# Reproducibilidad: Total (doRNG, Semillas Dinámicas, Entorno Aislado)
###############################################################################

rm(list = ls()) # Limpiar entorno

# =====================================================================
# 1. CONFIGURACIÓN GLOBAL Y SEMILLAS (GRUPO 3)
# =====================================================================
# CONFIG_n <- 60
# CONFIG_semilla_maestra <- 3060 # Semillas en la serie 3000

# CONFIG_n <- 120
# CONFIG_semilla_maestra <- 3120

CONFIG_n <- 250
CONFIG_semilla_maestra <- 3250

CONFIG_B <- 1000
CONFIG_niveles_phi <- c(0.3, 0.5, 0.9)
CONFIG_distribuciones <- c("mezcla", "mw")
CONFIG_ruta_salida <- "./resultados_simulacion/Grupo3_Mezclas"

# =====================================================================
# 2. TABLA INTERNA DE GAMMAS ÓPTIMOS (GRUPO 2)
# =====================================================================
tabla_gammas <- data.frame(
  n = c(rep(60, 3), rep(120, 3), rep(250, 3)),
  phi = rep(c(0.3, 0.5, 0.9), times = 3),
  gamma_optimo = c(
    0.48125,	0.41875,	0.4625, # Valores para Normal/Laplace n=60
    0.15,	0.475,	0.39375, # Valores para Normal/Laplace n=120
    0.1,	0.0975,	0.43125  # Valores para Normal/Laplace n=250
  )
)
		

obtener_gamma <- function(n_req, phi_req) {
  val <- tabla_gammas$gamma_optimo[tabla_gammas$n == n_req & tabla_gammas$phi == phi_req]
  if(length(val) == 0 || is.na(val)) stop("Gamma no encontrado para n=", n_req, " y phi=", phi_req)
  return(val)
}

# =====================================================================
# 3. INSTALACIÓN Y CARGA DE DEPENDENCIAS
# =====================================================================
paquetes <- c("matrixStats", "nor1mix", "doParallel", "foreach", "doRNG", "xtable", "tidyr", "writexl")
nuevos_paquetes <- paquetes[!(paquetes %in% installed.packages()[,"Package"])]
if(length(nuevos_paquetes)) install.packages(nuevos_paquetes)
invisible(lapply(paquetes, library, character.only = TRUE))

# =====================================================================
# 4. DEFINICIÓN DE DENSIDADES Y GENERADORES 
# =====================================================================
mw_list <- list(MW.nm1, MW.nm2, MW.nm3, MW.nm4, MW.nm5, MW.nm6, MW.nm7, MW.nm8, 
                MW.nm9, MW.nm10, MW.nm11, MW.nm12, MW.nm13, MW.nm14, MW.nm15)

GMix <- function(x) { 0.5 * pnorm(x, mean = -2, sd = 1) + 0.5 * pnorm(x, mean = 2, sd = 1) }
grid.mix.x <- seq(-8, 8, length.out = 10000)
grid.mix.u <- GMix(grid.mix.x)
QMix <- function(u) { approx(x = grid.mix.u, y = grid.mix.x, xout = u, rule = 2)$y }

generar_uniformes_ar1 <- function(nsim, phi, me = 0, de = 1) {
  if (phi == 0) { Z <- rnorm(nsim) } else { Z <- arima.sim(list(ar = phi), n = nsim) }
  sd_Z <- if (phi == 0) 1 else sqrt(1 / (1 - phi^2))
  return(pnorm(Z, mean = 0, sd = sd_Z))
}

generar_truncnorm_dep <- function(n, phi, a, b, m = 0, de = 1) {
  U <- generar_uniformes_ar1(nsim = n, phi = phi)
  F_a <- pnorm(a, mean = m, sd = de); F_b <- pnorm(b, mean = m, sd = de)
  return(qnorm((F_b - F_a) * U + F_a, mean = m, sd = de))
}

generar_muestra <- function(n, phi, tipo, mw_idx = 1, a = -2, b = 2, me = 0, de = 1) {
  U <- generar_uniformes_ar1(nsim = n, phi = phi, me, de)
  if (tipo == "truncnorm") {
    X <- generar_truncnorm_dep(n, phi, a, b, 0, 1); c_norm <- pnorm(b) - pnorm(a); gsup <- dnorm(0) / c_norm 
    y_teo_func <- function(x) { y <- dnorm(x) / c_norm; y[x < a | x > b] <- 0; return(y) }
  } else if (tipo == "normal-est") {
    X <- qnorm(U); gsup <- 1 / sqrt(2 * pi); y_teo_func <- dnorm; lim_inf <- -8; lim_sup <- 8
  } else if (tipo == "normal") {
    X <- qnorm(U, mean = me, sd = de); gsup <- 1 / (sqrt(2 * pi) * de); y_teo_func <- function(x) dnorm(x, mean = me, sd = de); lim_inf <- -8; lim_sup <- 8
  } else if (tipo == "lognormal") {
    X <- qlnorm(U, meanlog = 0, sdlog = 0.5); gsup <- exp(0.125) / (0.5 * sqrt(2 * pi)); y_teo_func <- function(x) dlnorm(x, meanlog = 0, sdlog = 0.5); lim_inf <- 0; lim_sup <- qlnorm(0.999999, 0, 0.5)
  } else if (tipo == "mezcla") {
    X <- QMix(U); gsup <- 0.5 * dnorm(-2, -2, 1) + 0.5 * dnorm(-2, 2, 1); y_teo_func <- function(x) 0.5 * dnorm(x, -2, 1) + 0.5 * dnorm(x, 2, 1); lim_inf <- -8; lim_sup <- 8
  } else if (tipo == "mw") {
    mw_obj <- mw_list[[mw_idx]]; X <- qnorMix(U, mw_obj)
    opt <- optimize(function(x) dnorMix(x, mw_obj), interval = c(-5, 5), maximum = TRUE)
    gsup <- opt$objective; y_teo_func <- function(x) dnorMix(x, mw_obj); lim_inf <- -8; lim_sup <- 8
  } else if (tipo == "logistica") {
    s <- sqrt(3) / pi; X <- qlogis(U, location = 3, scale = s); gsup <- dlogis(3, location = 3, scale = s); y_teo_func <- function(x) dlogis(x, location = 3, scale = s); lim_inf <- -8; lim_sup <- 8
  } else if (tipo == "laplace") {
    b_param <- 1 / sqrt(2); U_centrada <- U - 0.5; X <- 3 - b_param * sign(U_centrada) * log(1 - 2 * abs(U_centrada)); gsup <- 1 / (2 * b_param); y_teo_func <- function(x) (1 / (2 * b_param)) * exp(-abs(x - 3) / b_param); lim_inf <- -8; lim_sup <- 8
  } else if (tipo == "gamma") {
    media_ln <- exp(0.125); var_ln <- (exp(0.25) - 1) * exp(0.25); shape_g <- (media_ln^2) / var_ln; rate_g <- media_ln / var_ln
    X <- qgamma(U, shape = shape_g, rate = rate_g); moda_g <- (shape_g - 1) / rate_g
    gsup <- dgamma(moda_g, shape = shape_g, rate = rate_g); y_teo_func <- function(x) dgamma(x, shape = shape_g, rate = rate_g); lim_inf <- 0; lim_sup <- 8
  } else if (tipo == "weibull") {
    shape_w <- 1.5; media_ln <- exp(0.125); scale_w <- media_ln / gamma(1 + 1/shape_w); X <- qweibull(U, shape = shape_w, scale = scale_w)
    moda_w <- scale_w * ((shape_w - 1)/shape_w)^(1/shape_w); gsup <- dweibull(moda_w, shape = shape_w, scale = scale_w); y_teo_func <- function(x) dweibull(x, shape = shape_w, scale = scale_w); lim_inf <- 0; lim_sup <- 8
  } else { stop("Tipo de densidad no válido.") }
  return(list(X = X, gsup = gsup, y_teo_func = y_teo_func, lim_inf = lim_inf, lim_sup = lim_sup))
}

obtener_puntos_evaluacion <- function(tipo, mw_idx = 1, me = 0, de = 1) {
  nombres <- c("Moda", "P05", "Q1", "Mediana", "Q3", "P95"); pts <- numeric(6)
  if (tipo == "normal") { pts[1] <- me; pts[2:6] <- qnorm(c(0.05, 0.25, 0.50, 0.75, 0.95), mean = me, sd = de)
  } else if (tipo == "lognormal") { pts[1] <- exp(0 - 0.5^2); pts[2:6] <- qlnorm(c(0.05, 0.25, 0.50, 0.75, 0.95), meanlog = 0, sdlog = 0.5)
  } else if (tipo == "logistica") { s <- sqrt(3) / pi; pts[1] <- 3; pts[2:6] <- qlogis(c(0.05, 0.25, 0.50, 0.75, 0.95), location = 3, scale = s)
  } else if (tipo == "laplace") { b_param <- 1 / sqrt(2); pts[1] <- 3; q_laplace <- function(p, m, b) { ifelse(p < 0.5, m + b * log(2 * p), m - b * log(2 * (1 - p))) }; pts[2:6] <- q_laplace(c(0.05, 0.25, 0.50, 0.75, 0.95), m = 3, b = b_param)
  } else if (tipo == "gamma") { media_ln <- exp(0.125); var_ln <- (exp(0.25) - 1) * exp(0.25); shape_g <- (media_ln^2) / var_ln; rate_g <- media_ln / var_ln; pts[1] <- (shape_g - 1) / rate_g; pts[2:6] <- qgamma(c(0.05, 0.25, 0.50, 0.75, 0.95), shape = shape_g, rate = rate_g)
  } else if (tipo == "weibull") { shape_w <- 1.5; media_ln <- exp(0.125); scale_w <- media_ln / gamma(1 + 1/shape_w); pts[1] <- scale_w * ((shape_w - 1)/shape_w)^(1/shape_w); pts[2:6] <- qweibull(c(0.05, 0.25, 0.50, 0.75, 0.95), shape = shape_w, scale = scale_w)
  } else if (tipo == "mezcla") { pts[1] <- 2; pts[2:6] <- QMix(c(0.05, 0.25, 0.50, 0.75, 0.95))
  } else if (tipo == "mw") { mw_obj <- mw_list[[mw_idx]]; opt <- optimize(function(x) dnorMix(x, mw_obj), interval = c(-5, 5), maximum = TRUE); pts[1] <- opt$maximum; pts[2:6] <- qnorMix(c(0.05, 0.25, 0.50, 0.75, 0.95), mw_obj) }
  names(pts) <- nombres; return(pts)
}

# =====================================================================
# 5. ESTIMADOR GL Y EVALUACIÓN
# =====================================================================
Vh_vectorizado <- function(H_vec, n, gsup, gamma, phi) {
  delta.n <- sqrt(log(n)); K2 <- sqrt(1 / (2 * sqrt(pi)))
  return(sqrt(2 * gamma * gsup) * K2 * (2) * (1 + delta.n) * ((log(n)^(-1/2)) / sqrt(n * H_vec)))
}

GL_matrix_estimator <- function(X, gridcal, H, V_H) {
  n <- length(X); dist_base <- outer(gridcal, X, "-")
  Matriz_GH <- matrix(0, nrow = length(gridcal), ncol = length(H))
  for(j in seq_along(H)) Matriz_GH[, j] <- rowSums(dnorm(dist_base / H[j])) / (n * H[j])
  
  Criterio_A_Matriz <- matrix(0, nrow = length(gridcal), ncol = length(H))
  for(i in seq_along(H)) {
    GHH_temp <- matrix(0, nrow = length(gridcal), ncol = length(H))
    for(j in seq_along(H)) { hs <- sqrt(H[i]^2 + H[j]^2); GHH_temp[, j] <- rowSums(dnorm(dist_base / hs)) / (n * hs) }
    Penalizado <- sweep(abs(GHH_temp - Matriz_GH), MARGIN = 2, STATS = V_H, FUN = "-")
    Penalizado[Penalizado < 0] <- 0; Criterio_A_Matriz[, i] <- rowMaxs(Penalizado)
  }
  Criterio_Final <- sweep(Criterio_A_Matriz, MARGIN = 2, STATS = V_H, FUN = "+")
  indices_minimos <- max.col(-Criterio_Final, ties.method = "first")
  return(list(densidad = Matriz_GH[cbind(1:length(gridcal), indices_minimos)], h_local = H[indices_minimos]))
}

evaluar_muestra_completa <- function(n, phi, tipo_densidad, mw_idx = 1, gamma_gl = 0.05, me = 0, de = 1) {
  info_muestra <- generar_muestra(n = n, phi = phi, tipo = tipo_densidad, mw_idx = mw_idx, me = me, de = de)
  X <- info_muestra$X; gsup <- info_muestra$gsup; y_teo_func <- info_muestra$y_teo_func
  lim_inf <- info_muestra$lim_inf; lim_sup <- info_muestra$lim_sup
  
  grid_global <- seq(lim_inf, lim_sup, length.out = 200); delta_x <- grid_global[2] - grid_global[1]
  y_teo_global <- y_teo_func(grid_global)
  grid_puntual <- obtener_puntos_evaluacion(tipo = tipo_densidad, mw_idx = mw_idx, me = me, de = de)
  y_teo_puntual <- y_teo_func(grid_puntual)
  
  d_nrd0 <- density(X, bw = "nrd0", from = lim_inf, to = lim_sup, n = 200)
  d_sj   <- density(X, bw = "SJ",   from = lim_inf, to = lim_sup, n = 200)
  d_ucv  <- tryCatch({ suppressWarnings(density(X, bw = "ucv", from = lim_inf, to = lim_sup, n = 200)) }, error = function(e) return(d_nrd0))
  
  estimador_nucleo_clasico <- function(X, puntos_eval, h) { sapply(puntos_eval, function(x_eval) mean(dnorm(x_eval, mean = X, sd = h))) }
  dens_nrd0_pt <- estimador_nucleo_clasico(X, grid_puntual, d_nrd0$bw)
  dens_ucv_pt  <- estimador_nucleo_clasico(X, grid_puntual, d_ucv$bw)
  dens_sj_pt   <- estimador_nucleo_clasico(X, grid_puntual, d_sj$bw)
  
  H <- exp(-seq(0, floor(log(n)) * (2/3), by = 0.1))
  V_H <- Vh_vectorizado(H, n, gsup, gamma_gl, phi)
  
  gl_global <- GL_matrix_estimator(X, gridcal = grid_global, H, V_H)$densidad
  gl_puntual <- GL_matrix_estimator(X, gridcal = grid_puntual, H, V_H)$densidad
  
  ise_nrd0 <- sum((d_nrd0$y - y_teo_global)^2) * delta_x; ise_ucv  <- sum((d_ucv$y  - y_teo_global)^2) * delta_x
  ise_sj   <- sum((d_sj$y   - y_teo_global)^2) * delta_x; ise_gl   <- sum((gl_global - y_teo_global)^2) * delta_x
  ISE_vector <- c(Silverman = ise_nrd0, UCV = ise_ucv, SJ = ise_sj, GL = ise_gl)
  
  SE_matriz <- cbind(ROT = (dens_nrd0_pt - y_teo_puntual)^2, UCV = (dens_ucv_pt - y_teo_puntual)^2, SJ = (dens_sj_pt - y_teo_puntual)^2, GL = (gl_puntual - y_teo_puntual)^2)
  rownames(SE_matriz) <- names(grid_puntual)
  
  return(list(ISE = ISE_vector, SE_Puntual = SE_matriz, Muestra = X))
}

# =====================================================================
# 6. MOTOR MONTE CARLO PARALELIZADO
# =====================================================================
ejecutar_montecarlo <- function(B, n, phi, tipo_densidad, mw_idx, gamma_gl, me, de, semilla) {
  cl <- makeCluster(max(1, parallel::detectCores() - 2))
  on.exit({ stopCluster(cl); closeAllConnections() }, add = TRUE)
  registerDoParallel(cl)
  
  cat(sprintf("-> Ejecutando %s (n=%d, phi=%.1f) | Semilla: %d\n", tools::toTitleCase(tipo_densidad), n, phi, semilla))
  
  set.seed(semilla) # Fijación determinista
  resultados_mc <- foreach(m = 1:B, .packages = c("stats", "matrixStats", "nor1mix"),
                           .export = c("evaluar_muestra_completa", "generar_muestra", "obtener_puntos_evaluacion", 
                                       "generar_uniformes_ar1", "GL_matrix_estimator", "Vh_vectorizado", "QMix",
                                       "grid.mix.x", "grid.mix.u", "mw_list", "generar_truncnorm_dep")) %dorng% {
                                         evaluar_muestra_completa(n, phi, tipo_densidad, mw_idx, gamma_gl, me, de)
                                       }
  
  matriz_ISE <- do.call(rbind, lapply(resultados_mc, function(res) res$ISE))
  MSE_Puntual <- Reduce("+", lapply(resultados_mc, function(res) res$SE_Puntual)) / B
  matriz_muestras <- do.call(rbind, lapply(resultados_mc, function(res) res$Muestra))
  
  return(list(MISE_Global = colMeans(matriz_ISE), MSE_Puntual = MSE_Puntual, Muestras = matriz_muestras, Matriz_ISE_Cruda = matriz_ISE))
}

# =====================================================================
# 7. EJECUCIÓN PRINCIPAL
# =====================================================================
resultados_lista <- list(); resultados_mse_lista <- list()
resultados_muestras_lista <- list(); resultados_anova_lista <- list()
fila <- 1

for (phi_actual in CONFIG_niveles_phi) {
  gamma_actual <- obtener_gamma(CONFIG_n, phi_actual)
  
  for (dist_actual in CONFIG_distribuciones) {
    mw_idx_actual <- ifelse(dist_actual == "mw", 10, 1)
    
    # Nomenclatura limpia
    if(dist_actual == "mezcla") { nombre_dist <- "Mezcla Bimodal"
    } else if(dist_actual == "mw") { nombre_dist <- "Modelo 10 (Garra)"
    } else { nombre_dist <- tools::toTitleCase(dist_actual) }
    
    # Semilla única indexada matemáticamente
    semilla_escenario <- CONFIG_semilla_maestra + (phi_actual * 10) + which(CONFIG_distribuciones == dist_actual)
    
    res <- ejecutar_montecarlo(B = CONFIG_B, n = CONFIG_n, phi = phi_actual, tipo_densidad = dist_actual, 
                               mw_idx = mw_idx_actual, me = 3, de = 1, gamma_gl = gamma_actual, semilla = semilla_escenario)
    
    resultados_lista[[fila]] <- data.frame(n = CONFIG_n, Phi = phi_actual, Distribucion = nombre_dist, 
                                           ROT = res$MISE_Global["Silverman"], UCV = res$MISE_Global["UCV"], SJ = res$MISE_Global["SJ"], GL = res$MISE_Global["GL"])
    
    df_mse_temp <- as.data.frame(res$MSE_Puntual)
    df_mse_temp$Punto <- rownames(res$MSE_Puntual); df_mse_temp$n <- CONFIG_n; df_mse_temp$Phi <- phi_actual; df_mse_temp$Distribucion <- nombre_dist
    resultados_mse_lista[[fila]] <- df_mse_temp[, c("n", "Phi", "Distribucion", "Punto", "ROT", "UCV", "SJ", "GL")]
    
    resultados_muestras_lista[[fila]] <- list(n = CONFIG_n, Phi = phi_actual, Distribucion = nombre_dist, Matriz_Datos = res$Muestras)
    resultados_anova_lista[[fila]] <- list(n = CONFIG_n, Phi = phi_actual, Distribucion = nombre_dist, Matriz_Datos = res$Matriz_ISE_Cruda)
    
    fila <- fila + 1
  }
}

df_final_mise <- do.call(rbind, resultados_lista)
df_final_mse <- do.call(rbind, resultados_mse_lista)

# =====================================================================
# 8. EXPORTACIÓN DINÁMICA Y REGISTRO
# =====================================================================
sufijo <- sprintf("%s_n%d", paste(CONFIG_distribuciones, collapse="_"), CONFIG_n)

# Nombres correctos (Corregido el error de sobrescritura de Excel)
write_xlsx(df_final_mise, path = file.path(CONFIG_ruta_salida, paste0("mise_", sufijo, ".xlsx")))
write_xlsx(df_final_mse, path = file.path(CONFIG_ruta_salida, paste0("mse_puntual_", sufijo, ".xlsx")))
save(resultados_muestras_lista, file = file.path(CONFIG_ruta_salida, paste0("muestras_", sufijo, ".RData")))
save(resultados_anova_lista, file = file.path(CONFIG_ruta_salida, paste0("anova_", sufijo, ".RData")))

capture.output(sessionInfo(), file = file.path(CONFIG_ruta_salida, paste0("session_info_", sufijo, ".txt")))
cat("\nSimulación finalizada. Resultados en:", normalizePath(CONFIG_ruta_salida), "\n")