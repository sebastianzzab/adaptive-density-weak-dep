###########################################################
# Parámetros Globales y Bibliotecas
###########################################################
# install.packages("nor1mix")
library(nor1mix)

n <- 1000          # Tamaño de la muestra
phi <- 0.9         # Autocorrelación AR(1)
gamma_fijo <- 0.19 # Gamma original del artículo para baseline

mw_list <- list(
  MW.nm1,  MW.nm2,  MW.nm3,  MW.nm4,  MW.nm5,
  MW.nm6,  MW.nm7,  MW.nm8,  MW.nm9,  MW.nm10,
  MW.nm11, MW.nm12, MW.nm13, MW.nm14, MW.nm15
)

###########################################################
# Funciones de Generación de Datos Dependientes
###########################################################
generar_uniformes_ar1 <- function(nsim, phi) {
  if (phi == 0) {
    Z <- rnorm(nsim)
  } else {
    Z <- arima.sim(list(ar=phi), n=nsim) 
  }
  sd_Z <- if (phi == 0) 1 else sqrt(1 / (1 - phi^2))
  return(pnorm(Z, mean=0, sd=sd_Z))
}

generar_dependiente_mw <- function(nsim, phi, mw_obj) {
  U <- generar_uniformes_ar1(nsim, phi)
  return(qnorMix(U, mw_obj))
}

###########################################################
# Funciones Originales GL (Adaptadas para usar Gamma Dinámico)
###########################################################
gh <- function(x, h, X) {
  n_len <- length(X)
  return(sum(dnorm((x-X)/h))/(n_len*h))
}

ghh <- function(x, h, hp, X) {
  hs <- sqrt(h^2 + hp^2)
  return(gh(x, hs, X))
}

Vh <- function(h, n, gsup, gamma) {
  delta.n <- sqrt(log(n))
  K1 <- 1
  K2 <- sqrt(1/(2*sqrt(pi)))
  return(sqrt(2*gamma*gsup)*K2*(K1+1)*(1+delta.n)*(((log(n))^(-1/2))/sqrt(n*h)))
}

Ahx <- function(x, h, H, X, n, gsup, gamma) {
  valores <- numeric(length(H))
  for(i in seq_along(H)) {
    hp <- H[i]
    valores[i] <- max(abs(ghh(x, h, hp, X) - gh(x, hp, X)) - Vh(hp, n, gsup, gamma), 0)
  }
  return(max(valores))
}

hx <- function(x, H, X, n, gsup, gamma) {
  criterio <- numeric(length(H))
  for(i in seq_along(H)) {
    h <- H[i]
    criterio[i] <- Ahx(x, h, H, X, n, gsup, gamma) + Vh(h, n, gsup, gamma)
  }
  hhat <- H[which.min(criterio)]
  return(list(hhat=hhat, ghat=gh(x, hhat, X)))
}

GLdens <- function(gridcal, H, X, n, gsup, gamma) {
  m <- length(gridcal)
  ghat <- numeric(m)
  for(i in 1:m) {
    res <- hx(gridcal[i], H, X, n, gsup, gamma)
    ghat[i] <- res$ghat
  }
  return(ghat)
}

###########################################################
# Función: Calibración de Gamma usando PCO
###########################################################
Obtener_Gamma_PCO <- function(H, h_min, X, gridcal, Delta) {
  n_len <- length(X)
  gh_min <- sapply(gridcal, function(z) gh(z, h_min, X))
  
  norm_diff_vec <- numeric(length(H))
  for(i in seq_along(H)) {
    gh_h <- sapply(gridcal, function(z) gh(z, H[i], X))
    norm_diff_vec[i] <- sum((gh_h - gh_min)^2) * Delta
  }
  
  gammas_grid <- seq(-2, 3, by=0.05)
  h_selected <- numeric(length(gammas_grid))
  norm2_hmin <- 1 / (2 * sqrt(pi) * h_min)
  
  for(j in seq_along(gammas_grid)) {
    gam <- gammas_grid[j]
    criterios_gam <- numeric(length(H))
    for(i in seq_along(H)) {
      h <- H[i]
      norm2_h <- 1 / (2 * sqrt(pi) * h)
      inner_k <- 1 / sqrt(2 * pi * (h^2 + h_min^2))
      norm2_diff_K <- norm2_h + norm2_hmin - 2 * inner_k
      
      pen_gam <- (gam * norm2_h - norm2_diff_K) / n_len
      criterios_gam[i] <- norm_diff_vec[i] + pen_gam
    }
    h_selected[j] <- H[which.min(criterios_gam)]
  }
  
  saltos <- diff(h_selected)
  gamma_min_hat <- gammas_grid[which.max(saltos)]
  return(gamma_min_hat + 1) # Regla óptima teórica PCO
}

###########################################################
# Bucle de Evaluación para las 15 Densidades MW
###########################################################
resultados_ise <- data.frame(MW=integer(), GL_Fijo=numeric(), GL_PCO=numeric(), CV=numeric(), Silverman=numeric())
u <- seq(0, floor(log(n))*(2/3), by = 0.1)
H <- sort(exp(-u))
h_min <- min(H)

cat("Iniciando evaluación de densidades...\n")

for(idx in 1:15) {
  cat(sprintf("Procesando Densidad MW %d...\n", idx))
  
  mw_obj <- mw_list[[idx]]
  X <- generar_dependiente_mw(n, phi, mw_obj)
  
  grid_min <- qnorMix(0.001, mw_obj)
  grid_max <- qnorMix(0.999, mw_obj)
  gridcal <- seq(grid_min, grid_max, length=100) # Reducido a 100 para viabilidad computacional de GL
  Delta <- diff(gridcal)[1]
  g.real <- dnorMix(gridcal, mw_obj)
  
  # Estimar Cota Superior gsup empírica exacta (f_max)
  gsup <- max(g.real)
  
  # 1. PCO Gamma Selection
  gamma_opt <- Obtener_Gamma_PCO(H, h_min, X, gridcal, Delta)
  
  # 2. GL Estimations
  ghat_gl_fijo <- GLdens(gridcal, H, X, n, gsup, gamma_fijo)
  ghat_gl_pco <- GLdens(gridcal, H, X, n, gsup, gamma_opt)
  
  # 3. CV & Silverman
  hcv <- bw.ucv(X)
  gh_cv <- sapply(gridcal, function(z) gh(z, hcv, X))
  
  h.sil <- bw.nrd0(X)
  gh_sil <- sapply(gridcal, function(z) gh(z, h.sil, X))
  
  # 4. ISE Calculation
  ise_gl_fijo <- sum((ghat_gl_fijo - g.real)^2) * Delta
  ise_gl_pco  <- sum((ghat_gl_pco - g.real)^2) * Delta
  ise_cv      <- sum((gh_cv - g.real)^2) * Delta
  ise_sil     <- sum((gh_sil - g.real)^2) * Delta
  
  resultados_ise <- rbind(resultados_ise, data.frame(
    MW = idx, GL_Fijo = ise_gl_fijo, GL_PCO = ise_gl_pco, 
    CV = ise_cv, Silverman = ise_sil
  ))
}

###########################################################
# Gráfico Diagnóstico Global: Comparación de ISE
###########################################################
print(resultados_ise)

# Barplot comparativo GL Fijo vs GL PCO
library(graphics)
matriz_barras <- t(as.matrix(resultados_ise[, c("GL_Fijo", "GL_PCO")]))
colnames(matriz_barras) <- paste("MW", 1:15)

barplot(matriz_barras, beside = TRUE, col = c("darkgray", "darkgreen"),
        legend.text = c("GL Fijo (0.19)", "GL PCO-Calibrado"),
        args.legend = list(x = "topleft", bty = "n"),
        main = "Comparativa de ISE: GL Clásico vs GL PCO-Adaptativo",
        ylab = "Integrated Squared Error (ISE)",
        las = 2, cex.names = 0.8)