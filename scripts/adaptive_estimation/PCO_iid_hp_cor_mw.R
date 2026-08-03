###########################################################
#                  Configuración de Parámetros
###########################################################
# install.packages("nor1mix") # Descomentar si no está instalado
library(nor1mix)

n <- 1000          
phi <- 0.0         # Cambiar a 0.0 para probar i.i.d., o 0.9 para AR(1)
mw_index <- 1      

###########################################################
# Configuración de las Densidades de Marron & Wand
###########################################################
mw_list <- list(
  MW.nm1,  MW.nm2,  MW.nm3,  MW.nm4,  MW.nm5,
  MW.nm6,  MW.nm7,  MW.nm8,  MW.nm9,  MW.nm10,
  MW.nm11, MW.nm12, MW.nm13, MW.nm14, MW.nm15
)

###########################################################
# FUNCIONES PARA GENERAR VARIABLES DEPENDIENTES
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
# Generación de la Muestra Objetivo
###########################################################
mw_obj_selected <- mw_list[[mw_index]]
X <- generar_dependiente_mw(n, phi, mw_obj_selected)

grid_min <- qnorMix(0.001, mw_obj_selected)
grid_max <- qnorMix(0.999, mw_obj_selected)
gridcal <- seq(grid_min, grid_max, length=200)
Delta <- diff(gridcal)[1]

g.real <- dnorMix(gridcal, mw_obj_selected)

###########################################################
# Estimación de Densidad por PCO - Diagnóstico y Calibración
###########################################################

u <- seq(0, floor(log(n))*(2/3), by = 0.1)
H <- sort(exp(-u)) 
h_min <- min(H)

gh <- function(x, h, X) {
  n_len <- length(X)
  return(sum(dnorm((x-X)/h))/(n_len*h))
}

gh_min <- sapply(gridcal, function(z) gh(z, h_min, X))

PCO_Diagnostic_Selection <- function(H, h_min, X, gridcal, Delta, gh_min, g.real) {
  n_len <- length(X)
  
  # 1. Precomputar estimadores, distancias empíricas e ISE (Eficiencia Computacional)
  norm_diff_vec <- numeric(length(H))
  estimadores <- matrix(0, nrow=length(gridcal), ncol=length(H))
  ISE_vec <- numeric(length(H))
  
  for(i in seq_along(H)) {
    h <- H[i]
    gh_h <- sapply(gridcal, function(z) gh(z, h, X))
    estimadores[, i] <- gh_h
    norm_diff_vec[i] <- sum((gh_h - gh_min)^2) * Delta
    ISE_vec[i] <- sum((gh_h - g.real)^2) * Delta
  }
  
  # 2. Grid de lambdas centrado (abarca valores negativos para capturar el salto)
  lambdas_grid <- seq(-2, 3, by=0.05)
  h_selected <- numeric(length(lambdas_grid))
  ise_selected <- numeric(length(lambdas_grid))
  
  norm2_hmin <- 1 / (2 * sqrt(pi) * h_min)
  
  # 3. Evaluar el criterio PCO para cada lambda
  for(j in seq_along(lambdas_grid)) {
    lam <- lambdas_grid[j]
    criterios_lam <- numeric(length(H))
    
    for(i in seq_along(H)) {
      h <- H[i]
      norm2_h <- 1 / (2 * sqrt(pi) * h)
      inner_k <- 1 / sqrt(2 * pi * (h^2 + h_min^2))
      norm2_diff_K <- norm2_h + norm2_hmin - 2 * inner_k
      
      pen_lam <- (lam * norm2_h - norm2_diff_K) / n_len
      criterios_lam[i] <- norm_diff_vec[i] + pen_lam
    }
    
    idx_opt <- which.min(criterios_lam)
    h_selected[j] <- H[idx_opt]
    ise_selected[j] <- ISE_vec[idx_opt]
  }
  
  # 4. Encontrar el salto crítico de la heurística
  saltos <- diff(h_selected)
  idx_salto <- which.max(saltos)
  lambda_min_hat <- lambdas_grid[idx_salto]
  
  # 5. CORRECCIÓN TEÓRICA: La penalidad óptima añade una varianza (+1) a la mínima
  lambda_opt <- lambda_min_hat + 1
  
  # 6. Selección final usando lambda_opt
  criterios_final <- numeric(length(H))
  for(i in seq_along(H)) {
    h <- H[i]
    norm2_h <- 1 / (2 * sqrt(pi) * h)
    inner_k <- 1 / sqrt(2 * pi * (h^2 + h_min^2))
    norm2_diff_K <- norm2_h + norm2_hmin - 2 * inner_k
    
    pen_lam_opt <- (lambda_opt * norm2_h - norm2_diff_K) / n_len
    criterios_final[i] <- norm_diff_vec[i] + pen_lam_opt
  }
  idx_final <- which.min(criterios_final)
  
  # 7. Generar el Plot Diagnóstico de la Transición de Fase
  par(mfrow=c(2,1), mar=c(4,4,2,1))
  
  # Plot A: Selección de Ventana vs Lambda
  plot(lambdas_grid, h_selected, type="s", col="blue", lwd=2,
       xlab=expression(lambda), ylab=expression(hat(h)),
       main="Transición de Fase: Selección de Ventana")
  abline(v = lambda_min_hat, col="red", lty=2)
  abline(v = lambda_opt, col="darkgreen", lty=2)
  legend("topleft", legend=c("Salto (min)", "Óptimo (+1)"), 
         col=c("red", "darkgreen"), lty=2, bty="n")
  
  # Plot B: ISE vs Lambda
  plot(lambdas_grid, ise_selected, type="l", col="purple", lwd=2,
       xlab=expression(lambda), ylab="ISE",
       main="Error Cuadrático Integrado (ISE) vs Lambda")
  abline(v = lambda_min_hat, col="red", lty=2)
  abline(v = lambda_opt, col="darkgreen", lty=2)
  
  par(mfrow=c(1,1)) # Resetear layout
  
  return(list(
    hhat = H[idx_final],
    ghat = estimadores[, idx_final],
    lambda_min = lambda_min_hat,
    lambda_opt = lambda_opt
  ))
}

PCO_result <- PCO_Diagnostic_Selection(H, h_min, X, gridcal, Delta, gh_min, g.real)

###########################################################
# Comparativa y Resultados Finales
###########################################################
hcv <- bw.ucv(X)
gh.cv <- sapply(gridcal, function(z) gh(z, hcv, X))

h.sil <- bw.nrd0(X)
gh.sil <- sapply(gridcal, function(z) gh(z, h.sil, X))

cat("=== Diagnóstico PCO (Slope Heuristics) ===\n")
cat("Lambda mínimo detectado:", PCO_result$lambda_min, "\n")
cat("Lambda óptimo aplicado:", PCO_result$lambda_opt, "\n\n")

ISE.PCO <- sum((PCO_result$ghat - g.real)^2) * Delta
ISE.CV  <- sum((gh.cv - g.real)^2) * Delta
ISE.SIL <- sum((gh.sil - g.real)^2) * Delta

print(data.frame(
  Metodo = c("PCO (Corregido)", "CV", "Silverman"), 
  ISE = c(ISE.PCO, ISE.CV, ISE.SIL)
))