###########################################################
#                  Configuración de Parámetros
###########################################################
# install.packages("nor1mix") # Descomentar si no está instalado
library(nor1mix)

n <- 1000          # Tamaño de la muestra a simular
phi <- 0.9         # Coeficiente de autocorrelación (dependencia temporal)
mw_index <- 1      # Índice de la densidad de Marron & Wand a simular (1 a 15)

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

# Generación de malla dinámica basada en los cuantiles de la mezcla
grid_min <- qnorMix(0.001, mw_obj_selected)
grid_max <- qnorMix(0.999, mw_obj_selected)
gridcal <- seq(grid_min, grid_max, length=200)
Delta <- diff(gridcal)[1]

###########################################################
# Estimación de Densidad por PCO - Calibración Data-Driven
###########################################################

u <- seq(0, floor(log(n))*(2/3), by = 0.1)
H <- exp(-u)
# Ordenamos de menor a mayor para facilitar la lectura del salto
H <- sort(H) 
h_min <- min(H)

gh <- function(x, h, X) {
  n_len <- length(X)
  return(sum(dnorm((x-X)/h))/(n_len*h))
}

gh_min <- sapply(gridcal, function(z) gh(z, h_min, X))

PCO_DataDriven_Selection <- function(H, h_min, X, gridcal, Delta, gh_min) {
  n_len <- length(X)
  
  # 1. Precomputar estimadores y distancias empíricas una sola vez (Eficiencia)
  norm_diff_vec <- numeric(length(H))
  estimadores <- matrix(0, nrow=length(gridcal), ncol=length(H))
  
  for(i in seq_along(H)) {
    h <- H[i]
    gh_h <- sapply(gridcal, function(z) gh(z, h, X))
    estimadores[, i] <- gh_h
    norm_diff_vec[i] <- sum((gh_h - gh_min)^2) * Delta
  }
  
  # 2. Definir grid de lambdas para explorar la Transición de Fase
  lambdas_grid <- seq(0.01, 3, by=0.01)
  h_selected <- numeric(length(lambdas_grid))
  
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
      
      # Penalidad PCO exacta según el artículo
      pen_lam <- (lam * norm2_h - norm2_diff_K) / n_len
      criterios_lam[i] <- norm_diff_vec[i] + pen_lam
    }
    # Guardar el h óptimo seleccionado para este lambda
    h_selected[j] <- H[which.min(criterios_lam)]
  }
  
  # 4. Encontrar el salto crítico (Slope Heuristics)
  # Buscamos el lambda donde h_selected da el salto más grande (sale del overfitting)
  saltos <- diff(h_selected)
  idx_salto <- which.max(saltos)
  lambda_min_hat <- lambdas_grid[idx_salto]
  
  # 5. Aplicar la regla del doble: \lambda_{opt} = 2 * \lambda_{min}
  lambda_opt <- 2 * lambda_min_hat
  
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
  
  return(list(
    hhat = H[idx_final],
    ghat = estimadores[, idx_final],
    lambda_min = lambda_min_hat,
    lambda_opt = lambda_opt
  ))
}

# Ejecutar estimación Data-Driven
PCO_result <- PCO_DataDriven_Selection(H, h_min, X, gridcal, Delta, gh_min)

###########################################################
# Comparativa con Métodos Clásicos
###########################################################
hcv <- bw.ucv(X)
gh.cv <- sapply(gridcal, function(z) gh(z, hcv, X))

h.sil <- bw.nrd0(X)
gh.sil <- sapply(gridcal, function(z) gh(z, h.sil, X))

g.real <- dnorMix(gridcal, mw_obj_selected)

###########################################################
# Resultados, Visualización y Errores (ISE)
###########################################################
cat("=== Calibración Data-Driven (Slope Heuristics) ===\n")
cat("Lambda mínimo estimado:", PCO_result$lambda_min, "\n")
cat("Lambda óptimo aplicado:", PCO_result$lambda_opt, "\n\n")

plot(gridcal, PCO_result$ghat, type="l", lwd=2, col=4,
     ylim=c(0, max(PCO_result$ghat, gh.cv, gh.sil, g.real)), 
     ylab="Densidad", xlab="x", main=paste("PCO Data-Driven - MW Densidad", mw_index))
lines(gridcal, gh.cv, col=3, lwd=2)
lines(gridcal, gh.sil, col=6, lwd=2)
lines(gridcal, g.real, col=2, lwd=2)
legend("topright", legend=c("PCO (Data-Driven)", "CV", "Silverman", "Real"), col=c(4, 3, 6, 2), lwd=2)

ISE.PCO <- sum((PCO_result$ghat - g.real)^2) * Delta
ISE.CV  <- sum((gh.cv - g.real)^2) * Delta
ISE.SIL <- sum((gh.sil - g.real)^2) * Delta

print(data.frame(
  Metodo = c("PCO (Data-Driven)", "CV", "Silverman"), 
  ISE = c(ISE.PCO, ISE.CV, ISE.SIL)
))