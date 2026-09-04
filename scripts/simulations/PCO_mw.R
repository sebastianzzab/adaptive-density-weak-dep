###########################################################
#                  Configuración de Parámetros
###########################################################
# install.packages("nor1mix") # Descomentar si no está instalado
library(nor1mix)

n <- 200          # Tamaño de la muestra a simular
phi <- 0        # Coeficiente de autocorrelación (dependencia temporal)
# phi <- 0.9         # Coeficiente de autocorrelación (dependencia temporal)
mw_index <- 9      # Índice de la densidad de Marron & Wand a simular (1 a 15)

###########################################################
# Configuración de las Densidades de Marron & Wand
###########################################################
# El paquete nor1mix provee objetos predefinidos MW.nm1 hasta MW.nm15.
mw_list <- list(
  MW.nm1,  MW.nm2,  MW.nm3,  MW.nm4,  MW.nm5,
  MW.nm6,  MW.nm7,  MW.nm8,  MW.nm9,  MW.nm10,
  MW.nm11, MW.nm12, MW.nm13, MW.nm14, MW.nm15
)

###########################################################
# FUNCION PARA GENERAR VARIABLES UNIFORMES DEPENDIENTES
###########################################################
generar_uniformes_ar1 <- function(nsim, phi) {
  # Generación del proceso AR(1)
  if (phi == 0) {
    Z <- rnorm(nsim)
  } else {
    Z <- arima.sim(list(ar=phi), n=nsim) 
  }
  
  # Cálculo de la varianza teórica
  sd_Z <- if (phi == 0) 1 else sqrt(1 / (1 - phi^2))
  
  # Retornar variables uniformes U
  return(pnorm(Z, mean=0, sd=sd_Z))
}

###########################################################
# FUNCIÓN GENERADORA VARIABLES CON NORMAL TRUNCADA DEPENDIENTE AR(1)
###########################################################
gnormtrunc_dep <- function(nsim, phi, a, b, m, de) {
  U <- generar_uniformes_ar1(nsim, phi)
  
  F_a <- pnorm(a, mean = m, sd = de)
  F_b <- pnorm(b, mean = m, sd = de)
  p <- F_b - F_a
  
  return(qnorm(p * U + F_a, mean=m, sd=de))
}

###########################################################
# FUNCIÓN GENERADORA MEDIANTE TRANSFORMADA INVERSA GENERALIZADA 
###########################################################
generar_dependiente_mw <- function(nsim, phi, mw_obj) {
  U <- generar_uniformes_ar1(nsim, phi)
  return(qnorMix(U, mw_obj))
}

###########################################################
# Generación de la Muestra Objetivo
###########################################################
mw_obj_selected <- mw_list[[mw_index]]

# Ejecución de la generación
X <- generar_dependiente_mw(n, phi, mw_obj_selected)

# Generación de malla dinámica basada en los cuantiles de la mezcla seleccionada
grid_min <- qnorMix(0.001, mw_obj_selected)
grid_max <- qnorMix(0.999, mw_obj_selected)
gridcal <- seq(grid_min, grid_max, length=200)
Delta <- diff(gridcal)[1]

###########################################################
# Estimación de Densidad por Método PCO 
###########################################################

# Definición de la familia de ventanas candidatas (H)
u <- seq(0, floor(log(n))*(2/3), by = 0.1)
H <- exp(-u)
h_min <- min(H)

# Estimador de Núcleo (Kernel) Gaussiano básico
gh <- function(x, h, X) {
  n_len <- length(X)
  return(sum(dnorm((x-X)/h))/(n_len*h))
}

# Producto interno exacto de dos núcleos Gaussianos para el cálculo de penalidad
inner_prod_K <- function(h1, h2) {
  return(1 / sqrt(2 * pi * (h1^2 + h2^2)))
}

# Evaluación del estimador de "overfitting" en la malla
gh_min <- sapply(gridcal, function(z) gh(z, h_min, X))

# Función principal para seleccionar el ancho de banda usando PCO
PCO_Selection <- function(H, h_min, X, gridcal, Delta, gh_min) {
  n_len <- length(X)
  criterios <- numeric(length(H))
  estimadores <- matrix(0, nrow=length(gridcal), ncol=length(H))
  
  for(i in seq_along(H)) {
    h <- H[i]
    
    # 1. Calcular estimador para ventana h
    gh_h <- sapply(gridcal, function(z) gh(z, h, X))
    estimadores[, i] <- gh_h
    
    # 2. Calcular la norma L2 empírica (ISE discreto) respecto al overfitting
    norm_diff <- sum((gh_h - gh_min)^2) * Delta
    
    # 3. Calcular penalidad óptima teórica
    pen_opt <- 2 * inner_prod_K(h, h_min) / n_len
    
    # 4. Criterio PCO
    criterios[i] <- norm_diff + pen_opt
  }
  
  # Seleccionar el estimador que minimiza el criterio
  idx_optimo <- which.min(criterios)
  
  return(list(
    hhat = H[idx_optimo],
    ghat = estimadores[, idx_optimo],
    criterios = criterios
  ))
}

# Ejecutar estimación PCO
PCO_result <- PCO_Selection(H, h_min, X, gridcal, Delta, gh_min)

###########################################################
# Comparativa con Métodos Clásicos
###########################################################

# 1. Cross-Validation (Validación Cruzada)
hcv <- bw.ucv(X)
gh.cv <- sapply(gridcal, function(z) gh(z, hcv, X))

# 2. Regla de Silverman
h.sil <- bw.nrd0(X)
gh.sil <- sapply(gridcal, function(z) gh(z, h.sil, X))

# Determinar la densidad real exacta para el cálculo del error (ISE)
g.real <- dnorMix(gridcal, mw_obj_selected)

###########################################################
# Resultados, Visualización y Errores (ISE)
###########################################################

# Gráfico comparativo de todas las estimaciones
plot(gridcal, PCO_result$ghat, type="l", lwd=2, col=4,
     ylim=c(0, max(PCO_result$ghat, gh.cv, gh.sil, g.real)), 
     ylab="Densidad", xlab="x", main=paste("Estimación MW Densidad", mw_index))
lines(gridcal, gh.cv, col=3, lwd=2)
lines(gridcal, gh.sil, col=6, lwd=2)
lines(gridcal, g.real, col=2, lwd=2)
legend("topright", legend=c("PCO", "CV", "Silverman", "Real"), col=c(4, 3, 6, 2), lwd=2)

# Cálculo del Error Cuadrático Integrado Aproximado (ISE)
ISE.PCO <- sum((PCO_result$ghat - g.real)^2) * Delta
ISE.CV  <- sum((gh.cv - g.real)^2) * Delta
ISE.SIL <- sum((gh.sil - g.real)^2) * Delta

print(data.frame(
  Metodo = c("PCO", "CV", "Silverman"), 
  ISE = c(ISE.PCO, ISE.CV, ISE.SIL)
))