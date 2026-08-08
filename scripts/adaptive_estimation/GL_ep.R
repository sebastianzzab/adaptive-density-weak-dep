###########################################################
# Calibración de Gamma - Método GL (Riesgo Local)
# Optimizado con Monte Carlo y Paralelización
# Distribución: Normal Estándar
###########################################################

library(parallel)

###########################################################
#                  Parámetros de Simulación
###########################################################
n <- 120           # Tamaño de la muestra
phi <- 0.9         # Coeficiente de autocorrelación AR(1)
B <- 100           # Iteraciones de Monte Carlo

# Cuadrícula fina para la búsqueda del parámetro gamma
GAM <- seq(0.05, 3.0, by = 0.025)

# Configuración de la malla de evaluación
gridcal <- seq(-8, 8, length=200)
Delta <- diff(gridcal)[1]

# Densidad Real (Normal) y cota superior teórica
g.real <- dnorm(gridcal)
gsup <- 1 / sqrt(2 * pi)

# Familia de ventanas H candidatas
u <- seq(0, floor(log(n))*(2/3), by = 0.1)
H <- exp(-u)

###########################################################
#   Funciones del Estimador GL (Vectorizadas)
###########################################################

gh <- function(x, h, X) {
  n_len <- length(X)
  return(sum(dnorm((x - X)/h)) / (n_len * h))
}

ghh <- function(x, h, hp, X) {
  hs <- sqrt(h^2 + hp^2)
  return(gh(x, hs, X))
}

Vh_vectorizado <- function(H_vec, n, gsup, gamma) {
  delta.n <- sqrt(log(n))
  K1 <- 1
  K2 <- sqrt(1 / (2 * sqrt(pi)))
  return(sqrt(2 * gamma * gsup) * K2 * (K1 + 1) * (1 + delta.n) * (((log(n))^(-1/2)) / sqrt(n * H_vec)))
}

Ahx <- function(x, h, H, X, V_H) {
  valores <- numeric(length(H))
  for(i in seq_along(H)) {
    hp <- H[i]
    valores[i] <- max(abs(ghh(x, h, hp, X) - gh(x, hp, X)) - V_H[i], 0)
  }
  return(max(valores))
}

hx <- function(x, H, X, V_H) {
  criterio <- numeric(length(H))
  for(i in seq_along(H)) {
    h <- H[i]
    criterio[i] <- Ahx(x, h, H, X, V_H) + V_H[i]
  }
  idx_min <- which.min(criterio)
  hhat <- H[idx_min]
  return(list(hhat = hhat, ghat = gh(x, hhat, X)))
}

GLdens <- function(gridcal, H, X, n, gsup, gamma) {
  m <- length(gridcal)
  ghat <- numeric(m)
  V_H <- Vh_vectorizado(H, n, gsup, gamma)
  for(i in 1:m) {
    x <- gridcal[i]
    res <- hx(x, H, X, V_H)
    ghat[i] <- res$ghat
  }
  return(list(ghat = ghat))
}

########################################################
#   Ejecución: Grid Search y Monte Carlo
########################################################

MISE <- numeric(length(GAM))

# Configuración del cluster paralelo (protegiendo 1 núcleo del sistema) (modificado para colab)
# num_cores <- max(1, detectCores() - 1) # comentada esta linea
num_cores <- detectCores() # Usará los 2 núcleos completos en la versión gratuita (hipotesis de Gemini)
cat(sprintf("Iniciando procesamiento en paralelo con %d núcleos...\n", num_cores))
cat(sprintf("Simulando Densidad Normal | n = %d | phi = %.2f | B = %d iteraciones\n\n", n, phi, B))

Tiempo <- system.time({
  for(I in seq_along(GAM)) {
    gamma_actual <- GAM[I]
    
    run_monte_carlo <- function(b) {
      # 1. Generación AR(1) dependiente
      if (phi == 0) {
        Z <- rnorm(n)
      } else {
        Z <- arima.sim(model = list(ar = phi), n = n)
      }
      sd.Z <- if (phi == 0) 1 else sqrt(1 / (1 - phi^2))
      U <- pnorm(Z, mean = 0, sd = sd.Z)
      
      # 2. Mapeo a densidad Normal (Transformación cuantil)
      X <- qnorm(U)
      
      # 3. Estimación GL
      GLgamma <- GLdens(gridcal, H, X, n, gsup, gamma_actual)
      
      # 4. Cálculo de Error (ISE)
      sum((GLgamma$ghat - g.real)^2) * Delta
    }
    
    # Manejo multiplataforma para ejecución concurrente
    if (.Platform$OS.type == "unix") {
      ISE_resultados <- mclapply(1:B, run_monte_carlo, mc.cores = num_cores)
    } else {
      cl <- makeCluster(num_cores)
      clusterExport(cl, c("phi", "n", "gridcal", "H", "gsup", "gamma_actual",
                          "g.real", "Delta", "GLdens", "hx", "Ahx",
                          "Vh_vectorizado", "ghh", "gh"))
      ISE_resultados <- parLapply(cl, 1:B, run_monte_carlo)
      stopCluster(cl)
    }
    
    # Extraer y promediar los resultados para obtener el MISE
    ISE <- unlist(ISE_resultados)
    MISE[I] <- mean(ISE)
    
    cat(sprintf("Evaluando gamma = %.2f | MISE = %.6f\n", gamma_actual, MISE[I]))
  }
})

cat("\n========================================\n")
cat("Tiempo total de ejecución:\n")
print(Tiempo)

########################################################
#   Presentación Gráfica de Resultados (Mejorada)
########################################################

# Identificar el óptimo global
idx_optimo <- which.min(MISE)
gamma_optimo <- GAM[idx_optimo]
mise_minimo <- MISE[idx_optimo]

cat(sprintf("\n=> RESULTADO: El gamma óptimo estimado es %.2f (MISE = %.6f)\n", gamma_optimo, mise_minimo))
cat("========================================\n")

# # Gráfico de la curva MISE con mejoras estéticas
# plot(GAM, MISE, type = "b", pch = 19, col = "darkblue", lwd = 2,
#      xlab = expression(gamma), ylab = "MISE (Error Cuadrático Medio Integrado)",
#      main = expression(paste("Optimización de ", gamma, " en Estimador GL (Proceso AR(1) con phi=0.9 y n=1000)")),
#      cex.main = 1.2, cex.lab = 1.1)
# grid(col = "lightgray", lty = "dotted", lwd = 1)

# Gráfico de la curva MISE con mejoras estéticas
plot(GAM, MISE, type = "b", pch = 19, col = "darkblue", lwd = 2,
     xlab = expression(gamma), ylab = "MISE (Error Cuadrático Medio Integrado)",
     cex.main = 1.2, cex.lab = 1.1)
grid(col = "lightgray", lty = "dotted", lwd = 1)

# Añadir líneas de referencia y resaltar el punto óptimo
abline(v = gamma_optimo, col = "red", lty = 2, lwd = 2)
points(gamma_optimo, mise_minimo, col = "red", cex = 2, lwd = 2, pch = 19)

# Etiqueta de texto posicional para el valor óptimo
# Se coloca ligeramente por encima del mínimo para evitar superposición con la curva
# Etiqueta de texto posicional para el valor óptimo
desplazamiento_y <- (max(MISE) - min(MISE)) * 0.05
text(gamma_optimo, mise_minimo + desplazamiento_y,
     labels = bquote(gamma ~ "óptimo =" ~ .(sprintf("%.2f", gamma_optimo))),
     pos = 4, col = "red", font = 2, cex = 1.1)