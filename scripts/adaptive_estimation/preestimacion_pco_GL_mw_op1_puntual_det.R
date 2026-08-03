###########################################################
# 1. Configuración Global y Paquetes
###########################################################
if (!require(nor1mix)) install.packages("nor1mix")
library(nor1mix)

set.seed(42)       # Para reproducibilidad
M <- 5             # Iteraciones de Monte Carlo
n <- 250           # Tamaño de muestra por iteración (reducido por costo computacional del GL local)
phi <- 0.5         # Dependencia AR(1)
gamma_fijo <- 0.19 # Parámetro Gamma fijo (baseline del código original)

# Selección de la densidad objetivo de Marron-Wand (ej. MW 2: Bimodal)
target_density <- MW.nm2

# Familia de ventanas candidatas (H)
u <- seq(0, floor(log(n))*(2/3), by = 0.1)
H <- sort(exp(-u))
h_min <- min(H)

###########################################################
# 2. Funciones de Generación de Datos
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
# 3. Funciones del Estimador GL Clásico (Riesgo Local)
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
# 4. Funciones de Calibración PCO a GL
###########################################################
Obtener_Lambda_PCO <- function(H, h_min, X, gridcal, Delta) {
  n_len <- length(X)
  gh_min <- sapply(gridcal, function(z) gh(z, h_min, X))
  
  norm_diff_vec <- numeric(length(H))
  for(i in seq_along(H)) {
    gh_h <- sapply(gridcal, function(z) gh(z, H[i], X))
    norm_diff_vec[i] <- sum((gh_h - gh_min)^2) * Delta
  }
  
  lambdas_grid <- seq(-2, 3, by=0.05)
  h_selected <- numeric(length(lambdas_grid))
  norm2_hmin <- 1 / (2 * sqrt(pi) * h_min)
  
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
    h_selected[j] <- H[which.min(criterios_lam)]
  }
  
  saltos <- diff(h_selected)
  lambda_min_hat <- lambdas_grid[which.max(saltos)]
  return(lambda_min_hat + 1) # Regla óptima teórica PCO
}

Transformar_Lambda_a_Gamma <- function(lambda, n, gsup) {
  numerador <- lambda * log(n)
  denominador <- 4 * gsup * (1 + sqrt(log(n)))^2
  return(numerador / denominador)
}

###########################################################
# 5. Ejecución de la Simulación Monte Carlo
###########################################################
# Inicialización de matriz para guardar resultados
resultados <- data.frame(
  Iteracion = 1:M,
  Gamma_PCO = numeric(M),
  ISE_GL_Fijo = numeric(M),
  ISE_GL_PCO  = numeric(M)
)

# Definición de la malla de evaluación
grid_min <- qnorMix(0.001, target_density)
grid_max <- qnorMix(0.999, target_density)
gridcal <- seq(grid_min, grid_max, length=80) # Resolución equilibrada para el costo del GL
Delta <- diff(gridcal)[1]
g.real <- dnorMix(gridcal, target_density)
gsup <- max(g.real)

cat(sprintf("Iniciando Monte Carlo (M=%d) - Densidad MW | n=%d | phi=%.2f\n", M, n, phi))

# Configuración del área de gráficos (2 filas, 3 columnas)
par(mfrow=c(2,3), mar=c(4,4,3,1))

for (m in 1:M) {
  cat(sprintf("Iteración %d... ", m))
  
  # 1. Generación de muestra
  X_m <- generar_dependiente_mw(n, phi, target_density)
  
  # 2. Calibración dinámica del Gamma vía PCO
  lambda_opt <- Obtener_Lambda_PCO(H, h_min, X_m, gridcal, Delta)
  gamma_pco <- Transformar_Lambda_a_Gamma(lambda_opt, n, gsup)
  
  # 3. Estimación GL con Gamma Fijo vs Gamma PCO
  ghat_fijo <- GLdens(gridcal, H, X_m, n, gsup, gamma_fijo)
  ghat_pco  <- GLdens(gridcal, H, X_m, n, gsup, gamma_pco)
  
  # 4. Cálculo del ISE
  ise_fijo <- sum((ghat_fijo - g.real)^2) * Delta
  ise_pco  <- sum((ghat_pco - g.real)^2) * Delta
  
  # 5. Guardado de datos
  resultados$Gamma_PCO[m] <- gamma_pco
  resultados$ISE_GL_Fijo[m] <- ise_fijo
  resultados$ISE_GL_PCO[m]  <- ise_pco
  
  # 6. Gráfico de la iteración actual
  max_y <- max(c(ghat_fijo, ghat_pco, g.real))
  plot(gridcal, g.real, type="l", col="black", lwd=2, ylim=c(0, max_y),
       main=paste("Iteración", m), xlab="x", ylab="Densidad")
  lines(gridcal, ghat_fijo, col="darkgray", lwd=2, lty=2)
  lines(gridcal, ghat_pco, col="darkgreen", lwd=2, lty=1)
  
  if(m == 1) {
    legend("topright", legend=c("Real", "GL Fijo", "GL PCO"), 
           col=c("black", "darkgray", "darkgreen"), 
           lwd=2, lty=c(1,2,1), bty="n", cex=0.8)
  }
  cat("Completada.\n")
}

# 6. Cálculo del MISE
mise_fijo <- mean(resultados$ISE_GL_Fijo)
mise_pco  <- mean(resultados$ISE_GL_PCO)

###########################################################
# 6. Gráfico Resumen y Salida de Datos
###########################################################
# Gráfico comparativo de ISE
plot(1:M, resultados$ISE_GL_Fijo, type="b", col="darkgray", lwd=2, pch=19,
     ylim=c(0, max(c(resultados$ISE_GL_Fijo, resultados$ISE_GL_PCO))),
     main="Evolución del ISE por Iteración", xlab="Iteración", ylab="ISE")
lines(1:M, resultados$ISE_GL_PCO, type="b", col="darkgreen", lwd=2, pch=19)
abline(h=mise_fijo, col="darkgray", lty=3, lwd=2)
abline(h=mise_pco, col="darkgreen", lty=3, lwd=2)
legend("topleft", legend=c("ISE GL Fijo", "ISE GL PCO", "MISE Fijo", "MISE PCO"),
       col=c("darkgray", "darkgreen", "darkgray", "darkgreen"),
       lwd=2, lty=c(1,1,3,3), pch=c(19,19,NA,NA), bty="n", cex=0.8)

# Resetear parámetros gráficos
par(mfrow=c(1,1))

cat("\n=== Resultados de la Simulación Monte Carlo ===\n")
print(resultados)
cat("\n-----------------------------------------------\n")
cat(sprintf("MISE (GL Fijo %g): %.6f\n", gamma_fijo, mise_fijo))
cat(sprintf("MISE (GL Calibrado PCO): %.6f\n", mise_pco))
cat("-----------------------------------------------\n")