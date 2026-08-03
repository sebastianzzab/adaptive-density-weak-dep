###########################################################
# 1. Configuración Global y Paquetes
###########################################################
if (!require(nor1mix)) install.packages("nor1mix")
library(nor1mix)

set.seed(42)       # Para reproducibilidad
M <- 5             # Iteraciones de Monte Carlo
n <- 250           # Tamaño de muestra nominal
phi <- 0.5         # Dependencia AR(1)
gamma_fijo <- 0.19 # Parámetro Gamma original 

# Cálculo del Tamaño de Muestra Efectivo (AR(1))
n_eff <- n * ((1 - phi) / (1 + phi))

target_density <- MW.nm2

# Familia de ventanas candidatas (H)
u <- seq(0, floor(log(n))*(2/3), by = 0.1)
H <- sort(exp(-u))

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
# 3. Funciones del Estimador GL (Con Inyección de n_penalizacion)
###########################################################
# El estimador empírico SIEMPRE usa la longitud real de los datos (n)
gh <- function(x, h, X) {
  n_len <- length(X)
  return(sum(dnorm((x-X)/h))/(n_len*h))
}

ghh <- function(x, h, hp, X) {
  hs <- sqrt(h^2 + hp^2)
  return(gh(x, hs, X))
}

# La función de penalidad usa el tamaño de muestra (nominal o efectivo)
Vh <- function(h, n_penalizacion, gsup, gamma) {
  delta.n <- sqrt(log(n_penalizacion))
  K1 <- 1
  K2 <- sqrt(1/(2*sqrt(pi)))
  return(sqrt(2*gamma*gsup)*K2*(K1+1)*(1+delta.n)*(((log(n_penalizacion))^(-1/2))/sqrt(n_penalizacion*h)))
}

Ahx <- function(x, h, H, X, n_penalizacion, gsup, gamma) {
  valores <- numeric(length(H))
  for(i in seq_along(H)) {
    hp <- H[i]
    valores[i] <- max(abs(ghh(x, h, hp, X) - gh(x, hp, X)) - Vh(hp, n_penalizacion, gsup, gamma), 0)
  }
  return(max(valores))
}

hx <- function(x, H, X, n_penalizacion, gsup, gamma) {
  criterio <- numeric(length(H))
  for(i in seq_along(H)) {
    h <- H[i]
    criterio[i] <- Ahx(x, h, H, X, n_penalizacion, gsup, gamma) + Vh(h, n_penalizacion, gsup, gamma)
  }
  hhat <- H[which.min(criterio)]
  return(list(hhat=hhat, ghat=gh(x, hhat, X)))
}

GLdens <- function(gridcal, H, X, n_penalizacion, gsup, gamma) {
  m <- length(gridcal)
  ghat <- numeric(m)
  for(i in 1:m) {
    res <- hx(gridcal[i], H, X, n_penalizacion, gsup, gamma)
    ghat[i] <- res$ghat
  }
  return(ghat)
}

###########################################################
# 4. Ejecución de la Simulación Monte Carlo
###########################################################
resultados <- data.frame(
  Iteracion = 1:M,
  ISE_GL_Nominal = numeric(M),
  ISE_GL_Ajustado  = numeric(M)
)

grid_min <- qnorMix(0.001, target_density)
grid_max <- qnorMix(0.999, target_density)
gridcal <- seq(grid_min, grid_max, length=80) 
Delta <- diff(gridcal)[1]
g.real <- dnorMix(gridcal, target_density)
gsup <- max(g.real)

cat(sprintf("Iniciando Monte Carlo (M=%d) | phi=%.2f\n", M, phi))
cat(sprintf("n Nominal: %d | n Efectivo (n_eff): %.2f\n\n", n, n_eff))

par(mfrow=c(2,3), mar=c(4,4,3,1))

for (m in 1:M) {
  cat(sprintf("Iteración %d... ", m))
  
  X_m <- generar_dependiente_mw(n, phi, target_density)
  
  # Estimación GL usando n nominal (Clásico)
  ghat_nominal <- GLdens(gridcal, H, X_m, n, gsup, gamma_fijo)
  
  # Estimación GL usando n efectivo (Ajustado por Dependencia)
  ghat_ajustado <- GLdens(gridcal, H, X_m, n_eff, gsup, gamma_fijo)
  
  # Cálculo de Errores (ISE)
  ise_nominal <- sum((ghat_nominal - g.real)^2) * Delta
  ise_ajustado <- sum((ghat_ajustado - g.real)^2) * Delta
  
  resultados$ISE_GL_Nominal[m] <- ise_nominal
  resultados$ISE_GL_Ajustado[m] <- ise_ajustado
  
  # Gráficos Individuales
  max_y <- max(c(ghat_nominal, ghat_ajustado, g.real))
  plot(gridcal, g.real, type="l", col="black", lwd=2, ylim=c(0, max_y),
       main=paste("Iteración", m), xlab="x", ylab="Densidad")
  lines(gridcal, ghat_nominal, col="darkgray", lwd=2, lty=2)
  lines(gridcal, ghat_ajustado, col="blue", lwd=2, lty=1)
  
  if(m == 1) {
    legend("topright", legend=c("Real", "GL Nominal", "GL Ajustado"), 
           col=c("black", "darkgray", "blue"), 
           lwd=2, lty=c(1,2,1), bty="n", cex=0.8)
  }
  cat("Completada.\n")
}

# Cálculo del MISE
mise_nominal <- mean(resultados$ISE_GL_Nominal)
mise_ajustado <- mean(resultados$ISE_GL_Ajustado)

###########################################################
# 5. Gráfico Resumen y Salida de Datos
###########################################################
plot(1:M, resultados$ISE_GL_Nominal, type="b", col="darkgray", lwd=2, pch=19,
     ylim=c(0, max(c(resultados$ISE_GL_Nominal, resultados$ISE_GL_Ajustado))),
     main="Evolución del ISE por Iteración", xlab="Iteración", ylab="ISE")
lines(1:M, resultados$ISE_GL_Ajustado, type="b", col="blue", lwd=2, pch=19)
abline(h=mise_nominal, col="darkgray", lty=3, lwd=2)
abline(h=mise_ajustado, col="blue", lty=3, lwd=2)
legend("topleft", legend=c("ISE GL Nominal", "ISE GL Ajustado", "MISE Nominal", "MISE Ajustado"),
       col=c("darkgray", "blue", "darkgray", "blue"),
       lwd=2, lty=c(1,1,3,3), pch=c(19,19,NA,NA), bty="n", cex=0.8)

par(mfrow=c(1,1))

cat("\n=== Resultados de la Simulación Monte Carlo ===\n")
print(resultados)
cat("\n-----------------------------------------------\n")
cat(sprintf("MISE (GL Nominal): %.6f\n", mise_nominal))
cat(sprintf("MISE (GL Ajustado): %.6f\n", mise_ajustado))
cat("-----------------------------------------------\n")