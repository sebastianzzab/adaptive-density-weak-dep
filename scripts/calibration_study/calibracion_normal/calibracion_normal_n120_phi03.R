###########################################################
# Calibración de Gamma para Normal - Método GL
###########################################################

# Instalar y cargar matrixStats si no está disponible (esencial para velocidad)
if (!require(matrixStats)) install.packages("matrixStats", repos = "http://cran.us.r-project.org")
library(parallel)
library(matrixStats)

# Establecer directorio actual en la ubicación del archivo de trabajo
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

directorio_salida <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###########################################################
#                  Parámetros de Simulación
###########################################################
n <- 120           # Tamaño de la muestra
phi <- 0.3         # Coeficiente de autocorrelación AR(1)
B <- 300           # Iteraciones de Monte Carlo

# Selector dinámico de densidad
densidad <- "normal" 

# Cuadrícula fina para la búsqueda del parámetro gamma
GAM <-c(seq(0.0025,0.0275,by=0.0025),seq(0.030,0.100,by=0.0025),seq(0.125,0.5,by=0.00625))
# GAM <- c(
#   seq(0.20, 0.40, by = 0.01),   # Subida rápida hacia la zona de interés
#   seq(0.405, 0.550, by = 0.005), # Alta resolución alrededor de tu óptimo actual (0.4375)
#   seq(0.56, 1.00, by = 0.02)    # Extensión de seguridad para asegurar que la curva suba
# )
# GAM <- seq(0.350, 0.500, by = 0.005)
# GAM <- seq(0.400, 0.500, by = 0.0005)

###########################################################
# Funciones auxiliares para la mezcla
###########################################################
GMix <- function(x) {
  return(0.5*pnorm(x, mean=-2, sd=1) + 0.5*pnorm(x, mean=2, sd=1))
}

grid.mix.x <- seq(-8, 8, length.out=10000)
grid.mix.u <- GMix(grid.mix.x)

QMix <- function(u) {
  approx(x = grid.mix.u, y = grid.mix.x, xout = u, rule = 2)$y
}

###########################################################
# Inicialización Dinámica del Entorno Global
###########################################################
if(densidad == "normal") {
  gsup <- 1/sqrt(2*pi)
  gridcal <- seq(-8, 8, length=200) 
  g.real <- dnorm(gridcal)
} else if(densidad == "lognormal") {
  gsup <- exp(0.125)/(0.5*sqrt(2*pi))
  gridcal <- seq(0, qlnorm(0.999999, 0, 0.5), length=200) 
  g.real <- dlnorm(gridcal, meanlog=0, sdlog=0.5)
} else if(densidad == "mezcla") {
  gsup <- 0.5*dnorm(-2, -2, 1) + 0.5*dnorm(-2, 2, 1)
  gridcal <- seq(-8, 8, length=200) 
  g.real <- 0.5*dnorm(gridcal, -2, 1) + 0.5*dnorm(gridcal, 2, 1)
} else {
  stop("Densidad no implementada")
}

Delta <- diff(gridcal)[1]

# Familia de ventanas H candidatas
u <- seq(0, floor(log(n))*(2/3), by = 0.1)
H <- exp(-u)

###########################################################
#   Función de Penalización Vectorizada
###########################################################
Vh_vectorizado <- function(H_vec, n, gsup, gamma, phi) {
  if  (phi == 0) {
    delta.n <- 0
  } else {
    delta.n <- sqrt(log(n))
  }
  K1 <- 1
  K2 <- sqrt(1 / (2 * sqrt(pi)))
  return(sqrt(2 * gamma * gsup) * K2 * (K1 + 1) * (1 + delta.n) * (((log(n))^(-1/2)) / sqrt(n * H_vec)))
}

########################################################
#   Ejecución: Grid Search y Monte Carlo
########################################################
MISE <- numeric(length(GAM))

# num_cores <- detectCores() 
num_cores <- max(1, detectCores() - 2)
cat(sprintf("Iniciando procesamiento paralelo MATRICIAL con %d núcleos...\n", num_cores))
cat(sprintf("Simulando Densidad: %s | n = %d | phi = %.2f | B = %d iteraciones\n\n", toupper(densidad), n, phi, B))

Tiempo <- system.time({
  
  for(I in seq_along(GAM)) {
    gamma_actual <- GAM[I]
    
    # Precalculamos V(h) una sola vez para este gamma (fuera del Monte Carlo)
    V_H <- Vh_vectorizado(H, n, gsup, gamma_actual, phi)
    
    run_monte_carlo <- function(b) {
      # 1. Generación AR(1) dependiente
      if (phi == 0) {
        Z <- rnorm(n)
      } else {
        Z <- arima.sim(model = list(ar = phi), n = n)
      }
      sd.Z <- if (phi == 0) 1 else sqrt(1 / (1 - phi^2))
      U <- pnorm(Z, mean = 0, sd = sd.Z)
      
      # 2. Mapeo a la densidad seleccionada
      if(densidad == "normal") { X <- qnorm(U)
      } else if(densidad == "lognormal") { X <- qlnorm(U, meanlog=0, sdlog=0.5)
      } else if(densidad == "mezcla") { X <- QMix(U) }
      
      # ============================================================
      # 3. MOTOR MATRICIAL DEL ESTIMADOR GL
      # ============================================================
      
      # Matriz base de distancias (100 puntos de malla x tamaño de muestra n)
      # Esto elimina miles de restas repetitivas
      dist_base <- outer(gridcal, X, "-")
      
      # Precomputar matriz de estimaciones para cada ventana h (Matriz GH)
      Matriz_GH <- matrix(0, nrow = length(gridcal), ncol = length(H))
      for(j in seq_along(H)) {
        Matriz_GH[, j] <- rowSums(dnorm(dist_base / H[j])) / (n * H[j])
      }
      
      # Matriz para almacenar el criterio A(h,x) de toda la malla simultáneamente
      Criterio_A_Matriz <- matrix(0, nrow = length(gridcal), ncol = length(H))
      
      for(i in seq_along(H)) {
        h_actual <- H[i]
        
        # Construir matriz de comparaciones sobresuavizadas (GHH) en memoria dinámica
        GHH_temp <- matrix(0, nrow = length(gridcal), ncol = length(H))
        for(j in seq_along(H)) {
          hp_actual <- H[j]
          hs <- sqrt(h_actual^2 + hp_actual^2)
          GHH_temp[, j] <- rowSums(dnorm(dist_base / hs)) / (n * hs)
        }
        
        # Calcular el penalizado restando matrices completas (Álgebra Lineal)
        Diferencia <- abs(GHH_temp - Matriz_GH)
        Penalizado <- sweep(Diferencia, MARGIN = 2, STATS = V_H, FUN = "-")
        Penalizado[Penalizado < 0] <- 0 
        
        # Encontrar el máximo de Lepski por fila (gracias a matrixStats en C)
        Criterio_A_Matriz[, i] <- rowMaxs(Penalizado)
      }
      
      # Selección automática de la ventana final 
      Criterio_Final <- sweep(Criterio_A_Matriz, MARGIN = 2, STATS = V_H, FUN = "+")
      indices_minimos <- max.col(-Criterio_Final, ties.method = "first")
      
      # Extraer las estimaciones finales óptimas
      ghat_final <- Matriz_GH[cbind(1:length(gridcal), indices_minimos)]
      
      # ============================================================
      # 4. Cálculo de Error (ISE)
      # ============================================================
      return(sum((ghat_final - g.real)^2) * Delta)
    }
    
    # Manejo multiplataforma
    if (.Platform$OS.type == "unix") {
      ISE_resultados <- mclapply(1:B, run_monte_carlo, mc.cores = num_cores)
    } else {
      cl <- makeCluster(num_cores)
      clusterExport(cl, c("phi", "n", "gridcal", "H", "gsup", "V_H",
                          "g.real", "Delta", "densidad", "rowMaxs",
                          "QMix", "grid.mix.u", "grid.mix.x"))
      ISE_resultados <- parLapply(cl, 1:B, run_monte_carlo)
      stopCluster(cl)
    }
    
    ISE <- unlist(ISE_resultados)
    MISE[I] <- mean(ISE)
    
    cat(sprintf("Evaluando gamma = %.3f | MISE = %.6f\n", gamma_actual, MISE[I]))
  }
})

cat("\n========================================\n")
cat("Tiempo total de ejecución:\n")
print(Tiempo)

########################################################
#   Presentación Gráfica de Resultados
########################################################
idx_optimo <- which.min(MISE)
gamma_optimo <- GAM[idx_optimo]
mise_minimo <- MISE[idx_optimo]

cat(sprintf("\n=> RESULTADO: El gamma óptimo estimado es %.3f (MISE = %.6f)\n", gamma_optimo, mise_minimo))
cat("========================================\n")

###############################################################################
if(n == 60 & phi==0) {
  nombre_imagen <- "Resultados_Normal_n60_phi0.png"
} else if(n == 60 & phi==0.5) {
  nombre_imagen <- "Resultados_Normal_n60_phi05.png"
} else if(n == 60 & phi==0.7) {
  nombre_imagen <- "Resultados_Normal_n60_phi07.png"
} else if(n == 60 & phi==0.8) {
  nombre_imagen <- "Resultados_Normal_n60_phi08.png"
} else if(n == 60 & phi==0.85) {
  nombre_imagen <- "Resultados_Normal_n60_phi085.png"
} else if(n == 60 & phi==0.9) {
  nombre_imagen <- "Resultados_Normal_n60_phi09.png"
} else if(n == 120 & phi==0) {
  nombre_imagen <- "Resultados_Normal_n120_phi0.png"
} else if(n == 120 & phi==0.5) {
  nombre_imagen <- "Resultados_Normal_n120_phi05.png"
} else if(n == 120 & phi==0.7) {
  nombre_imagen <- "Resultados_Normal_n120_phi07.png"
} else if(n == 120 & phi==0.8) {
  nombre_imagen <- "Resultados_Normal_n120_phi08.png"
} else if(n == 120 & phi==0.85) {
  nombre_imagen <- "Resultados_Normal_n120_phi085.png"
} else if(n == 120 & phi==0.9) {
  nombre_imagen <- "Resultados_Normal_n120_phi09.png"
} else if(n == 250 & phi==0) {
  nombre_imagen <- "Resultados_Normal_n250_phi0.png"
} else if(n == 250 & phi==0.5) {
  nombre_imagen <- "Resultados_Normal_n250_phi05.png"
} else if(n == 250 & phi==0.7) {
  nombre_imagen <- "Resultados_Normal_n250_phi07.png"
} else if(n == 250 & phi==0.8) {
  nombre_imagen <- "Resultados_Normal_n250_phi08.png"
} else if(n == 250 & phi==0.85) {
  nombre_imagen <- "Resultados_Normal_n250_phi085.png"
}else if(n == 250 & phi==0.9) {
  nombre_imagen <- "Resultados_Normal_n250_phi09.png"
}else {
  nombre_imagen <- "Resultados_Normal.png"
}


png(
  filename = file.path(directorio_salida, nombre_imagen),
  width = 1500,
  height = 500,
  res = 120
)

plot(GAM, MISE, type = "b", pch = 19, col = "darkblue", lwd = 2,
     xlab = expression(gamma), ylab = "MISE (Error Cuadrático Medio Integrado)",
     # main = bquote("Optimización Matricial:" ~ gamma ~ "en GL (AR(1) phi="*.(phi)*" y n="*.(n)*")"),
     cex.main = 1.2, cex.lab = 1.1)
grid(col = "lightgray", lty = "dotted", lwd = 1)

abline(v = gamma_optimo, col = "red", lty = 2, lwd = 2)
points(gamma_optimo, mise_minimo, col = "red", cex = 2, lwd = 2, pch = 19)

desplazamiento_y <- (max(MISE) - min(MISE)) * 0.05
text(gamma_optimo, mise_minimo + desplazamiento_y,
     labels = bquote(gamma ~ "óptimo =" ~ .(sprintf("%.3f", gamma_optimo))),
     pos = 4, col = "red", font = 2, cex = 1.1)
grid()
dev.off()

######################################################################
Resultados <- list(GAM = GAM,MISE = MISE, GAM_OPT = gamma_optimo, 
                   MISE_MIN = mise_minimo, DES = desplazamiento_y)
###########################################################
# Resultados
###########################################################
if(n == 60 & phi==0) {
  save(Resultados,file = "Resultados_Normal_n60_phi0.RData")
} else if(n == 60 & phi==0.3) {
  save(Resultados,file = "Resultados_Normal_n60_phi03.RData")
} else if(n == 60 & phi==0.5) {
  save(Resultados,file = "Resultados_Normal_n60_phi05.RData")
} else if(n == 60 & phi==0.7) {
  save(Resultados,file = "Resultados_Normal_n60_phi07.RData")
} else if(n == 60 & phi==0.85) {
  save(Resultados,file = "Resultados_Normal_n60_phi085.RData")
} else if(n == 60 & phi==0.9) {
  save(Resultados,file = "Resultados_Normal_n60_phi09.RData")
} else if(n == 120 & phi==0) {
  save(Resultados,file = "Resultados_Normal_n120_phi0.RData")
} else if(n == 120 & phi==0.3) {
  save(Resultados,file = "Resultados_Normal_n120_phi03.RData")
} else if(n == 120 & phi==0.5) {
  save(Resultados,file = "Resultados_Normal_n120_phi05.RData")
} else if(n == 120 & phi==0.7) {
  save(Resultados,file = "Resultados_Normal_n120_phi07.RData")
} else if(n == 120 & phi==0.85) {
  save(Resultados,file = "Resultados_Normal_n120_phi085.RData")
} else if(n == 120 & phi==0.9) {
  save(Resultados,file = "Resultados_Normal_n120_phi09.RData")
} else if(n == 250 & phi==0) {
  save(Resultados,file = "Resultados_Normal_n250_phi0.RData")
} else if(n == 250 & phi==0.3) {
  save(Resultados,file = "Resultados_Normal_n250_phi03.RData")
} else if(n == 250 & phi==0.5) {
  save(Resultados,file = "Resultados_Normal_n250_phi05.RData")
} else if(n == 250 & phi==0.7) {
  save(Resultados,file = "Resultados_Normal_n250_phi07.RData")
} else if(n == 250 & phi==0.8) {
  save(Resultados,file = "Resultados_Normal_n250_phi08.RData")
} else if(n == 250 & phi==0.85) {
  save(Resultados,file = "Resultados_Normal_n250_phi085.RData")
} else if(n == 250 & phi==0.9) {
  save(Resultados,file = "Resultados_Normal_n250_phi09.RData")
}else {
  stop("Entradas Invalidas")
}