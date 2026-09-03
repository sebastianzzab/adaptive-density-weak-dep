rm(list = ls())

setwd("~/Desktop/tesis_sebastian/adaptive-density-weak-dep/scripts/calibration_study/Normal")

list.files(pattern = "RData")



load("Resultados_Normal_n60_phi0.RData")
Resultados_Normal_n60_phi0<-Resultados
load("Resultados_Normal_n60_phi03.RData")
Resultados_Normal_n60_phi03<-Resultados
load("Resultados_Normal_n60_phi05.RData")
Resultados_Normal_n60_phi05<-Resultados
load("Resultados_Normal_n60_phi07.RData")
Resultados_Normal_n60_phi07<-Resultados
load("Resultados_Normal_n60_phi08.RData")
Resultados_Normal_n60_phi08<-Resultados
load("Resultados_Normal_n60_phi085.RData")
Resultados_Normal_n60_phi085<-Resultados
load("Resultados_Normal_n60_phi09.RData")
Resultados_Normal_n60_phi09<-Resultados



#################################################
# Verificar que todas las mallas GAM coincidan
#################################################




#################################################
# Malla de gamma
#################################################

GAM <- Resultados_Normal_n60_phi0$GAM


#################################################
# Perdidas relativas
#################################################


PR_M_n60_phi0 <-
  Resultados_Normal_n60_phi0$MISE /
  min(Resultados_Normal_n60_phi0$MISE)

PR_M_n60_phi03 <-
  Resultados_Normal_n60_phi03$MISE /
  min(Resultados_Normal_n60_phi03$MISE)

PR_M_n60_phi05 <-
  Resultados_Normal_n60_phi05$MISE /
  min(Resultados_Normal_n60_phi05$MISE)

PR_M_n60_phi07 <-
  Resultados_Normal_n60_phi07$MISE /
  min(Resultados_Normal_n60_phi07$MISE)

PR_M_n60_phi08 <-
  Resultados_Normal_n60_phi08$MISE /
  min(Resultados_Normal_n60_phi08$MISE)

PR_M_n60_phi085 <-
  Resultados_Normal_n60_phi085$MISE /
  min(Resultados_Normal_n60_phi085$MISE)

PR_M_n60_phi09 <-
  Resultados_Normal_n60_phi09$MISE /
  min(Resultados_Normal_n60_phi09$MISE)


#################################################
# Perdida Relativa Promedio (PRP)
#################################################

PRP <- (
  PR_M_n60_phi0 +
    PR_M_n60_phi05 +
    PR_M_n60_phi07 +
    PR_M_n60_phi08 +
    PR_M_n60_phi085 +
    PR_M_n60_phi09
) / 6


#################################################
# Calibracion global
#################################################

indice.opt <- which.min(PRP)

gamma.opt <- GAM[indice.opt]

PRP.min <- min(PRP)

cat("\n")
cat("Gamma optimo global =", gamma.opt, "\n")
cat("PRP minima =", PRP.min, "\n")


#################################################
# Tabla resumen
#################################################

Tabla_PRP <- data.frame(
  GAM      = GAM,
  PR_M_n60_phi0  = PR_M_n60_phi0,
  PR_M_n60_phi0  = PR_M_n60_phi03,
  PR_M_n60_phi05 = PR_M_n60_phi05,
  PR_M_n60_phi07 = PR_M_n60_phi07,
  PR_M_n60_phi08 = PR_M_n60_phi08,
  PR_M_n60_phi085 = PR_M_n60_phi085,
  PR_M_n60_phi09 = PR_M_n60_phi09,
  PRP      = PRP
)

print(Tabla_PRP)


#################################################
# Region estable (5%)
#################################################

indice.estable <-
  which(
    PRP <= 1.05 * min(PRP)
  )

GAM.estables <- GAM[indice.estable]

cat("\n")
cat("Valores estables de gamma:\n")
print(GAM.estables)


#################################################
# Grafico de las PR
#################################################

# ind.phi0  <- which.min(PR_M_n60_phi0)
ind.phi03  <- which.min(PR_M_n60_phi03)
ind.phi05 <- which.min(PR_M_n60_phi05)
# ind.phi07 <- which.min(PR_M_n60_phi07)
# ind.phi08 <- which.min(PR_M_n60_phi08)
# ind.phi085 <- which.min(PR_M_n60_phi085)
ind.phi09 <- which.min(PR_M_n60_phi09)

# gam.phi0  <- GAM[ind.phi0]
gam.phi03  <- GAM[ind.phi03]
gam.phi05 <- GAM[ind.phi05]
# gam.phi07 <- GAM[ind.phi07]
# gam.phi08 <- GAM[ind.phi08]
# gam.phi085 <- GAM[ind.phi085]
gam.phi09 <- GAM[ind.phi09]

# min.phi0  <- min(PR_M_n60_phi0)
min.phi03  <- min(PR_M_n60_phi03)
min.phi05 <- min(PR_M_n60_phi05)
# min.phi07 <- min(PR_M_n60_phi07)
# min.phi08 <- min(PR_M_n60_phi08)
# min.phi085 <- min(PR_M_n60_phi085)
min.phi09 <- min(PR_M_n60_phi09)

# max_phi <- max(gam.phi0, gam.phi05, gam.phi07, gam.phi08, gam.phi085, gam.phi09)
max_phi <- max(gam.phi03, gam.phi05, gam.phi09)

# 1. Matplot con Zoom y formato limpio
matplot(
  GAM,
  cbind(
    # PR_M_n60_phi0,
    PR_M_n60_phi03,
    PR_M_n60_phi05,
    # PR_M_n60_phi07,
    # PR_M_n60_phi08,
    # PR_M_n60_phi085,
    PR_M_n60_phi09
  ),
  type = "l",           # Solo líneas ("l") para eliminar el ruido de los puntos solapados
  lty  = 1:6,           # Tipos de línea distintos (sólida, punteada, guiones...)
  col  = 1:3,           # Colores distintos
  lwd  = 2,             # Grosor de línea
  xlim = c(0, max_phi + 0.05),    # ZOOM EJE X: Cortamos en 0.15 para ver la separación de los mínimos
  ylim = c(0.9, 2.5),   # ZOOM EJE Y: Enfocado en la base donde ocurren las caídas
  xlab = expression(gamma),
  ylab = "Pérdida Relativa (PR)",
  main = "",            # APA 7: El título va en el documento de texto, no dentro del gráfico
  bty  = "l"            # APA 7: Ejes en forma de "L", sin recuadro superior ni derecho
)

# 2. Leyenda movida al espacio vacío
legend(
  "topleft",            # Movida a la izquierda, donde las curvas nacen y dejan espacio libre arriba
  legend = c(
    # expression(paste(phi == 0)),
    expression(paste(phi == 0.3)),
    expression(paste(phi == 0.5)),
    # expression(paste(phi == 0.7)),
    # expression(paste(phi == 0.8)),
    # expression(paste(phi == 0.85)),
    expression(paste(phi == 0.9))
  ),
  col = 1:3,
  lty = 1:5,            # Coincide con los tipos de línea del matplot
  lwd = 2,
  cex = 0.9,
  bty = "n"             # Sin caja alrededor de la leyenda
)

# 3. Líneas guía (suavizadas para no distraer)
# abline(h = min.phi0,  col = 1, lty = 3, lwd = 1)
abline(h = min.phi03,  col = 1, lty = 3, lwd = 1)
abline(h = min.phi05, col = 2, lty = 3, lwd = 1)
# abline(h = min.phi07, col = 3, lty = 3, lwd = 1)
# abline(h = min.phi08, col = 4, lty = 3, lwd = 1)
# abline(h = min.phi085, col = 5, lty = 3, lwd = 1)
abline(h = min.phi09, col = 3, lty = 3, lwd = 1)

# abline(v = gam.phi0,  col = 1, lty = 3, lwd = 1)
abline(v = gam.phi03,  col = 1, lty = 3, lwd = 1)
abline(v = gam.phi05, col = 2, lty = 3, lwd = 1)
# abline(v = gam.phi07, col = 3, lty = 3, lwd = 1)
# abline(v = gam.phi08, col = 4, lty = 3, lwd = 1)
# abline(v = gam.phi085, col = 5, lty = 3, lwd = 1)
abline(v = gam.phi09, col = 3, lty = 3, lwd = 1)

# 4. Puntos marcando exclusivamente los mínimos
# points(gam.phi0,  min.phi0,  pch = 19, col = 1, cex = 1.5)
points(gam.phi03,  min.phi03,  pch = 19, col = 1, cex = 1.5)
points(gam.phi05, min.phi05, pch = 19, col = 2, cex = 1.5)
# points(gam.phi07, min.phi07, pch = 19, col = 3, cex = 1.5)
# points(gam.phi08, min.phi08, pch = 19, col = 4, cex = 1.5) 
# points(gam.phi085, min.phi085, pch = 19, col = 5, cex = 1.5) 
points(gam.phi09, min.phi09, pch = 19, col = 3, cex = 1.5)



#################################################
# Grafico de la PRP
#################################################

#################################################
# Grafico de la PRP (Mejorado para APA 7)
#################################################

plot(
  GAM,
  PRP,
  type = "l",           # Cambiamos a "l" (línea) para eliminar los puntos solapados
  lwd  = 2,             # Línea más gruesa para mejor visualización
  col  = "black",
  xlim = c(0, 0.5),    # ZOOM EJE X: Para ver el comportamiento real del valle
  ylim = c(1, 2.5),     # ZOOM EJE Y: Recortamos los valores extremos para enfocarnos en el mínimo
  xlab = expression(gamma),
  ylab = "Pérdida Relativa Promedio (PRP)",
  main = "",            # APA 7 exige que el título vaya en el caption del documento, no en la imagen
  bty  = "l"            # APA 7: Ejes en L, sin caja superior ni derecha
)

# Línea base del mínimo (roja)
abline(
  h = PRP.min,
  col = "red",
  lty = 2,
  lwd = 1.5
)

# Línea vertical del gamma óptimo (azul)
abline(
  v = gamma.opt,
  col = "blue",
  lty = 2,
  lwd = 1.5
)

# Punto marcando exactamente el mínimo
points(
  gamma.opt,
  PRP.min,
  pch = 19,
  col = "blue",
  cex = 1.5
)

# (Opcional) Sombrear la región estable del 5% que calculaste en tu código
rect(
  xleft = min(GAM.estables), 
  ybottom = 0, 
  xright = max(GAM.estables), 
  ytop = 1.05 * PRP.min, 
  col = rgb(0.1, 0.1, 0.1, alpha = 0.1), # Gris transparente
  border = NA
)


#################################################
# Guardar resultados
#################################################
setwd("~/Desktop/tesis_sebastian/adaptive-density-weak-dep/scripts/calibration_study/PRP_Normal")

Resultados_PRP_Normal_n60 <- list(
  GAM = GAM,
  # PR_M_n60_phi0 = PR_M_n60_phi0,
  PR_M_n60_phi03 = PR_M_n60_phi03,
  PR_M_n60_phi05 = PR_M_n60_phi05,
  # PR_M_n60_phi07 = PR_M_n60_phi07,
  # PR_M_n60_phi08 = PR_M_n60_phi08,
  # PR_M_n60_phi085 = PR_M_n60_phi085,
  PR_M_n60_phi09 = PR_M_n60_phi09,
  PRP = PRP,
  gamma.opt = gamma.opt,
  PRP.min = PRP.min,
  GAM.estables = GAM.estables,
  gam.phi03 = gam.phi03,
  gam.phi05 = gam.phi05,
  # gam.phi07 = gam.phi07,
  # gam.phi08 = gam.phi08,
  # gam.phi085 = gam.phi085,
  gam.phi09 = gam.phi09
)

save(
  Resultados_PRP_Normal_n60,
  file = "Resultados_PRP_Normal_n60.RData"
)