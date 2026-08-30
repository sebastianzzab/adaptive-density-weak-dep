# =====================================================================
# SCRIPT FIGURA DE DENSIDADES (DISEÑO 3x2)
# =====================================================================

# 1. Configuración de alta resolución (Opcional)
# pdf("Figura_Densidades.pdf", width = 8, height = 10)

# 2. Configuración de la cuadrícula y márgenes
# mar = c(abajo, izquierda, arriba, derecha). 
# Se ajustó a c(3, 4, 3, 1) para dar espacio al título arriba (3 líneas) y menos espacio abajo.
par(mfrow = c(3, 2), mar = c(3, 4, 3, 1), oma = c(0, 0, 0, 0))

# ---------------------------------------------------------------------
# (a) Normal
# ---------------------------------------------------------------------
x_norm <- seq(-3, 3, length.out = 500)
y_norm <- dnorm(x_norm, mean = 0, sd = 1)
plot(x_norm, y_norm, type = "l", xlab = "", ylab = "", 
     main = "(a) Normal", las = 0, cex.main = 1.2, font.main = 1)

# ---------------------------------------------------------------------
# (b) Laplace
# ---------------------------------------------------------------------
b_param <- 1 / sqrt(2) 
x_lap <- seq(-4, 4, length.out = 500)
y_lap <- (1 / (2 * b_param)) * exp(-abs(x_lap) / b_param)
plot(x_lap, y_lap, type = "l", xlab = "", ylab = "", 
     main = "(b) Laplace", las = 0, cex.main = 1.2, font.main = 1)

# ---------------------------------------------------------------------
# (c) Lognormal
# ---------------------------------------------------------------------
x_ln <- seq(0, 4, length.out = 500)
y_ln <- dlnorm(x_ln, meanlog = 0, sdlog = 0.5)
plot(x_ln, y_ln, type = "l", xlab = "", ylab = "", 
     main = "(c) Lognormal", las = 0, cex.main = 1.2, font.main = 1)

# ---------------------------------------------------------------------
# (d) Gamma
# ---------------------------------------------------------------------
media_ln <- exp(0.125)
var_ln <- (exp(0.25) - 1) * exp(0.25)
shape_g <- (media_ln^2) / var_ln
rate_g <- media_ln / var_ln

x_gam <- seq(0, 4, length.out = 500)
y_gam <- dgamma(x_gam, shape = shape_g, rate = rate_g)
plot(x_gam, y_gam, type = "l", xlab = "", ylab = "", 
     main = "(d) Gamma", las = 0, cex.main = 1.2, font.main = 1)

# ---------------------------------------------------------------------
# (e) Mezcla Bimodal
# ---------------------------------------------------------------------
x_mix <- seq(-6, 6, length.out = 500)
y_mix <- 0.5 * dnorm(x_mix, mean = -2, sd = 1) + 0.5 * dnorm(x_mix, mean = 2, sd = 1)
plot(x_mix, y_mix, type = "l", xlab = "", ylab = "", 
     main = "(e) Mezcla Bimodal", las = 0, cex.main = 1.2, font.main = 1)

# ---------------------------------------------------------------------
# (f) Modelo 10 (Garra)
# ---------------------------------------------------------------------
x_garra <- seq(-3, 3, length.out = 500)
y_garra <- 0.5 * dnorm(x_garra, mean = 0, sd = 1)
for(j in 0:4) {
  y_garra <- y_garra + 0.1 * dnorm(x_garra, mean = (j/2 - 1), sd = 0.1)
}
plot(x_garra, y_garra, type = "l", xlab = "", ylab = "", 
     main = "(f) Modelo 10 (Garra)", las = 0, cex.main = 1.2, font.main = 1)

# 3. Restaurar parámetros originales de la ventana gráfica
par(mfrow = c(1, 1))

# dev.off() # Descomentar si usaste pdf()