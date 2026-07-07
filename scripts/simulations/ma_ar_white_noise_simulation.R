# ==============================================================================
# SIMULACIÓN DE PROCESOS MA(1), AR(1) Y RUIDO BLANCO
# ==============================================================================

# 1. Fijamos la semilla para asegurar que los resultados sean reproducibles
set.seed(123)

# Configuramos la ventana gráfica para ver Serie, ACF y #pacf juntos (3 filas, 1 columna)
par(mfrow = c(2, 1))

# ==============================================================================
# A) RUIDO BLANCO
# ==============================================================================
sim_rb <- arima.sim(n = 500, model = list(), sd = 1)

plot.ts(sim_rb, col = "black", main = "Ruido Blanco")
acf(sim_rb, main = "ACF - Ruido Blanco", lag=100)
#pacf(sim_rb, main = "#pacf - Ruido Blanco")

# ==============================================================================
# B) SIMULACIONES MA(1)
# ==============================================================================

# Ejemplo 1: MA(1) con theta = 0.45 y media = 0.3
sim_ma1 <- arima.sim(n = 500, model = list(ma = 0.45), sd = 0.2) + 0.3
plot.ts(sim_ma1, col = "blue", main = "MA(1): theta = 0.45, media = 0.3")
acf(sim_ma1, main = "ACF - MA(1) [0.45]", lag=100)
#pacf(sim_ma1, main = "#pacf - MA(1) [0.45]")

# Ejemplo 2: MA(1) con theta = 0.65 y media = 0.7
sim_ma2 <- arima.sim(n = 500, model = list(ma = 0.65), sd = 0.2) + 0.7
plot.ts(sim_ma2, col = "blue", main = "MA(1): theta = 0.65, media = 0.7")
acf(sim_ma2, main = "ACF - MA(1) [0.65]", lag=100)
#pacf(sim_ma2, main = "#pacf - MA(1) [0.65]")

# ==============================================================================
# C) SIMULACIONES AR(1)
# ==============================================================================
# Nota: La media teórica real de un AR(1) se calcula como: Intercepto / (1 - phi).
# Usaremos tu intercepto de 0.58 para calcular dónde se centrará la serie.

n= 500
# Ejemplo 1: AR(1) con phi = 0.25
media_ar1 <- 0.58 / (1 - 0.25)
sim_ar1 <- arima.sim(n = 500, model = list(ar = 0.25), sd = 0.10) + media_ar1
plot.ts(sim_ar1, col = "darkgreen", main = "AR(1): phi = 0.25")
acf(sim_ar1, main = "ACF - AR(1) [0.25]", lag=100, lwd=2)
#pacf(sim_ar1, main = "#pacf - AR(1) [0.25]")
plot(seq(1, n, 1), sim_ar1, xlab="", ylab="") # Diagrama de dispersion de los
# datos generados {X_i}_{i=1}^n

# Ejemplo 2: AR(1) con phi = 0.65
media_ar2 <- 0.58 / (1 - 0.65)
sim_ar2 <- arima.sim(n = 500, model = list(ar = 0.65), sd = 0.10) + media_ar2
plot.ts(sim_ar2, col = "darkgreen", main = "AR(1): phi = 0.65")
acf(sim_ar2, main = "ACF - AR(1) [0.65]", lag=100, lwd=2)
#pacf(sim_ar2, main = "#pacf - AR(1) [0.65]")

# Ejemplo 3: AR(1) con phi = 0.95
media_ar3 <- 0.58 / (1 - 0.95)
sim_ar3 <- arima.sim(n = 500, model = list(ar = 0.95), sd = 0.10) + media_ar3
plot.ts(sim_ar3, col = "darkred", main = "AR(1): phi = 0.95")
acf(sim_ar3, main = "ACF - AR(1) [0.95]", lag=100, lwd=2)
#pacf(sim_ar3, main = "#pacf - AR(1) [0.95]")

# # Seleccionando h bajo distintos métodos
# h1 <- bw.SJ(sim_ar3) # Sheather y Jones (1991)
# h2 <- bw.nrd0(sim_ar3) # Silverman(1986)
# h3 <- bw.nrd(sim_ar3)  # Scott(1992)
# h4 <- bw.ucv(sim_ar3)  # CV insesgada
# h5 <- bw.bcv(sim_ar3)  # CV sesgada  # In bw.bcv(x) : minimum occurred at one end of the range
# 
# npden1 <- density(sim_ar3, bw = h1) # 3.931768
# npden2 <- density(sim_ar3, bw = h2) # 3.847892
# npden3 <- density(sim_ar3, bw = h3) # 4.531962
# npden4 <- density(sim_ar3, bw = h4) # 4.861868
# npden5 <- density(sim_ar3, bw = h5) # 6.680812
# 
# # plot(npden)
# hist(sim_ar3, freq = FALSE, main = "Comparación de métodos clásicos para seleccionar h",
#      xlab = "Valores de x", ylab = "Densidad",
#      #xlab = paste("Bandwidth =", formatC(c(h1,h2,h3,h4,h5))), lty = 2,
#      border = "darkgray")
# #lines(c(npden1,npden2, npden3, npden4, npden5), lwd = 2)
# lines(npden1, lwd = 2, col = "blue")
# lines(npden2, lwd = 2, col = "red")
# lines(npden3, lwd = 2, col = "green")
# lines(npden4, lwd = 2, col = "darkblue")
# lines(npden5, lwd = 2, col = "darkgreen")
# rug(sim_ar3, col = "darkgray")
# 
# # Añadir leyenda
# legend("topright", 
#        legend = c("SJ", "Silverman", "Scott", "UCV", "BCV"), 
#        col = c("blue", "red", "green","darkblue", "darkgreen"), 
#        lwd = 2, 
#        bty = "n", # Sin caja alrededor de la leyenda
#        cex = 0.9)



# ==============================================================================
# 2. Restauramos la ventana gráfica a su estado original (1 gráfico por ventana)
# ==============================================================================
par(mfrow = c(1, 1))

# 1. Configurar el lienzo: 2 filas y 3 columnas
# Ajustamos los márgenes (abajo, izquierda, arriba, derecha)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# --- FILA 1: Trayectorias de los procesos AR(1) ---
plot.ts(sim_ar1, col = "darkgreen", main = "Trayectoria AR(1): phi = 0.25", 
        ylab = "Valor", xlab = "Tiempo")

plot.ts(sim_ar2, col = "darkgreen", main = "Trayectoria AR(1): phi = 0.65", 
        ylab = "Valor", xlab = "Tiempo")

plot.ts(sim_ar3, col = "darkred", main = "Trayectoria AR(1): phi = 0.95", 
        ylab = "Valor", xlab = "Tiempo")

# --- FILA 2: Correlogramas (ACF) ---
# Mantenemos lwd = 3 para que las barras sean bien visibles
acf(sim_ar1, main = "ACF [0.25]", lwd = 3, col = "darkblue", 
    ylab = "Autocorrelación")

acf(sim_ar2, main = "ACF [0.65]", lwd = 3, col = "darkblue", 
    ylab = "Autocorrelación")

acf(sim_ar3, main = "ACF [0.95]", lwd = 3, col = "darkred", 
    ylab = "Autocorrelación")

# 2. Restaurar los parámetros gráficos a su estado original (1x1)
par(mfrow = c(1, 1))


# 1. Abrir el dispositivo PDF (con height = 4.5 para evitar el error de espacio en LaTeX)
pdf("trayectorias_acf_ar1.pdf", width = 12, height = 4.5)

# 2. Configurar el lienzo y márgenes
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# --- FILA 1: Trayectorias ---
plot.ts(sim_ar1, col = "dodgerblue4", main = "Trayectoria AR(1): phi = 0.25", 
        ylab = "Valor", xlab = "Tiempo")

plot.ts(sim_ar2, col = "forestgreen", main = "Trayectoria AR(1): phi = 0.65", 
        ylab = "Valor", xlab = "Tiempo")

plot.ts(sim_ar3, col = "firebrick", main = "Trayectoria AR(1): phi = 0.95", 
        ylab = "Valor", xlab = "Tiempo")

# --- FILA 2: Correlogramas ---
acf(sim_ar1, main = "ACF [0.25]", lwd = 3, col = "darkorange3", 
    ylab = "Autocorrelación")

acf(sim_ar2, main = "ACF [0.65]", lwd = 3, col = "darkorchid4", 
    ylab = "Autocorrelación")

acf(sim_ar3, main = "ACF [0.95]", lwd = 3, col = "darkcyan", 
    ylab = "Autocorrelación")

# 3. Restaurar parámetros y cerrar el dispositivo para guardar el archivo
par(mfrow = c(1, 1))
dev.off()


# ==========================================
# === SIMULACIONES MA(2) ===
# ==========================================

# Ejemplo 1: MA(2) con theta1 = 0.6, theta2 = 0.3 y media = 0.5
sim_ma2_a <- arima.sim(n = 500, model = list(ma = c(0.6, 0.3)), sd = 0.2) + 0.5
plot.ts(sim_ma2_a, col = "darkgreen", main = "MA(2): theta = c(0.6, 0.3), media = 0.5")
acf(sim_ma2_a, main = "ACF - MA(2) [0.6, 0.3]")
#pacf(sim_ma2_a, main = "PACF - MA(2) [0.6, 0.3]")

# Ejemplo 2: MA(2) con coeficientes mixtos (theta1 = 0.8, theta2 = -0.5) y media = 0
sim_ma2_b <- arima.sim(n = 500, model = list(ma = c(0.8, -0.5)), sd = 0.2) + 0.0
plot.ts(sim_ma2_b, col = "darkgreen", main = "MA(2): theta = c(0.8, -0.5), media = 0")
acf(sim_ma2_b, main = "ACF - MA(2) [0.8, -0.5]")
#pacf(sim_ma2_b, main = "PACF - MA(2) [0.8, -0.5]")

