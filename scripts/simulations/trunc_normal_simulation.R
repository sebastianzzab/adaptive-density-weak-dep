# Generacion de numeros aleatorios con distribucion normal truncada [-c,c]
# Fecha: 30/06/2026
# Por: Sebastian Zabala
# Ubicacion: IVIC

# Paquete para comparar luego los numeros generados con los teoricos
# install.packages("truncnorm")
library(truncnorm)

# Estableciendo semilla para reproducibilidad
set.seed(1234)

nsim <- 10000
U <- runif(nsim)
# U = pnorm(Z, mean = 0, sd = 1) # caso dependiente

# U=F(Z), donde Z~N(0, 1/(1-0.75^2))
# y F es su funcion de distribucion.
# U~U(0,1) y es independiente.

c<- 2; # Soporte para la normal truncada, independiente. En este caso igual a 2
p<- pnorm(c) - pnorm(-c) # Area bajo la normal (0,1) truncada
# en c, tal densidad se denota por g
# y ademas g=G'.

x<- qnorm(p*U+pnorm(-c), mean = 0, sd = 1) 
# X=G^(-1)(U)=(G^(-1)oF)(Z).
# Normal truncada dependiente.

# Histograma 
hist(x, breaks = "FD", freq = FALSE, col = "lightblue", 
     main = "Histograma vs Normal Truncada(-2,2) Teórica",
     xlab = "Valores de X", ylab = "Densidad")

curve(dtruncnorm(x, a=-c, b=c, mean = 0, sd = 1), col = "darkred", add = TRUE)

legend("topright", legend = c("Simulado (Hist)", "Teorico (Funcion)"), 
       col = c("lightblue", "red"), lwd = 2, bty = "n")
# 
# # Ajuste de leyenda para reflejar el color del histograma (fill) y la línea (col)
# legend("topright", legend = c("Simulado (Hist)", "Teórico (Función)"), 
#        fill = c("lightblue", NA), border = c("black", NA),
#        col = c(NA, "darkred"), lwd = c(NA, 2), bty = "n")


# Prueba de Kolmogorov-Smirnov (KS)
# ------------------------------------------------------
# H_0: Los datos siguen una distribucion Normal Truncada (-c,c)
ks.test(x, "ptruncnorm", a = -2, b = 2, mean = 0, sd = 1)

# Prueba de Ljung-Box (independencia)
#------------------------------------------------------
# H_0: Los datos son independientes (ruido blanco)
Box.test(x, type = "Ljung-Box")

# Gráfico de dispersion retardado
plot(x[-length(x)], x[-1], xlab = "X_t", ylab = "X_t+1")

# Correlacion
cor(x[-length(x)], x[-1]) # -0.03698361

gnormtrunc <- function(nsim, a, b, m, de) {
  # Generar variables uniformes
  U <- runif(nsim) 
  # Calcular probabilidades acumuladas en los límites a y b con los parámetros dados
  F_a <- pnorm(a, mean = m, sd = de)
  F_b <- pnorm(b, mean = m, sd = de)
  # Diferencia de probabilidades (área truncada)
  p <- F_b - F_a
  # Inversa de la función de distribución (qnorm) ajustada a la escala
  x <- qnorm(p * U + F_a, mean = m, sd = de) 
  return(x)
}

# =====================================================================
# CASO 1: Normal Estándar truncada en el centro [-1, 1]
# =====================================================================
a1 <- -1; b1 <- 1; m1 <- 0; de1 <- 1
x1 <- gnormtrunc(nsim = 10000, a = a1, b = b1, m = m1, de = de1)

# Histograma simulado
hist(x1, breaks = 50, col = "lightblue", prob = TRUE, border = "white",
     main = "Simulación vs Teoría: N(0,1) en [-1, 1]", xlab = "Valor")

# Curva teórica
curve(dtruncnorm(x, a = a1, b = b1, mean = m1, sd = de1), 
      col = "darkred", lwd = 2, add = TRUE)

# Leyenda
legend("topright", legend = c("Simulado (Hist)", "Teórico (Función)"), 
       col = c("lightblue", "darkred"), lwd = 2, bty = "n", fill=c("lightblue", NA), border=c("black", NA))

ks.test(x1, "ptruncnorm", a = a1, b = b1, mean = m1, sd = de1) # p-value = 0.2087

# =====================================================================
# CASO 2: Normal Personalizada truncada N(100, 15) en [85, 130]
# =====================================================================
a2 <- 85; b2 <- 130; m2 <- 100; de2 <- 15
x2 <- gnormtrunc(nsim = 10000, a = a2, b = b2, m = m2, de = de2)

# Histograma simulado
hist(x2, breaks = 50, col = "lightblue", prob = TRUE, border = "white",
     main = "Simulación vs Teoría: N(100,15) en [85, 130]", xlab = "Valor")

# Curva teórica
curve(dtruncnorm(x, a = a2, b = b2, mean = m2, sd = de2), 
      col = "darkred", lwd = 2, add = TRUE)

# Leyenda
legend("topright", legend = c("Simulado (Hist)", "Teórico (Función)"), 
       col = c("lightblue", "darkred"), lwd = 2, bty = "n", fill=c("lightblue", NA), border=c("black", NA))

ks.test(x2, "ptruncnorm", a = a2, b = b2, mean = m2, sd = de2) # p-value = 0.5771
# =====================================================================
# CASO 3: Extremo (Cola Derecha) N(0, 1) en [2, 5]
# =====================================================================
a3 <- 2; b3 <- 5; m3 <- 0; de3 <- 1
x3 <- gnormtrunc(nsim = 10000, a = a3, b = b3, m = m3, de = de3)

# Histograma simulado
hist(x3, breaks = 50, col = "lightblue", prob = TRUE, border = "white",
     main = "Simulación vs Teoría: N(0,1) en [2, 5]", xlab = "Valor")

# Curva teórica
curve(dtruncnorm(x, a = a3, b = b3, mean = m3, sd = de3), 
      col = "darkred", lwd = 2, add = TRUE)

# Leyenda
legend("topright", legend = c("Simulado (Hist)", "Teórico (Función)"), 
       col = c("lightblue", "darkred"), lwd = 2, bty = "n", fill=c("lightblue", NA), border=c("black", NA))

ks.test(x3, "ptruncnorm", a = a3, b = b3, mean = m3, sd = de3) # p-value = 0.5641


# Determinar que g(x) es una f.d valida
# g(x) es la f.d de una normal(0,1) truncada en (0,1)
g.x <- function(x, a0, b0) {
  p <- pnorm(b0) - pnorm(a0)
  ifelse(x < a0 | x > b0, 0, dnorm(x)) / p
}

a <- integrate(f = g.x, lower = 0, upper = 1, a0 = 0, b0 = 1)$value
a # 1