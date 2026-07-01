# Generacion de numeros aleatorios con distribucion normal truncada [-c,c]
# Fecha: 30/06/2026
# Por: Sebastian Zabala
# Ubicacion: IVIC

# Paquete para comparar luego los numeros generados con los teoricos
# install.packages("truncnorm")
library(truncnorm)

# Estableciendo semilla para reproducibilidad
set.seed(1234)

nsim <- 1000
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

