# Generacion de numeros aleatorios con distribucion normal truncada [-c,c]
# Fecha: 30/06/2026
# Por: Sebastian Zabala
# Ubicacion: IVIC

# Paquete para comparar luego los numeros generados con los teoricos
# install.packages("truncnorm")
library(truncnorm)

# Estableciendo semilla para reproducibilidad
set.seed(1234)

nsim <- 10^5
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

hist(x, breaks = "FD", freq = FALSE, col = "lightblue", 
     main = "Histograma datos con Normal Truncada(-2,2)") # Histograma de datos generados
curve(dtruncnorm(x, a=-c, b=c), col = "darkred", add = TRUE)


# gnorm_trunc <- function(n,c) {
#   X <- 
#   U <- runif(n)
#   p<- pnorm(c) - pnorm(-c)
#   X<- qnorm(p*U[i]+pnorm(-c), mean = 0, sd = 1) 
# }