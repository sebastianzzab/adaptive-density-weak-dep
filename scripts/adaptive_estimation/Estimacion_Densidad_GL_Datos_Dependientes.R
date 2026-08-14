###############################################################################
# Script: Estimacion_Densidad_GL_Datos_Dependientes.R
# Descripción: Simulación de series de tiempo con distribuciones marginales 
#              específicas (Normal, Lognormal, Mezcla) manteniendo una 
#              dependencia temporal AR(1). Posteriormente, se estima la densidad 
#              marginal utilizando el método adaptativo local de Goldenshluger-Lepski (GL)
#              y se compara su rendimiento (vía ISE) contra métodos clásicos 
#              (Validación Cruzada, Silverman y Sheather-Jones).
###############################################################################

###########################################################
#                 Parámetros Generales
###########################################################
n <- 100               # Tamaño de la muestra
phi <- 0.9             # Coeficiente autorregresivo AR(1) (Controla la dependencia temporal)
densidad <- "lognormal" # Tipo de densidad marginal objetivo ("normal", "lognormal", "mezcla")
gamma <- 0.35           # Constante de calibración para la penalización en el método GL

###########################################################
#### Parte 0: Funciones auxiliares para la distribución Mezcla
###########################################################

# Función de Distribución Acumulada (FDA) de una mezcla bimodal gaussiana
# (50% N(-2, 1) y 50% N(2, 1))
GMix <- function(x) {
  return(0.5 * pnorm(x, mean = -2, sd = 1) + 0.5 * pnorm(x, mean = 2, sd = 1))
}

# Malla para aproximar la función cuantil (inversa de la FDA) de la mezcla.
# Se construye una sola vez por eficiencia computacional.
grid.mix.x <- seq(-8, 8, length.out = 10000)
grid.mix.u <- GMix(grid.mix.x)

# Función Cuantil aproximada numéricamente por interpolación lineal.
# Necesaria para el método de la transformada inversa.
QMix <- function(u) {
  approx(x = grid.mix.u, y = grid.mix.x, xout = u, rule = 2)$y
}

###########################################################
#### Parte 1: Generación de Uniformes dependientes U(0,1)
###########################################################

# Paso 1: Generación de un proceso latente AR(1) estacionario
Z <- arima.sim(model = list(ar = phi), n = n)

# Paso 2: Cálculo de la desviación estándar teórica del proceso marginal AR(1)
sd.Z <- sqrt(1 / (1 - phi^2))

# Paso 3: Transformada Integral de Probabilidad (Uso de Cópula Gaussiana)
# Al aplicar la FDA normal a Z, obtenemos U.
U <- pnorm(Z, mean = 0, sd = sd.Z)

# Conclusión de Parte 1: 
# 'U' ahora tiene una distribución marginal Uniforme(0,1) exacta, 
# pero hereda la fuerte dependencia temporal (autocorrelación) del proceso AR(1).

###########################################################
#### Parte 2: Generación de la muestra con densidad objetivo
###########################################################

GenerarMuestra <- function(U, densidad) {
  # Se utiliza el Método de la Transformada Inversa: X = F^-1(U)
  # Esto garantiza que X tenga la densidad marginal deseada conservando
  # la estructura de dependencia de U.
  
  if(densidad == "normal") {
    MP <- 1
    X <- qnorm(U) # Inversa Normal estándar
    gsup <- 1 / sqrt(2 * pi) # Cota superior teórica (supremo) de la densidad
    gridcal <- seq(-8, 8, length = 200) # Malla espacial para evaluar la estimación
    Delta <- diff(gridcal)[1] # Tamaño del paso de la malla (para integrales numéricas)
    
  } else if(densidad == "lognormal") {
    MP <- 2
    X <- qlnorm(U, meanlog = 0, sdlog = 0.5)
    gsup <- exp(0.125) / (0.5 * sqrt(2 * pi)) # Supremo analítico de la lognormal
    gridcal <- seq(0, qlnorm(0.999999, 0, 0.5), length = 200) 
    Delta <- diff(gridcal)[1]
    
  } else if(densidad == "mezcla") {
    MP <- 3
    X <- QMix(U) # Inversa de la mezcla (aproximada numéricamente)
    gsup <- 0.5 * dnorm(-2, -2, 1) + 0.5 * dnorm(-2, 2, 1) # Supremo (ocurre en las modas)
    gridcal <- seq(-8, 8, length = 200) 
    Delta <- diff(gridcal)[1]
    
  } else {
    stop("Error: Densidad no implementada. Elija 'normal', 'lognormal' o 'mezcla'.")
  }
  
  return(list(X = X, gsup = gsup, MP = MP, gridcal = gridcal, Delta = Delta))
}

# Ejecución de la simulación
Muestra <- GenerarMuestra(U, densidad)
X <- Muestra$X
gsup <- Muestra$gsup
MP <- Muestra$MP
gridcal <- Muestra$gridcal
Delta <- Muestra$Delta

###########################################################
#### Parte 3: Verificación visual y estadística (EDA)
###########################################################

# Comprobación de momentos básicos
cat("Media empírica:", mean(X), "\nDesviación estándar empírica:", sd(X), "\n")

# Gráfico 1: Histograma vs Densidad Teórica
hist(X, probability = TRUE, breaks = 20, main = paste("Histograma Marginal -", densidad))
if(densidad == "normal") curve(dnorm(x), add = TRUE, col = 2, lwd = 2)
if(densidad == "lognormal") curve(dlnorm(x, meanlog = 0, sdlog = 0.5), add = TRUE, col = 2, lwd = 2)
if(densidad == "mezcla") curve(0.5 * dnorm(x, -2, 1) + 0.5 * dnorm(x, 2, 1), add = TRUE, col = 2, lwd = 2)

# Gráfico 2 y 3: Serie temporal y Función de Autocorrelación (evidencia de dependencia)
plot(X, type = "l", main = paste("Serie simulada (Dependiente) -", densidad), xlab = "Tiempo", ylab = "Valor")
acf(X, lag.max = 100, main = paste("Autocorrelación (ACF) -", densidad))


###########################################################
#### Parte 4: Método Adaptativo de Goldenshluger-Lepski (GL)
###########################################################

###########################################################
# Familia de ventanas (anchos de banda) candidatos H
###########################################################
# Se construye una grilla geométrica decreciente de ventanas.
u <- seq(0, floor(log(n)) * (2/3), by = 0.1)
H <- exp(-u)

###########################################################
# Estimador de Densidad por Núcleo (KDE) Estándar Gaussiano
###########################################################
gh <- function(x, h, X) {
  n <- length(X)
  return(sum(dnorm((x - X) / h)) / (n * h))
}

###########################################################
# Estimador Convolucionado (Sobresuavizado)
# Utiliza la propiedad de convolución de densidades normales.
# El ancho de banda resultante es la hipotenusa de h y hp.
###########################################################
ghh <- function(x, h, hp, X) {
  hs <- sqrt(h^2 + hp^2)
  return(gh(x, hs, X))
}

###########################################################
# Término de Penalización V(h)
# Basado en la varianza del estimador y adaptado para series 
# dependientes. 'gamma' funciona como factor de sintonización.
###########################################################
Vh <- function(h, n, gsup, gamma) {
  delta.n <- sqrt(log(n))
  K1 <- 1
  K2 <- sqrt(1 / (2 * sqrt(pi))) # Norma L2 del núcleo gaussiano
  
  termino_principal <- sqrt(2 * gamma * gsup) * K2 * (K1 + 1) * (1 + delta.n)
  tasa_convergencia <- (log(n)^(-1/2)) / sqrt(n * h)
  
  return(termino_principal * tasa_convergencia)
}

###########################################################
# Criterio A(h,x): Discrepancia Máxima
# Mide la diferencia máxima entre el estimador auxiliar cruzado
# y el estimador base, restando la penalización.
###########################################################
Ahx <- function(x, h, H, X, n, gsup, gamma) {
  valores <- numeric(length(H))
  for(i in seq_along(H)) {
    hp <- H[i]
    # Se toma la parte positiva de la diferencia penalizada
    valores[i] <- max(abs(ghh(x, h, hp, X) - gh(x, hp, X)) - Vh(hp, n, gsup, gamma), 0)
  }
  return(max(valores))
}

###########################################################
# Selección Local del Ancho de Banda (h hat)
# Minimiza el balance entre sesgo heurístico A(h,x) y varianza V(h)
###########################################################
hx <- function(x, H, X, n, gsup, gamma) {
  criterio <- numeric(length(H))
  for(i in seq_along(H)) {
    h <- H[i]
    criterio[i] <- Ahx(x, h, H, X, n, gsup, gamma) + Vh(h, n, gsup, gamma)
  }
  
  indice_optimo <- which.min(criterio)
  hhat <- H[indice_optimo]
  
  return(list(hhat = hhat, indice = indice_optimo, 
              criterio = criterio, ghat = gh(x, hhat, X)))
}

###########################################################
# Estimación Final GL sobre la Malla de Calibración
###########################################################
GLdens <- function(gridcal, H, X, n, gsup, gamma) {
  m <- length(gridcal)
  s <- length(H)
  
  criterio <- matrix(0, nrow = m, ncol = s)
  hhat <- numeric(m)
  indice <- numeric(m)
  ghat <- numeric(m)
  
  # Selecciona el ancho de banda punto por punto (Local)
  for(i in 1:m) {
    x <- gridcal[i]
    res <- hx(x, H, X, n, gsup, gamma)
    
    criterio[i,] <- res$criterio
    indice[i] <- res$indice
    hhat[i] <- res$hhat
    ghat[i] <- res$ghat
  }
  return(list(ghat = ghat, hhat = hhat, indice = indice, criterio = criterio))
}

# Ejecutar el método GL
GLgamma <- GLdens(gridcal, H, X, n, gsup, gamma)

###############################################################################
#### Parte 5: Comparación de Métodos de Selección de Ancho de Banda (KDE)
###############################################################################

# 1. Validación Cruzada Insesgada (Unbiased Cross-Validation - UCV)
hcv <- bw.ucv(X)
gh.cv <- numeric(length(gridcal))
for(i in seq_along(gridcal)) { gh.cv[i] <- gh(gridcal[i], hcv, X) }

# 2. Regla empírica de Silverman (Rule of Thumb - NRD0)
hrot <- bw.nrd0(X)
ghrot <- numeric(length(gridcal))
for(i in seq_along(gridcal)) { ghrot[i] <- gh(gridcal[i], hrot, X) }

# 3. Método de Sheather-Jones (SJ)
hsj <- bw.SJ(X)
ghsj <- numeric(length(gridcal))
for(i in seq_along(gridcal)) { ghsj[i] <- gh(gridcal[i], hsj, X) }

# Generar la densidad real teórica para comparar
if(MP == 1) { g.real <- dnorm(gridcal) }
if(MP == 2) { g.real <- dlnorm(gridcal, meanlog = 0, sdlog = 0.5) }
if(MP == 3) { g.real <- 0.5 * dnorm(gridcal, -2, 1) + 0.5 * dnorm(gridcal, 2, 1) }

###########################################################
# Gráficos Comparativos
###########################################################

# Gráfico principal de estimación de densidades
plot(gridcal, GLgamma$ghat, type = "l", lwd = 2, col = 4,
     ylim = c(0, max(GLgamma$ghat, gh.cv, ghrot, g.real)),
     ylab = "Densidad f(x)", xlab = "x", main = "Comparación de Estimadores KDE")

lines(gridcal, gh.cv, col = 3, lwd = 2)
lines(gridcal, ghsj, col = 7, lwd = 2)
lines(gridcal, ghrot, col = 6, lwd = 2)
lines(gridcal, g.real, col = 2, lwd = 3) # Real destacada
legend("topright", legend = c("GL", "CV", "SJ", "Silverman", "Densidad Real"),
       col = c(4, 3, 7, 6, 2), lwd = 2)

# Gráfico de anchos de banda locales seleccionados por GL
plot(gridcal, GLgamma$hhat, type = "l", lwd = 2, col = "darkblue",
     ylab = "Ancho de banda h(x)", xlab = "x", 
     main = "Selección de Ancho de Banda Local (GL)")

###########################################################
# Cálculo del Error Cuadrático Integrado (ISE)
# Aproximación de la integral mediante Sumas de Riemann (f_hat - f)^2 * Delta
###########################################################

ISE.GL  <- sum((GLgamma$ghat - g.real)^2) * Delta
ISE.CV  <- sum((gh.cv - g.real)^2) * Delta
ISE.SJ  <- sum((ghsj - g.real)^2) * Delta
ISE.SIL <- sum((ghrot - g.real)^2) * Delta

# Mostrar resultados
cat("\n--- Error Cuadrático Integrado (ISE) ---\n")
cat("ISE Goldenshluger-Lepski (GL):", ISE.GL, "\n")
cat("ISE Cross-Validation (CV):    ", ISE.CV, "\n")
cat("ISE Sheather-Jones (SJ):      ", ISE.SJ, "\n")
cat("ISE Silverman (Rot):          ", ISE.SIL, "\n")

