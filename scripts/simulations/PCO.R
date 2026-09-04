###########################################################
#                  Configuración de Parámetros
###########################################################
n <- 1000          # Tamaño de la muestra a simular[cite: 1]
phi <- 0.9         # Coeficiente de autocorrelación (dependencia temporal)[cite: 1]
# densidad <- "normal" # Tipo de densidad objetivo ('normal', 'lognormal' o 'mezcla')[cite: 1]
densidad <- "mezcla"
###########################################################
#### Parte 0: Funciones para la Densidad Mezcla
###########################################################

# Función de distribución acumulada (CDF) de la mezcla de dos normales[cite: 1]
GMix <- function(x)
{
  return(0.5*pnorm(x, mean=-2, sd=1) + 0.5*pnorm(x, mean=2, sd=1))
}

# Creación de una malla de puntos y sus probabilidades para invertir la CDF[cite: 1]
grid.mix.x <- seq(-8, 8, length.out=10000)
grid.mix.u <- GMix(grid.mix.x)

# Función cuantil (inversa de la CDF) aproximada mediante interpolación[cite: 1]
QMix <- function(u)
{
  approx(x = grid.mix.u, y = grid.mix.x, xout = u, rule = 2)$y
}

###########################################################
#### Parte 1: Generación de Uniformes con Dependencia
###########################################################

# Simulación de un proceso autorregresivo AR(1) para introducir dependencia[cite: 1]
Z <- arima.sim(model=list(ar=phi), n=n)

# Cálculo de la desviación estándar teórica del proceso estacionario[cite: 1]
sd.Z <- sqrt(1/(1 - phi^2))

# Transformación Integral de Probabilidad para obtener U(0,1) con dependencia[cite: 1]
U <- pnorm(Z, mean=0, sd=sd.Z)

###########################################################
#### Parte 2: Generación de la Muestra Objetivo
###########################################################

# Función para transformar las uniformes a la distribución deseada[cite: 1]
GenerarMuestra <- function(U, densidad)
{
  if(densidad=="normal")
  {
    MP <- 1
    X <- qnorm(U)
    gsup <- 1/sqrt(2*pi) 
    gridcal <- seq(-8, 8, length=200) 
    Delta <- diff(gridcal)[1]
  }
  else if(densidad=="lognormal")
  {
    MP <- 2
    X <- qlnorm(U, meanlog=0, sdlog=0.5)
    gsup <- exp(0.125)/(0.5*sqrt(2*pi))
    gridcal <- seq(0, qlnorm(0.999999, 0, 0.5), length=200)
    Delta <- diff(gridcal)[1]
  }
  else if(densidad=="mezcla")
  {
    MP <- 3
    X <- QMix(U)
    gsup <- 0.5*dnorm(-2, -2, 1) + 0.5*dnorm(-2, 2, 1)
    gridcal <- seq(-8, 8, length=200)
    Delta <- diff(gridcal)[1]
  }
  else
  {
    stop("Densidad no implementada")
  }
  
  return(list(X=X, gsup=gsup, MP=MP, gridcal=gridcal, Delta=Delta))
}

# Ejecución de la generación[cite: 1]
Muestra <- GenerarMuestra(U, densidad)
X <- Muestra$X
MP <- Muestra$MP
gridcal <- Muestra$gridcal
Delta <- Muestra$Delta

###########################################################
#### Parte 3: Análisis Descriptivo y Visualización
###########################################################

# Visualización del histograma comparado con la densidad teórica[cite: 1]
hist(X, probability = TRUE, breaks = 20, main="Histograma vs Teórica")

if(densidad=="normal") {curve(dnorm(x), add = TRUE, col = 2, lwd = 2)}
if(densidad=="lognormal") {curve(dlnorm(x, meanlog=0, sdlog=0.5), add=TRUE, col=2, lwd=2)}
if(densidad=="mezcla") {curve(0.5*dnorm(x, -2, 1) + 0.5*dnorm(x, 2, 1), add=TRUE, col=2, lwd=2)}

###########################################################
#### Parte 4: Estimación de Densidad por Método PCO 
###########################################################

# Definición de la familia de ventanas candidatas (H)[cite: 1]
u <- seq(0, floor(log(n))*(2/3), by = 0.1)
H <- exp(-u)
h_min <- min(H)

# Estimador de Núcleo (Kernel) Gaussiano básico[cite: 1]
gh <- function(x, h, X)
{
  n_len <- length(X)
  return(sum(dnorm((x-X)/h))/(n_len*h))
}

# Producto interno exacto de dos núcleos Gaussianos para el cálculo de penalidad
inner_prod_K <- function(h1, h2) {
  return(1 / sqrt(2 * pi * (h1^2 + h2^2)))
}

# Evaluación del estimador de "overfitting" en la malla
gh_min <- sapply(gridcal, function(z) gh(z, h_min, X))

# Función principal para seleccionar el ancho de banda usando PCO
PCO_Selection <- function(H, h_min, X, gridcal, Delta, gh_min) {
  n_len <- length(X)
  criterios <- numeric(length(H))
  estimadores <- matrix(0, nrow=length(gridcal), ncol=length(H))
  
  for(i in seq_along(H)) {
    h <- H[i]
    
    # 1. Calcular estimador para ventana h
    gh_h <- sapply(gridcal, function(z) gh(z, h, X))
    estimadores[, i] <- gh_h
    
    # 2. Calcular la norma L2 empírica (ISE discreto) respecto al overfitting
    norm_diff <- sum((gh_h - gh_min)^2) * Delta
    
    # 3. Calcular penalidad óptima
    pen_opt <- 2 * inner_prod_K(h, h_min) / n_len
    
    # 4. Criterio PCO
    criterios[i] <- norm_diff + pen_opt
  }
  
  # Seleccionar el estimador que minimiza el criterio
  idx_optimo <- which.min(criterios)
  
  return(list(
    hhat = H[idx_optimo],
    ghat = estimadores[, idx_optimo],
    criterios = criterios
  ))
}

# Ejecutar estimación PCO
PCO_result <- PCO_Selection(H, h_min, X, gridcal, Delta, gh_min)

###########################################################
#### Parte 5: Comparativa con Métodos Clásicos
###########################################################

# 1. Cross-Validation (Validación Cruzada)[cite: 1]
hcv <- bw.ucv(X)
gh.cv <- sapply(gridcal, function(z) gh(z, hcv, X))

# 2. Regla de Silverman[cite: 1]
h.sil <- bw.nrd0(X)
gh.sil <- sapply(gridcal, function(z) gh(z, h.sil, X))

# Determinar la densidad real para el cálculo del error (ISE)[cite: 1]
if(MP==1) {g.real <- dnorm(gridcal)}
if(MP==2) {g.real <- dlnorm(gridcal, meanlog=0, sdlog=0.5)}
if(MP==3) {g.real <- 0.5*dnorm(gridcal, -2, 1) + 0.5*dnorm(gridcal, 2, 1)}

###########################################################
#### Parte 6: Resultados y Errores (ISE)
###########################################################

# Gráfico comparativo de todas las estimaciones[cite: 1]
plot(gridcal, PCO_result$ghat, type="l", lwd=2, col=4,
     ylim=c(0, max(PCO_result$ghat, gh.cv, gh.sil, g.real)), ylab="Densidad")
lines(gridcal, gh.cv, col=3, lwd=2)
lines(gridcal, gh.sil, col=6, lwd=2)
lines(gridcal, g.real, col=2, lwd=2)
legend("topright", legend=c("PCO", "CV", "Silverman", "Real"), col=c(4, 3, 6, 2), lwd=2)

# Cálculo del Error Cuadrático Integrado Aproximado (ISE)[cite: 1]
ISE.PCO <- sum((PCO_result$ghat - g.real)^2) * Delta
ISE.CV  <- sum((gh.cv - g.real)^2) * Delta
ISE.SIL <- sum((gh.sil - g.real)^2) * Delta

print(data.frame(Metodo=c("PCO", "CV", "Silverman"), ISE=c(ISE.PCO, ISE.CV, ISE.SIL)))