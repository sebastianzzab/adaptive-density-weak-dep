# KDE: MÉTODOS CLÁSICOS DE SELECCIÓN DEL PARÁMETRO VENTANA (h) EN R
# Link: https://rubenfcasal.github.io/book_remuestreo/npden-r.html
# =========================================================================================================================
# R tiene varias funciones para seleccionar automáticamente el ancho de banda (h) para la estimación de densidad por núcleo.
# bw.nrd0, bw.nrd, bw.ucv, bw.bcv y bw.SJ son métodos predefinidos que implementan diferentes enfoques para seleccionar h.
# =========================================================================================================================
# bw.nrd0 implementa una regla práctica para elegir el ancho de banda de un estimador de densidad de 
# núcleo gaussiano. Por defecto, es 0,9 veces el mínimo de la desviación estándar y el rango intercuartil 
# dividido por 1,34 veces el tamaño de la muestra elevado a la potencia de menos un quinto 
# (= la "regla práctica" de Silverman, Silverman (1986, página 48, ecuación (3.31))) a menos que 
# los cuartiles coincidan, en cuyo caso se garantiza un resultado positivo.
# =========================================================================================================================
# bw.nrd es la variante más común propuesta por Scott (1992), con un factor de 1,06.
# =========================================================================================================================
# bw.ucv y bw.bcv implementan la validación cruzada insesgada y sesgada, respectivamente.
# =========================================================================================================================
# bw.SJ implementa los métodos de Sheather y Jones (1991) para seleccionar el ancho de banda mediante la 
# estimación piloto de derivadas.

# El algoritmo del método "ste" resuelve una ecuación (mediante la raíz unitaria) y, por ello, amplía el 
# intervalo c(inferior, superior) cuando los límites no han sido especificados por el usuario y no incluyen 
# la raíz.

# Los últimos tres métodos utilizan distancias agrupadas por pares: su complejidad es O(n²) hasta n = nb/2 y O(n) 
# a partir de entonces. Debido a este agrupamiento, los resultados difieren ligeramente cuando x se traslada o 
# cambia de signo.
# =========================================================================================================================

# =========================================================================================================================
# EFECTO DE LA SELECCIÓN DEL PARÁMETRO VENTANA (h) PARA LA ESTIMACIÓN DE LA  
# DENSIDAD DE LAS PRECIPITACIONES ANUALES EN CIERTAS CIUDADES ESTADOUNIDENSES
# =========================================================================================================================
# Definiendo el objeto con mis datos
x <- precip 

# Seleccionando h bajo distintos métodos
h1 <- bw.SJ(x) # Sheather y Jones (1991)
h2 <- bw.nrd0(x) # Silverman(1986)
h3 <- bw.nrd(x)  # Scott(1992)
h4 <- bw.ucv(x)  # CV insesgada
h5 <- bw.bcv(x)  # CV sesgada  #OBSERVACIÓN: In bw.bcv(x) : minimum occurred at one end of the range

# Ajustando las densidades para los distintos h obtenidos:
npden1 <- density(x, bw = h1) # h = 3.931768
npden2 <- density(x, bw = h2) # h = 3.847892
npden3 <- density(x, bw = h3) # h = 4.531962
npden4 <- density(x, bw = h4) # h = 4.861868
npden5 <- density(x, bw = h5) # h = 6.680812

# ventanas <- c(bw.SJ(x), bw.nrd0(x), bw.nrd(x), bw.ucv(x), bw.bcv(x))
# densidades <- c(density(x, bw = h1), density(x, bw = h2), density(x, bw = h3), density(x, bw = h4), density(x, bw = h5))
# colores = c("blue", "red", "green", "darkblue", "darkgreen")

hist(x, freq = FALSE, main = "Comparación de métodos clásicos para seleccionar h",
     xlab = "Valores de x", ylab = "Densidad",
     #xlab = paste("Bandwidth =", formatC(c(h1,h2,h3,h4,h5))), lty = 2,
     border = "darkgray", xlim = c(0, 80), ylim = c(0, 0.04))
#lines(c(npden1,npden2, npden3, npden4, npden5), lwd = 2)
lines(npden1, lwd = 2, col = "blue")
lines(npden2, lwd = 2, col = "red")
lines(npden3, lwd = 2, col = "green")
lines(npden4, lwd = 2, col = "darkblue")
lines(npden5, lwd = 2, col = "darkgreen")
rug(x, col = "darkgray")

# Añadir leyenda
legend("topright", 
       legend = c("SJ", "Silverman", "Scott", "UCV", "BCV"), 
       col = c("blue", "red", "green","darkblue", "darkgreen"), 
       lwd = 2, 
       bty = "n", # Sin caja alrededor de la leyenda
       cex = 0.9)