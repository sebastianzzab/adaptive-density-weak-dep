###############################################################################
# Script: Visualización de Cajas y Bigotes (Boxplots) del ISE
# Autor: Sebastian Zabala
# Propósito: Evaluar la variabilidad y estabilidad de los selectores de ventana
###############################################################################

# Limpiar entorno
rm(list = ls())

# 1. Carga de paquetes
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")

library(ggplot2)
library(dplyr)
library(tidyr)

# =====================================================================
# 1. CONFIGURACIÓN Y LECTURA AUTOMÁTICA DE DATOS
# =====================================================================
# Define la ruta donde están tus archivos anova_*.RData
ruta_resultados <- "C:/Users/sebas/Desktop/TG-Seminario/adaptive-density-weak-dep/scripts/simulations/resultados_simulacion"
setwd(ruta_resultados)

cat("Buscando archivos de datos crudos (ANOVA) en el directorio...\n")
archivos_anova <- list.files(pattern = "^anova_.*\\.RData$", full.names = TRUE)

if(length(archivos_anova) == 0) {
  stop("No se encontraron archivos .RData que comiencen con 'anova_'. Verifica la ruta.")
}

# Ensamblaje masivo de todos los datos
df_ise_global <- data.frame()

for (archivo in archivos_anova) {
  # Cargamos en un entorno aislado para no sobrescribir variables
  entorno_temp <- new.env()
  load(archivo, envir = entorno_temp)
  
  # Extraemos la lista (busca cualquier objeto tipo lista que contenga los datos)
  nombres_obj <- ls(entorno_temp)
  lista_datos <- get(nombres_obj[1], envir = entorno_temp)
  
  # Extraemos las matrices y las convertimos en Data Frame
  for (escenario in lista_datos) {
    df_temp <- as.data.frame(escenario$Matriz_Datos)
    df_temp$n <- escenario$n
    df_temp$Phi <- escenario$Phi
    df_temp$Distribucion <- escenario$Distribucion
    
    df_ise_global <- bind_rows(df_ise_global, df_temp)
  }
}

# Convertimos a formato largo (Tidy Data) para ggplot2
df_long <- df_ise_global %>%
  pivot_longer(
    cols = c("Silverman", "UCV", "SJ", "GL"), 
    names_to = "Metodo", 
    values_to = "ISE"
  ) %>%
  # Aseguramos el orden de los factores para que GL salga de último en el gráfico
  mutate(
    Metodo = factor(Metodo, levels = c("Silverman", "UCV", "SJ", "GL"), labels = c("ROT", "UCV", "SJ", "GL")),
    Phi = as.factor(Phi),
    n = as.factor(n)
  )

cat("Datos ensamblados exitosamente. Total de iteraciones registradas:", nrow(df_long) / 4, "\n")

# =====================================================================
# 2. FUNCIÓN DE GRAFICACIÓN (BOXPLOTS)
# =====================================================================
colores_metodos <- c("ROT" = "#0072B2", "UCV" = "#D55E00", "SJ" = "#009E73", "GL" = "blue")

generar_boxplot_ise <- function(df, nombre_densidad, tamano_n) {
  
  # Filtramos los datos exactos que queremos plotear
  df_sub <- df %>% filter(Distribucion == nombre_densidad, n == tamano_n)
  
  if(nrow(df_sub) == 0) return(NULL) # Si no hay datos, saltar
  
  # Etiquetador personalizado para los paneles (Facets)
  etiquetas_phi <- as_labeller(function(x) paste0("Dependencia (Phi) = ", x))
  
  p <- ggplot(df_sub, aes(x = Metodo, y = ISE, fill = Metodo)) +
    # outlier.alpha rebaja la opacidad de los puntos atípicos para que no manchen todo el gráfico
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.2, outlier.color = "gray30") +
    facet_wrap(~ Phi, labeller = etiquetas_phi, scales = "fixed") + 
    theme_bw() +
    # ESCALA LOGARÍTMICA: Fundamental para el ISE debido a la asimetría positiva extrema
    scale_y_log10() + 
    labs(
      # title = paste("Variabilidad del Error Integrado Cuadrático (ISE) - Escala Log10"),
      # subtitle = paste("Distribución:", nombre_densidad, "| Tamaño de Muestra: n =", tamano_n),
      x = "Selector de Ancho de Ventana",
      y = expression(paste("ISE (Escala ", log[10], ")"))
    ) +
    scale_fill_manual(values = colores_metodos) +
    theme(
      legend.position = "none", # Ocultamos leyenda porque el eje X ya dice el método
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold")
    )
  
  # Nombre del archivo dinámico
  nombre_archivo <- sprintf("Boxplot_ISE_%s_n%s.png", gsub(" ", "_", nombre_densidad), tamano_n)
  ggsave(nombre_archivo, plot = p, width = 10, height = 6, dpi = 300)
  cat(" -> Guardado:", nombre_archivo, "\n")
}

# =====================================================================
# 3. GENERACIÓN EN BUCLE DE TODOS LOS GRÁFICOS
# =====================================================================
distribuciones_unicas <- unique(df_long$Distribucion)
tamanos_n_unicos <- unique(df_long$n)


# Crear subcarpeta para que los gráficos no desordenen tu directorio principal
ruta_graficos <- file.path(ruta_resultados, "Graficos_Boxplots")
if(!dir.exists(ruta_graficos)) dir.create(ruta_graficos)
setwd(ruta_graficos)

cat("\nGenerando Boxplots...\n")
for (densidad in distribuciones_unicas) {
  for (tamano in tamanos_n_unicos) {
    generar_boxplot_ise(df_long, densidad, tamano)
  }
}

cat("\n=======================================================\n")
cat("Todos los Boxplots han sido generados en:\n", normalizePath(ruta_graficos), "\n")