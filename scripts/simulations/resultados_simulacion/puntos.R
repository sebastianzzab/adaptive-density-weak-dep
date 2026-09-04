###############################################################################
# Script: Visualización del Error Cuadrático Medio Puntual (MSE)
# Autor: Sebastian Zabala
# Propósito: Evaluar topológicamente el error en colas, centro y picos
###############################################################################

# Limpiar entorno
rm(list = ls())

# 1. Carga de paquetes
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(readxl)) install.packages("readxl")

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# =====================================================================
# 1. LECTURA Y PREPARACIÓN DE DATOS
# =====================================================================
# Define la ruta donde están tus archivos Excel de MSE puntual
ruta_resultados <- "C:/Users/sebas/Desktop/TG-Seminario/adaptive-density-weak-dep/scripts/simulations/resultados_simulacion"
setwd(ruta_resultados)

cat("Buscando archivos de MSE puntual en el directorio...\n")
archivos_mse <- list.files(pattern = "^mse_puntual_.*\\.xlsx$", full.names = TRUE)

if(length(archivos_mse) == 0) {
  stop("No se encontraron archivos que comiencen con 'mse_puntual_'. Verifica la ruta.")
}

# Leer y ensamblar todos los archivos de MSE puntual
df_mse_global <- data.frame()
for (archivo in archivos_mse) {
  df_temp <- read_excel(archivo)
  df_mse_global <- bind_rows(df_mse_global, df_temp)
}

# Convertir a formato largo para ggplot2 y ordenar los factores lógicamente
df_long_mse <- df_mse_global %>%
  pivot_longer(
    cols = c("ROT", "UCV", "SJ", "GL"), 
    names_to = "Metodo", 
    values_to = "MSE"
  ) %>%
  mutate(
    # Ordenar los métodos para coherencia con los gráficos anteriores
    Metodo = factor(Metodo, levels = c("ROT", "UCV", "SJ", "GL")),
    Phi = as.factor(Phi),
    n = as.factor(n),
    # Ordenar los puntos de evaluación de izquierda a derecha topológicamente
    Punto = factor(Punto, levels = c("P05", "Q1", "Mediana", "Moda", "Q3", "P95"))
  )

# =====================================================================
# 2. FUNCIÓN DE GRAFICACIÓN (PERFILES DE MSE)
# =====================================================================
colores_metodos <- c("ROT" = "#0072B2", "UCV" = "#D55E00", "SJ" = "#009E73", "GL" = "blue")

generar_grafico_perfil_mse <- function(df, nombre_densidad, tamano_n, nombre_archivo) {
  
  df_sub <- df %>% filter(Distribucion == nombre_densidad, n == tamano_n)
  
  if(nrow(df_sub) == 0) {
    cat("No hay datos para", nombre_densidad, "con n =", tamano_n, "\n")
    return(NULL)
  }
  
  # Etiquetador para los paneles
  etiquetas_phi <- as_labeller(function(x) paste0("Dependencia (\u03d5) = ", x))
  
  p <- ggplot(df_sub, aes(x = Punto, y = MSE, color = Metodo, group = Metodo)) +
    geom_line(linewidth = 1.2, alpha = 0.8) + 
    geom_point(size = 3.5, shape = 21, fill = "white", stroke = 1.5) +
    facet_wrap(~ Phi, labeller = etiquetas_phi) + 
    theme_bw() +
    # Escala logarítmica fundamental para no aplastar los errores del centro
    scale_y_log10() + 
    labs(
      x = "Puntos Representativos de Evaluación",
      y = expression(paste("MSE Puntual (Escala ", log[10], ")"))
      # Títulos y subtítulos eliminados por Norma APA 7
    ) +
    scale_color_manual(values = colores_metodos) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
      axis.title = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank() # Limpiar líneas de cuadrícula menores
    )
  
  ggsave(nombre_archivo, plot = p, width = 10, height = 5.5, dpi = 300)
  cat(" -> Guardado:", nombre_archivo, "\n")
}

# =====================================================================
# 3. GENERACIÓN DE LOS 3 GRÁFICOS ESTRATÉGICOS
# =====================================================================
ruta_graficos <- file.path(ruta_resultados, "Graficos_MSE")
if(!dir.exists(ruta_graficos)) dir.create(ruta_graficos)
setwd(ruta_graficos)

cat("\nGenerando Gráficos de Perfil Topológico del MSE...\n")

# CASO 1: Regiones de Alta Densidad (Centro) -> Normal, n = 120
generar_grafico_perfil_mse(df_long_mse, "Normal", 120, "Perfil_MSE_Normal_n120.png")

# CASO 2: Zonas de Escasez de Datos (Colas) -> Lognormal, n = 60
generar_grafico_perfil_mse(df_long_mse, "Lognormal", 60, "Perfil_MSE_Lognormal_n60.png")

# CASO 3: Estructuras Multimodales Fina (Picos) -> Modelo 10 (Garra), n = 250
generar_grafico_perfil_mse(df_long_mse, "Modelo 10 (Garra)", 250, "Perfil_MSE_Garra_n250.png")

cat("\n=======================================================\n")
cat("Los 3 gráficos estratégicos han sido generados en:\n", normalizePath(ruta_graficos), "\n")