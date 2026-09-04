###############################################################################
# Script: Generación Masiva de Gráficos MSE Puntual para Anexos
# Autor: Sebastian Zabala
# Propósito: Generar los 15 gráficos restantes excluyendo los 3 principales
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

# Convertir a formato largo y ordenar factores topológicamente
df_long_mse <- df_mse_global %>%
  pivot_longer(
    cols = c("ROT", "UCV", "SJ", "GL"), 
    names_to = "Metodo", 
    values_to = "MSE"
  ) %>%
  mutate(
    Metodo = factor(Metodo, levels = c("ROT", "UCV", "SJ", "GL")),
    Phi = as.factor(Phi),
    n = as.numeric(n), # Lo pasamos a numérico para el filtro lógico posterior
    Punto = factor(Punto, levels = c("P05", "Q1", "Mediana", "Moda", "Q3", "P95"))
  )

# =====================================================================
# 2. FUNCIÓN DE GRAFICACIÓN (PERFILES DE MSE)
# =====================================================================
colores_metodos <- c("ROT" = "#0072B2", "UCV" = "#D55E00", "SJ" = "#009E73", "GL" = "blue")

generar_grafico_perfil_mse <- function(df, nombre_densidad, tamano_n, nombre_archivo) {
  
  df_sub <- df %>% filter(Distribucion == nombre_densidad, n == tamano_n)
  
  if(nrow(df_sub) == 0) return(NULL)
  
  etiquetas_phi <- as_labeller(function(x) paste0("Dependencia (\u03d5) = ", x))
  
  p <- ggplot(df_sub, aes(x = Punto, y = MSE, color = Metodo, group = Metodo)) +
    geom_line(linewidth = 1.2, alpha = 0.8) + 
    geom_point(size = 3.5, shape = 21, fill = "white", stroke = 1.5) +
    facet_wrap(~ Phi, labeller = etiquetas_phi) + 
    theme_bw() +
    scale_y_log10() + 
    labs(
      x = "Puntos Representativos de Evaluación",
      y = expression(paste("MSE Puntual (Escala ", log[10], ")"))
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
      panel.grid.minor = element_blank() 
    )
  
  ggsave(nombre_archivo, plot = p, width = 10, height = 5.5, dpi = 300)
  cat(" -> Guardado:", nombre_archivo, "\n")
}

# =====================================================================
# 3. GENERACIÓN MASIVA CON FILTRO DE EXCLUSIÓN
# =====================================================================
# Crear nueva carpeta exclusiva para los Anexos
ruta_anexos <- file.path(ruta_resultados, "Graficos_MSE_Anexos")
if(!dir.exists(ruta_anexos)) dir.create(ruta_anexos)
setwd(ruta_anexos)

cat("\nGenerando Gráficos de Perfil MSE para Anexos...\n")

distribuciones_unicas <- unique(df_long_mse$Distribucion)
tamanos_n_unicos <- unique(df_long_mse$n)

# Contador de gráficos generados
contador <- 0

for (densidad in distribuciones_unicas) {
  for (tamano in tamanos_n_unicos) {
    
    # REGLA DE EXCLUSIÓN: Saltamos los 3 gráficos principales ya generados
    if ((densidad == "Normal" && tamano == 120) ||
        (densidad == "Lognormal" && tamano == 60) ||
        (densidad == "Modelo 10 (Garra)" && tamano == 250)) {
      next # Salta a la siguiente iteración
    }
    
    # Limpiamos el nombre para que el archivo no tenga espacios problemáticos
    nombre_limpio <- gsub(" ", "_", densidad)
    nombre_archivo <- sprintf("Perfil_MSE_Anexo_%s_n%d.png", nombre_limpio, tamano)
    
    # Generar el gráfico
    generar_grafico_perfil_mse(df_long_mse, densidad, tamano, nombre_archivo)
    contador <- contador + 1
  }
}

cat("\n=======================================================\n")
cat(sprintf("Proceso finalizado. Se generaron %d gráficos suplementarios.\n", contador))
cat("Ruta de los Anexos:\n", normalizePath(ruta_anexos), "\n")