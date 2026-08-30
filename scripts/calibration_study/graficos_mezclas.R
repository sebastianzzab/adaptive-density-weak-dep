# =====================================================================
# SCRIPT DE VISUALIZACIÓN: IMPORTACIÓN DE EXCEL Y GRÁFICOS
# =====================================================================

# 1. Carga de paquetes necesarios
if (!require(readxl)) install.packages("readxl")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")
if (!require(dplyr)) install.packages("dplyr")

library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)

# Directorio donde guardaste los Excel
setwd("~/Desktop/tesis_sebastian/adaptive-density-weak-dep/scripts/calibration_study/Nuevo")

# =====================================================================
# 2. IMPORTACIÓN Y CONSOLIDACIÓN DE DATOS (MISE y MSE)
# =====================================================================

# A) Lectura de archivos MISE (Inyectamos la columna 'n' de forma segura con mutate)
cat("Importando archivos MISE...\n")
mise_n60  <- read_excel("mise_simulacion_mezclas_n60.xlsx") %>% mutate(n = 60)
mise_n120 <- read_excel("mise_simulacion_mezclas_n120.xlsx") %>% mutate(n = 120)
mise_n250 <- read_excel("mise_simulacion_mezclas_n250.xlsx") %>% mutate(n = 250)

# Variables con nombres muy distintos para evitar confusiones
tabla_global_mise <- bind_rows(
  mise_n60, 
  mise_n120, 
  mise_n250
)

# B) Lectura de archivos MSE Puntual (Inyectamos la columna 'n')
cat("Importando archivos MSE Puntual...\n")
mse_n60  <- read_excel("mse_simulacion_mezclas_n60.xlsx") %>% mutate(n = 60)
mse_n120 <- read_excel("mse_simulacion_mezclas_n120.xlsx") %>% mutate(n = 120)
mse_n250 <- read_excel("mse_simulacion_mezclas_n250.xlsx") %>% mutate(n = 250)

tabla_puntual_mse <- bind_rows(
  mse_n60, 
  mse_n120, 
  mse_n250
)

# =====================================================================
# 3. FUNCIONES DE GRAFICACIÓN (ggplot2)
# =====================================================================

colores_metodos <- c("ROT" = "#0072B2", "UCV" = "#D55E00", "SJ" = "#009E73", "GL" = "blue")

# Función 1: Gráfico del MISE
generar_grafico_mise <- function(df_mise, nombre_densidad) {
  
  df_sub <- df_mise %>% filter(Distribucion == nombre_densidad)
  
  df_long <- df_sub %>%
    pivot_longer(cols = c(ROT, UCV, SJ, GL), names_to = "Metodo", values_to = "MISE")
  
  df_long$MISE <- as.numeric(df_long$MISE)
  
  p <- ggplot(df_long, aes(x = as.factor(Phi), y = MISE, color = Metodo, group = Metodo)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    facet_wrap(~ n, labeller = label_both) + 
    theme_bw() +
    labs(
      title = paste("MISE Global vs Dependencia y Tamaño de Muestra"),
      subtitle = paste("Distribución:", nombre_densidad),
      x = expression(phi ~ "(Nivel de Dependencia AR(1))"),
      y = "Error Integrado Cuadrático Medio (MISE)"
    ) +
    scale_color_manual(values = colores_metodos) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 11)
    )
  
  print(p)
  
  nombre_archivo <- paste0("Grafico_MISE_", gsub(" ", "_", nombre_densidad), ".png")
  ggsave(nombre_archivo, plot = p, width = 10, height = 5, dpi = 300)
  cat("Guardado:", nombre_archivo, "\n")
}

# Función 2: Gráfico del MSE Puntual
generar_grafico_mse_puntual <- function(df_mse, nombre_densidad) {
  
  df_sub <- df_mse %>% filter(Distribucion == nombre_densidad)
  
  df_long <- df_sub %>%
    pivot_longer(cols = c(ROT, UCV, SJ, GL), names_to = "Metodo", values_to = "MSE")
  
  df_long$MSE <- as.numeric(df_long$MSE)
  
  orden_puntos <- c("P05", "Q1", "Moda", "Mediana", "Q3", "P95")
  df_long$Punto <- factor(df_long$Punto, levels = orden_puntos)
  
  p <- ggplot(df_long, aes(x = as.factor(Phi), y = MSE, color = Metodo, group = Metodo)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(Punto ~ n, scales = "free_y", labeller = label_both) + 
    theme_bw() +
    labs(
      title = paste("MSE Puntual vs Dependencia y Tamaño de Muestra"),
      subtitle = paste("Distribución:", nombre_densidad),
      x = expression(phi),
      y = "Error Cuadrático Medio (MSE)"
    ) +
    scale_color_manual(values = colores_metodos) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "gray90"),
      strip.text.y = element_text(angle = 0, face = "bold"), 
      strip.text.x = element_text(face = "bold")
    )
  
  print(p)
  
  nombre_archivo <- paste0("Grafico_MSE_", gsub(" ", "_", nombre_densidad), ".png")
  ggsave(nombre_archivo, plot = p, width = 10, height = 12, dpi = 300)
  cat("Guardado:", nombre_archivo, "\n")
}

# =====================================================================
# 4. EJECUCIÓN DE GRÁFICOS
# =====================================================================
cat("\nGenerando gráficos para Mezcla Bimodal...\n")
generar_grafico_mise(tabla_global_mise, "Mezcla Bimodal")
generar_grafico_mse_puntual(tabla_puntual_mse, "Mezcla Bimodal")

cat("\nGenerando gráficos para Modelo 10 (Garra)...\n")
generar_grafico_mise(tabla_global_mise, "Modelo 10 (Garra)")
generar_grafico_mse_puntual(tabla_puntual_mse, "Modelo 10 (Garra)")