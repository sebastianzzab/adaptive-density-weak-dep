# ==============================================================================
# DEPENDENCIA DEL HISTOGRAMA AL ANCLA Y AL ANCHO DE VENTANA
# Histograma como estimador de densidad (APA 7.ª ed.)
# ==============================================================================
library(ggplot2)
library(patchwork)

# 1. Simulación de datos
set.seed(2023) 
n <- 40        
datos <- data.frame(x = rnorm(n, mean = 0, sd = 1.2))

# Parámetros comunes
y_limit <- 0.5        # Límite fijo eje Y para densidad
ancla_fija <- 0        # Origen constante para gráficos de ancho

# Colores sobrios
color_barras <- "#E8E8E8"      # Gris muy claro para relleno
color_borde <- "#333333"       # Gris oscuro para bordes

# ========================
# 2. Función para ANCLA
# ========================
plot_clean <- function(data, ancla, label_text) {
  ggplot(data, aes(x = x)) +
    geom_histogram(
      binwidth = 1.0,
      boundary = ancla,
      fill = color_barras,
      color = color_borde,   
      linewidth = 0.4,
      aes(y = after_stat(density))  # Densidad en lugar de frecuencia
    ) +
    geom_rug(sides = "b", 
             alpha = 0.5, 
             linewidth = 0.3,
             color = "#666666") +
    coord_cartesian(ylim = c(0, y_limit), xlim = c(-4, 4)) +
    labs(tag = label_text, x = NULL, y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      plot.tag.position = "topleft",
      plot.tag = element_text(size = 10, face = "italic", 
                              color = "black",
                              margin = margin(b = 6, l = 2)),
      panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      axis.line = element_line(color = "#333333", linewidth = 0.4),
      axis.ticks = element_line(color = "#333333", linewidth = 0.4),
      plot.margin = margin(t = 5, r = 8, b = 5, l = 8)
    )
}

# 3. Generar 4 gráficos (ancla)
p1 <- plot_clean(datos, ancla = 0.15, label_text = "(a) Origen: 0.15") + 
  labs(y = "Densidad")
p2 <- plot_clean(datos, ancla = 0.25, label_text = "(b) Origen: 0.25") + 
  labs(y = "Densidad")
p3 <- plot_clean(datos, ancla = 0.50, label_text = "(c) Origen: 0.50") + 
  labs(y = "Densidad", x = "x")
p4 <- plot_clean(datos, ancla = 0.75, label_text = "(d) Origen: 0.75") + 
  labs(y = "Densidad", x = "x")

# Composición final
fig_ancla <- (p1 + p2) / (p3 + p4)

ggsave("histograma_ancla.png", plot = fig_ancla, width = 7, height = 5, dpi = 300)

# ========================
# 4. Función para ANCHO DE VENTANA
# ========================
plot_bandwidth <- function(data, h, label_text) {
  ggplot(data, aes(x = x)) +
    geom_histogram(
      binwidth = h,
      boundary = ancla_fija,
      fill = color_barras,
      color = color_borde,   
      linewidth = 0.4,
      aes(y = after_stat(density))  # Densidad en lugar de frecuencia
    ) +
    geom_rug(sides = "b", 
             alpha = 0.5, 
             linewidth = 0.3,
             color = "#666666") +
    coord_cartesian(xlim = c(-4, 4)) +
    labs(tag = label_text, x = NULL, y = "Densidad") +
    theme_classic(base_size = 10) +
    theme(
      plot.tag.position = "topleft",
      plot.tag = element_text(size = 10, face = "italic", 
                              color = "black",
                              margin = margin(b = 6, l = 2)),
      panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      axis.line = element_line(color = "#333333", linewidth = 0.4),
      axis.ticks = element_line(color = "#333333", linewidth = 0.4),
      plot.margin = margin(t = 5, r = 8, b = 5, l = 8)
    )
}

# 5. Generar 4 gráficos (ancho de ventana)
b1 <- plot_bandwidth(datos, h = 0.25, label_text = "(a) Ancho: 0.25")
b2 <- plot_bandwidth(datos, h = 0.50, label_text = "(b) Ancho: 0.50")
b3 <- plot_bandwidth(datos, h = 1.00, label_text = "(c) Ancho: 1.00")
b4 <- plot_bandwidth(datos, h = 2.00, label_text = "(d) Ancho: 2.00")

# Ajuste de ejes redundantes
b1 <- b1 + labs(x = "x")
b2 <- b2 + labs(x = "x")
b3 <- b3 + labs(x = "x")
b4 <- b4 + labs(x = "x")

# Composición final
fig_ancho <- (b1 + b2) / (b3 + b4)

ggsave("histograma_ancho.png", plot = fig_ancho, width = 7, height = 5, dpi = 300)
