
################################################################################
### Comparación de comunidades de aves en dos sitios de bosque tropical ########
###################      a dos altitudes distintas    ##########################
################################################################################

# Proyecto final
# Autor: Jorge Alberto Valle Cruz
# Curso: Técnicas de campo para el estudio de fauna: teoría, aplicación y análisis
# Estación de Biología Chamela, Jalisco. Del 26 de octubre al 15 de noviembre del 2025



##### PREPARACIÓN DE ARCHIVOS #####

# Cargar archivos
veg <- read.csv("input/veg.csv", row.names = 1)
com <- read.csv("input/aves.csv", row.names = 1)
# Comprobar que los nombres de fila son iguales y están en el mismo orden "TRUE"
all(rownames(com) == rownames(veg))
# Convertir 'Sitio' a factor 
str(veg)
veg$Sitio <- as.factor(veg$Sitio)
str(veg)
# Ultima comprobación de datos
str(com)
all(rownames(com) == rownames(veg)) # filas en le mismo orden, debe de dar "TRUE"



                        ######## A N A L I S I S ########

##### 1. Caracterizacion descriptiva de la altitud de los sitios #####

# Cargar librerias 

library(dplyr)
library(ggplot2)

# Paleta consistente de colores (se usará a lo largo de todo el analisis)
colores <- c("A" = "#E69F00", "B" = "#56B4E9")

# tipo numérico
veg$Altitud <- as.numeric(veg$Altitud)

# Resumen por sitio
stats_por_sitio <- veg %>%
  group_by(Sitio) %>%
  summarise(
    n = n(),
    alt_min = min(Altitud, na.rm = TRUE),
    alt_q25 = quantile(Altitud, 0.25, na.rm = TRUE),
    alt_median = median(Altitud, na.rm = TRUE),
    alt_mean = mean(Altitud, na.rm = TRUE),
    alt_q75 = quantile(Altitud, 0.75, na.rm = TRUE),
    alt_max = max(Altitud, na.rm = TRUE),
    .groups = "drop"
  )

# Extremos globales y diferencias
punto_max <- veg %>% slice_max(Altitud, n = 1, with_ties = FALSE)
punto_min <- veg %>% slice_min(Altitud, n = 1, with_ties = FALSE)
delta_alt <- punto_max$Altitud - punto_min$Altitud

sitio_media_mas_alta <- stats_por_sitio %>% slice_max(alt_mean, n = 1)

print(stats_por_sitio)
as.data.frame(stats_por_sitio)

# Etiquetas
cat(sprintf(
  "Sitio con mayor media: %s (%.1f m)\nAltura máx global: %.1f m (sitio %s)\nAltura mín global: %.1f m (sitio %s)\nDiferencia máx-mín: %.1f m\n",
  sitio_media_mas_alta$Sitio, sitio_media_mas_alta$alt_mean,
  punto_max$Altitud, as.character(punto_max$Sitio),
  punto_min$Altitud, as.character(punto_min$Sitio),
  delta_alt
))

# Diferencia: máx de Sitio B - mín de Sitio A
veg$Altitud <- as.numeric(veg$Altitud)

alt_max_B <- max(veg$Altitud[veg$Sitio == "B"], na.rm = TRUE)
alt_min_A <- min(veg$Altitud[veg$Sitio == "A"], na.rm = TRUE)
delta_Bmax_Amin <- alt_max_B - alt_min_A

cat(sprintf("Máx Sitio B: %.1f m | Mín Sitio A: %.1f m | Diferencia (Bmax - Amin): %.1f m\n",
            alt_max_B, alt_min_A, delta_Bmax_Amin))


# Graficar (boxplot + puntos + medias)

p_alt <- ggplot(veg, aes(x = Sitio, y = Altitud, color = Sitio)) +
  geom_boxplot(aes(fill = Sitio), alpha = 0.2, width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, size = 2.2, alpha = 0.85) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2, fill = "white", color = "black") +

  # Marcar extremos globales
  geom_text(
    data = punto_max,
    aes(label = sprintf("Máx: %.1f m", Altitud)),
    vjust = -0.8, size = 3, color = colores[as.character(punto_max$Sitio)], show.legend = FALSE
  ) +
  geom_text(
    data = punto_min,
    aes(label = sprintf("Mín: %.1f m", Altitud)),
    vjust = 1.6, size = 3, color = colores[as.character(punto_min$Sitio)], show.legend = FALSE
  ) +
  scale_color_manual(values = colores) +
  scale_fill_manual(values = colores) +
  labs(
    title = "Altitudes por sitio",
    x = "Sitio",
    y = "Altitud (msnm)"
  ) +
  theme_minimal()

print(p_alt)
ggsave("output/Altitudes_por_sitio.tif", plot = p_alt,
       width = 18, height = 12, units = "cm", dpi = 600, bg = "white")



##### 2. Análisis de Riqueza (Curvas de Acumulación) #####

# Cargar librerías
library(vegan)

# Crear una matriz solo con los puntos del Sitio A
comm_A <- com[veg$Sitio == "A", ]

# Crear una matriz solo con los puntos del Sitio B
comm_B <- com[veg$Sitio == "B", ]


### 2.1 Calcular las Curvas de Acumulación ###

# 'method = "random" para la no independencia de los puntos en cada sitio

# para resultados reproducibles 
set.seed(123)

sac_A <- specaccum(comm_A, method = "random", permutations = 9999)
sac_B <- specaccum(comm_B, method = "random", permutations = 9999)

### 2.2 Graficar las Curvas ###

# abrir dispositivo gráfico para guardar la figura
tiff("output/curvas_acumulacion.tif", width = 18, height = 12, units = "cm", res = 600, compression = "lzw")

# Calcular el límite máximo del eje Y para ajustar tamaño del panel gráfico
ylim_max <- max(sac_A$richness, sac_B$richness) * 1.1 # 10% de espacio extra

# Primera curva (Sitio A)
plot(sac_A,
     xlab = "Puntos de Conteo (Esfuerzo)",
     ylab = "Riqueza de Especies Acumulada",
     main = "Curvas de Acumulación de Especies (Sitio A vs B)",
     col = "#E69F00",  # Color para la línea del Sitio A
     lwd = 3,       # Grosor de la línea
     ci.type = "polygon", # I.C. del 95% como un polígono
     ci.col = adjustcolor("#E69F00", alpha.f = 0.15), # con 15% de opacidad
     ylim = c(0, ylim_max) # Fijar el límite del eje Y
)

# Añadir la segunda curva (Sitio B)
plot(sac_B,
     add = TRUE,    # agregar al gráfico anterior
     col = "#56B4E9",   
     lwd = 3,
     ci.type = "polygon",
     ci.col = adjustcolor("#56B4E9", alpha.f = 0.15)) # con 15% de opacidad


# Añadir la leyenda
legend("bottomright",
       legend = c("Sitio A", "Sitio B"),
       col = c("#E69F00", "#56B4E9"),
       lwd = 2,
       bty = "n" # Sin caja alrededor de la leyenda
)

# Cerrar el dispositivo para guardar la figura
dev.off()





##### 3. Análisis de Composición #####

# Cargar librerías
library(vegan)
library(ggplot2)
library(ggrepel)
library(ggtext) 
library(vegan3d)
library(rgl)


### 3.1 Calcular la matriz de distancias de Jaccard ###
dist_aves <- vegdist(com, method = "jaccard", binary = TRUE)

#opcional visualizar y guardar como matriz
print(as.matrix(dist_aves)[1:10, 1:10])
dist_aves_tabla <- print(as.matrix(dist_aves)[1:10, 1:10])
write.csv(dist_aves_tabla, file = "output/dist_aves_tabla.csv")


### 3.2 Calcular la Ordenación NMDS 2D ###

# Para reproducibilidad 
set.seed(123) 

nmds_resultado <- metaMDS(com, 
                          distance = "jaccard",
                          binary = TRUE,
                          autotransform = FALSE,
                          k = 2,
                          trymax = 100,
                          trace = FALSE)
                          


# Revisar el "Stress"
print(nmds_resultado)

#opcional Evaluación del NMDS
nmds_resultado$stress
vegan::stressplot(nmds_resultado)     # Shepard + línea de ajuste
summary(vegan::goodness(nmds_resultado))  # ajuste por sitio

# scores de sitios y especies 
nmds_puntos <- as.data.frame(scores(nmds_resultado, display = "sites"))
plot_data_sitios <- cbind(nmds_puntos, Sitio = veg$Sitio)

nmds_especies <- as.data.frame(scores(nmds_resultado, display = "species"))
nmds_especies$Especie <- rownames(nmds_especies)


# Especies presentes en ambos sitios (presencia/ausencia)
pres_A <- colSums(comm_A > 0) > 0
pres_B <- colSums(comm_B > 0) > 0
especies_comunes <- names(which(pres_A & pres_B))

# Coordenadas NMDS solo para especies comunes
especies_com_df <- subset(nmds_especies, Especie %in% especies_comunes)

# Colocar etiquetas de nombres científicos
map <- read.csv("input/sp_comunes.csv", stringsAsFactors = FALSE)
dic <- setNames(map$Nombre_cientifico, map$Codigo)

especies_com_df$Especie <- ifelse(
  is.na(dic[especies_com_df$Especie]),
  especies_com_df$Especie,
  unname(dic[especies_com_df$Especie])
)

# Gráfico 2D con etiquetas solo de especies comunes
plot1 <- ggplot(plot_data_sitios, aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(color = Sitio), type = "t", level = 0.95, linewidth = 1.2, alpha = 0.7) +
  geom_point(aes(color = Sitio, shape = Sitio), size = 3.5, alpha = 0.8) +
  scale_color_manual(values = c("A" = "#E69F00", "B" = "#56B4E9")) +
  ggrepel::geom_text_repel(
    data = especies_com_df,
    aes(x = NMDS1, y = NMDS2, label = Especie),
    color = "grey20", fontface = "italic", size = 3, max.overlaps = 50, segment.alpha = 0.4
  ) +
  ggtitle("NMDS 2D con especies comunes a A y B",
          subtitle = paste("Stress:", round(nmds_resultado$stress, 3))) +
  xlab("Eje NMDS 1") + ylab("Eje NMDS 2") +
  theme_minimal() + coord_equal()

ggsave("output/NMDS2D_spcomunes.tif", plot = plot1, width = 18, height = 12, units = "cm", dpi = 600, bg = "white")



# Especies únicas por sitio
pres_A <- colSums(comm_A > 0) > 0
pres_B <- colSums(comm_B > 0) > 0
esp_A_unicas <- names(which(pres_A & !pres_B))
esp_B_unicas <- names(which(pres_B & !pres_A))

# Coordenadas NMDS de especies
# (nmds_especies arriba)
df_A <- subset(nmds_especies, Especie %in% esp_A_unicas)
df_B <- subset(nmds_especies, Especie %in% esp_B_unicas)

# Mapeo códigos a nombres científicos desde CSV
map <- read.csv("input/sp_unicas.csv", stringsAsFactors = FALSE)
dic <- setNames(map$Nombre_cientifico, map$Codigo)

df_A$Etiqueta <- ifelse(is.na(dic[df_A$Especie]), df_A$Especie, unname(dic[df_A$Especie]))
df_B$Etiqueta <- ifelse(is.na(dic[df_B$Especie]), df_B$Especie, unname(dic[df_B$Especie]))

# Etiquetas con cursiva y subrayado (color del sitio)
df_A$Etiqueta_html <- sprintf("<span style='font-style:italic; text-decoration:underline;'>%s</span>", df_A$Etiqueta)
df_B$Etiqueta_html <- sprintf("<span style='font-style:italic; text-decoration:underline;'>%s</span>", df_B$Etiqueta)

#Gráfico 2D con especies únicas subrayadas por color de sitio ---
plot2 <- ggplot(plot_data_sitios, aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(color = Sitio), type = "t", level = 0.95, linewidth = 1.2, alpha = 0.7) +
  geom_point(aes(color = Sitio, shape = Sitio), size = 3.5, alpha = 0.8) +
  scale_color_manual(values = c("A" = "#E69F00", "B" = "#56B4E9")) +
  
  # Etiquetas subrayadas para cada grupo (texto y subrayado en color del sitio)
  ggtext::geom_richtext(
    data = df_A, aes(label = Etiqueta_html, x = NMDS1, y = NMDS2),
    color = "#E69F00", size = 3, label.color = NA, fill = NA
  ) +
  ggtext::geom_richtext(
    data = df_B, aes(label = Etiqueta_html, x = NMDS1, y = NMDS2),
    color = "#56B4E9", size = 3, label.color = NA, fill = NA
  ) +
  
  ggtitle("NMDS 2D con especies únicas por sitio",
          subtitle = paste("Stress:", round(nmds_resultado$stress, 3))) +
  xlab("Eje NMDS 1") + ylab("Eje NMDS 2") +
  theme_minimal() + coord_equal()


ggsave("output/NMDS2D_spunicas.tif", plot = plot2, width = 18, height = 12, units = "cm", dpi = 600, bg = "white")


### 3.3 Calcular la Ordenación NMDS 3D ###

# Para reproducibilidad 
set.seed(123) 



nmds_3D_resultado <- metaMDS(com, 
                             distance = "jaccard",
                             binary = TRUE,
                             autotransform = FALSE,
                             k = 3, 
                             trymax = 100,
                             maxit = 200,       
                             trace = FALSE)

# Nuevo "Stress" (3D), Comparar este valor con el 2D.
print(nmds_3D_resultado)


# Crear el Gráfico 3D 

# Definir mismos colores
colores <- c("A" = "#E69F00", "B" = "#56B4E9")
col_sitios <- colores[veg$Sitio]

# Abrir ventana 3D

ordirgl(nmds_3D_resultado, 
        size = 4, # Tamaño de los puntos
        col = col_sitios,
        type = "p") # 'p' de 'points'

# Añadir elipses 3D
# Dibujar los elipsoides de confianza

orglellipse(
  nmds_3D_resultado,
  groups = veg$Sitio,
  display = "sites",
  kind = "se",       
  conf = 0.95,
  col = colores[levels(veg$Sitio)],
  alpha = 0.3
)

# Añadir leyenda

legend3d("topright", legend = c("Sitio A", "Sitio B"), pch = 16,
         col = colores, cex = 1.5)


# Título con stress
title3d(main = paste0("NMDS 3D (Stress = ", round(nmds_3D_resultado$stress, 3), ")"))

# Fondo blanco 
bg3d(color = "white")


##### 4. PERMANOVA #####
library(vegan)

### 4.1 PERMANOVA para el aves vs sitio ###

# Se usa 'dist_aves' y 'veg' estén en tu entorno de R)

# Beta-dispersión (homogeneidad de dispersión entre A y B)
bd <- betadisper(dist_aves, veg$Sitio)
print("--- Betadisper ANOVA ---"); print(anova(bd))
print("--- Betadisper Permutest (9999) ---"); print(permutest(bd, permutations = 9999))


# Ejecutar el PERMANOVA probando solo el efecto de 'Sitio'
permanova_sitio <- adonis2(dist_aves ~ Sitio, 
                           data = veg, 
                           permutations = 9999)

# Ver los resultados y guardar
print(permanova_sitio)
permanova_sitio_table <- as.matrix(permanova_sitio)
write.csv(permanova_sitio_table, file = "output/permanova_sitio_table.csv")



### 4.2 PERMANOVA para la vegetación ###

# Esta parte usa 'dist_aves' y 'veg'

## Preparación de datos

# Definir la lista de variables de vegetación 
vars_veg <- intersect(c("Cobertura","Densidad","AreaBasal","Altura","Tallos"), names(veg))

# Construir la fórmula con todas las variables de vegetación
formula_veg <- as.formula(paste("dist_aves ~", paste(vars_veg, collapse = " + ")))

# Imprimir la fórmula para verificar
# (La salida debería ser: dist_aves ~ Cobertura + Densidad + AreaBasal + Altura + Tallos)
print(formula_veg)

# Para reproducibilidad 
set.seed(123) 

# Diagnóstico de dispersión 
bd <- betadisper(dist_aves, veg$Sitio)
anova(bd)
permutest(bd, permutations = 9999)

### Ejecutar el PERMANOVA 
permanova_veg_solo <- adonis2(formula_veg,
                              data = veg,
                              permutations = 9999,
                              strata = veg$Sitio, # Controla la no independencia
                              by = "margin")      # Efecto único

# Ver los resultados
print(permanova_veg_solo)
permanova_veg_solo_table <- as.matrix(permanova_veg_solo)
write.csv(permanova_veg_solo_table, file = "output/permanova_veg_solo_table.csv")


### 4.3 PERMANOVA para la vegetación y sitio ###

# Construir la fórmula con Sitio y vegetación
formula_todo <- as.formula(paste("dist_aves ~ Sitio +", paste(vars_veg, collapse = " + ")))

# Imprimir la fórmula para verificar
# (La salida debería ser: dist_aves ~ Sitio + Cobertura + Densidad + ...)
print(formula_todo)

# Diagnóstico de dispersión 
bd <- betadisper(dist_aves, veg$Sitio)
print(anova(bd)); print(permutest(bd, permutations = 9999))


# Ejecutar el PERMANOVA
permanova_todo <- adonis2(formula_todo,
                          data = veg,
                          permutations = 9999,
                          by = "margin")

# Ver los resultados y guardar
print(permanova_todo)
permanova_todo_table <- as.matrix(permanova_todo)
write.csv(permanova_todo_table, file = "output/permanova_todo_table.csv")


### 4.3 PERMANOVA para la vegetación vs vegtación ###

# Crear una lista de las variables de vegetación
vars_veg <- intersect(c("Cobertura","Densidad","AreaBasal","Altura","Tallos"), names(veg))

# Crear una matriz con esas variables 
veg_matriz <- veg[, vars_veg]

# Estandarizar las variables
# Esto convierte todo a una media de 0 y desviación estándar de 1
veg_scaled <- scale(veg_matriz)

# Crear una matriz de distancia (Euclideana) a partir de la vegetación
dist_veg <- vegdist(veg_scaled, method = "euclidean")

# Homogeneidad de dispersión
bd_env <- betadisper(dist_veg, veg$Sitio)
anova(bd_env)
permutest(bd_env, permutations = 9999)

# Ejecutar el PERMANOVA
# La fórmula prueba si la distancia de la vegetación (dist_veg) es explicada por el 'Sitio'
permanova_veg <- adonis2(dist_veg ~ Sitio, 
                         data = veg, 
                         permutations = 9999)


# Ver los resultados y guardar
print(permanova_veg)
permanova_veg_table <- as.matrix(permanova_veg)
write.csv(permanova_veg_table, file = "output/permanova_veg_table.csv")


### 4.4 Prueba de t para cada variable de la vegetación ###

# Ejecutar una prueba t para cada variable de vegetación
# Usamos 'lapply' para "iterar" o "hacer un bucle" sobre tu lista 'vars_veg'
resultados_t <- lapply(vars_veg, function(variable) {
  
  # Ejecutar la prueba t: variable ~ Sitio (desde la tabla 'veg')
  # La fórmula significa "compara la 'variable' en función del 'Sitio'"
  t_test <- t.test(veg[[variable]] ~ veg$Sitio)
  
  # Devolver un resumen con la información clave
  data.frame(
    Variable = variable,
    Media_Sitio_A = t_test$estimate[1],
    Media_Sitio_B = t_test$estimate[2],
    p_valor = t_test$p.value
  )
})

# Combinar los resultados de la lista en una sola tabla
tabla_resultados_t <- do.call(rbind, resultados_t)

# Imprimir la tabla
print(tabla_resultados_t)

#  guardar
tabla_resultados_tt <- as.matrix(tabla_resultados_t)
write.csv(tabla_resultados_tt, file = "output/t_veg_table.csv")


##### 5. ANALISIS DISCRIMINANTE #####


library(vegan)

# (Asegúrate de tener tu matriz de comunidad 'com' (aves)
# y tu dataframe 'veg' (con la columna 'Sitio') en tu entorno)

#Ejecutar el análisis SIMPER
#    Compara la comunidad 'com' agrupándola por 'veg$Sitio'
simper_resultado <- simper(com, veg$Sitio)

# 2. Imprimir el resumen
#    Esto te mostrará la lista de especies ordenadas por
#    cuánto contribuyen a la diferencia entre A y B
summary(simper_resultado)
#  guardar
simper_resultado <- as.data.frame(simper_resultado)
write.csv(tabla_resultados_tt, file = "output/t_veg_table.csv")



### El flujo se divide en responder tres preguntas clave:

###¿Cuántas especies hay? (Riqueza)

#¿Son diferentes las comunidades? (Composición)

#¿Por qué son diferentes? (Relación con el ambiente)

## 1. Análisis de Riqueza (Curvas de Acumulación)
## 2. Análisis de Composición
## 3. Análisis Explicativo (Relación con el Ambiente)





