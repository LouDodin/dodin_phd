library(readxl)

# =====================================================
# Chemin du fichier
# =====================================================

file_path <- "bioinformatics/boxplot_cyto_vs_videodrop/Bilan extraction ADN.xlsx"

# =====================================================
# Import des données
# =====================================================

# Diamètre médian
diam_v1 <- read_excel(file_path, sheet = "Boxplot median diameter V1")
diam_v2 <- read_excel(file_path, sheet = "Boxplot median diameter V2")

# Nombre de particules
nb_v1 <- read_excel(file_path, sheet = "Boxplot number of particles V1")
nb_v2 <- read_excel(file_path, sheet = "Boxplot number of particles V2")

# =====================================================
# Fonction de création dataframe
# =====================================================

create_df <- function(data15, data20, data26) {
  
  data15 <- na.omit(data15)
  data20 <- na.omit(data20)
  data26 <- na.omit(data26)
  
  data.frame(
    value = c(data15, data20, data26),
    temperature = factor(
      c(rep("15°C", length(data15)),
        rep("20°C", length(data20)),
        rep("26°C", length(data26))),
      levels = c("15°C", "20°C", "26°C")
    )
  )
}

# =====================================================
# DataFrames diamètre
# =====================================================

df_diam_v1 <- create_df(diam_v1$data15,
                        diam_v1$data20,
                        diam_v1$data26)

df_diam_v2 <- create_df(diam_v2$data15,
                        diam_v2$data20,
                        diam_v2$data26)

# =====================================================
# DataFrames nombre particules
# =====================================================

df_nb_v1 <- create_df(nb_v1$data15,
                      nb_v1$data20,
                      nb_v1$data26)

df_nb_v2 <- create_df(nb_v2$data15,
                      nb_v2$data20,
                      nb_v2$data26)

# =====================================================
# Ratio diamètre / nombre
# =====================================================

ratio_v1_15 <- diam_v1$data15 / nb_v1$data15
ratio_v1_20 <- diam_v1$data20 / nb_v1$data20
ratio_v1_26 <- diam_v1$data26 / nb_v1$data26

ratio_v2_15 <- diam_v2$data15 / nb_v2$data15
ratio_v2_20 <- diam_v2$data20 / nb_v2$data20
ratio_v2_26 <- diam_v2$data26 / nb_v2$data26

df_ratio_v1 <- create_df(ratio_v1_15,
                         ratio_v1_20,
                         ratio_v1_26)

df_ratio_v2 <- create_df(ratio_v2_15,
                         ratio_v2_20,
                         ratio_v2_26)

# =====================================================
# Limites Y communes
# =====================================================

y_diam_min  <- min(df_diam_v1$value, df_diam_v2$value, na.rm = TRUE)
y_diam_max  <- max(df_diam_v1$value, df_diam_v2$value, na.rm = TRUE)

y_nb_min    <- min(df_nb_v1$value, df_nb_v2$value, na.rm = TRUE)
y_nb_max    <- max(df_nb_v1$value, df_nb_v2$value, na.rm = TRUE)

y_ratio_min <- min(df_ratio_v1$value, df_ratio_v2$value, na.rm = TRUE)
y_ratio_max <- max(df_ratio_v1$value, df_ratio_v2$value, na.rm = TRUE)

colors <- c("green", "yellow", "lightblue")

# =====================================================
# FIGURE 1 : DIAMETRE
# =====================================================

par(mfrow = c(1, 2), mar = c(5, 6, 4, 2))

boxplot(value ~ temperature,
        data = df_diam_v1,
        main = "V1 - Median diameter",
        ylab = "Median diameter (nm)",
        col = colors,
        border = "gray30",
        las = 1,
        ylim = c(y_diam_min, y_diam_max))

boxplot(value ~ temperature,
        data = df_diam_v2,
        main = "V2 - Median diameter",
        ylab = "Median diameter (nm)",
        col = colors,
        border = "gray30",
        las = 1,
        ylim = c(y_diam_min, y_diam_max))

# =====================================================
# FIGURE 2 : NOMBRE DE PARTICULES
# =====================================================

par(mfrow = c(1, 2), mar = c(5, 6, 4, 2))

boxplot(value ~ temperature,
        data = df_nb_v1,
        main = "V1 - Particles number",
        ylab = "Particles number",
        col = colors,
        border = "gray30",
        las = 1,
        ylim = c(y_nb_min, y_nb_max))

boxplot(value ~ temperature,
        data = df_nb_v2,
        main = "V2 - Particles number",
        ylab = "Particles number",
        col = colors,
        border = "gray30",
        las = 1,
        ylim = c(y_nb_min, y_nb_max))

# =====================================================
# FIGURE 3 : RATIO DIAMETRE / NOMBRE
# =====================================================

par(mfrow = c(1, 2), mar = c(5, 6, 4, 2))

boxplot(value ~ temperature,
        data = df_ratio_v1,
        main = "V1 - Median diameter / Particles number",
        ylab = "Ratio",
        col = colors,
        border = "gray30",
        las = 1,
        ylim = c(y_ratio_min, y_ratio_max))

boxplot(value ~ temperature,
        data = df_ratio_v2,
        main = "V2 - Median diameter / Particles number",
        ylab = "Ratio",
        col = colors,
        border = "gray30",
        las = 1,
        ylim = c(y_ratio_min, y_ratio_max))