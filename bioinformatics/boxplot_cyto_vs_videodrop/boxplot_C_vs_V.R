library(readxl)

data <- read_excel(
  "bioinformatics/boxplot_cyto_vs_videodrop/Bilan extraction ADN.xlsx",
  sheet = "Compare cyto videodrop"
)

# -----------------------
# CYTOMETER
# -----------------------
cyto15 <- na.omit(data$`cyto 15`)
cyto20 <- na.omit(data$`cyto 20`)
cyto26 <- na.omit(data$`cyto 26`)

df_cyto <- data.frame(
  concentration = c(cyto15, cyto20, cyto26),
  temperature = factor(rep(c("15°C", "20°C", "26°C"),
                           times = c(length(cyto15),
                                     length(cyto20),
                                     length(cyto26))))
)

# -----------------------
# VIDEODROP
# -----------------------
videodrop15 <- na.omit(data$`videodrop 15`)
videodrop20 <- na.omit(data$`videodrop 20`)
videodrop26 <- na.omit(data$`videodrop 26`)

df_videodrop <- data.frame(
  concentration = c(videodrop15, videodrop20, videodrop26),
  temperature = factor(rep(c("15°C", "20°C", "26°C"),
                           times = c(length(videodrop15),
                                     length(videodrop20),
                                     length(videodrop26))))
)

# -----------------------
# COMBINE
# -----------------------
df_cyto$method <- "Cytometer"
df_videodrop$method <- "Videodrop"

df_all <- rbind(df_cyto, df_videodrop)

df_all$group <- interaction(df_all$temperature, df_all$method)

df_all$group <- factor(df_all$group, levels = c(
  "15°C.Cytometer", "15°C.Videodrop",
  "20°C.Cytometer", "20°C.Videodrop",
  "26°C.Cytometer", "26°C.Videodrop"
))

# -----------------------
# COLORS
# -----------------------
colors <- rep(c("steelblue", "tomato"), 3)

# -----------------------
# PLOT SETTINGS
# -----------------------
par(
  mar = c(8, 8, 2, 2),
  mgp = c(6, 1.5, 0),
  cex.axis = 1.4,
  cex.lab = 1.8
)

boxplot(concentration ~ group, data = df_all,
        log = "y",
        yaxt = "n",   # on enlève les ticks automatiques Y
        xlab = "",
        ylab = "Concentration (parts/mL)",
        col = colors,
        border = "gray30",
        names = rep("", 6),
        las = 1,
        outline = FALSE)

# axe Y personnalisé (log scale)
axis(2, at = c(2e5, 1e6, 1e7, 1e8, 5e8),
     labels = c("2e5", "1e6", "1e7", "1e8", "5e8"),
     las = 1,
     cex.axis = 1.4)

# séparations verticales
abline(v = c(2.5, 4.5), lty = 2, col = "gray60", lwd = 1.5)

# labels des groupes
axis(1, at = c(1.5, 3.5, 5.5),
     labels = c("15°C", "20°C", "26°C"),
     tick = FALSE,
     cex.axis = 1.5)

legend("bottomleft",
       legend = c("Cytometer", "Videodrop"),
       fill = c("steelblue", "tomato"),
       bty = "n",
       cex = 1.3)