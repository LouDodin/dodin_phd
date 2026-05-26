library(readxl)

# Importer les données depuis Excel
data <- read_excel(
  "bioinformatics/boxplot_cyto_vs_videodrop/Bilan extraction ADN.xlsx",
  sheet = "Boxplot ratio Videodrop"
)

# Extraire chaque groupe en supprimant les NA
data15 <- na.omit(data$data15)
data20 <- na.omit(data$data20)
data26 <- na.omit(data$data26)

# Combiner dans un data frame
df <- data.frame(
  concentration = c(data15, data20, data26),
  temperature = factor(
    c(rep("15°C", length(data15)),
      rep("20°C", length(data20)),
      rep("26°C", length(data26))),
    levels = c("15°C", "20°C", "26°C")
  )
)

# Définir les couleurs
colors <- c("green", "yellow", "lightblue")

# 🔹 Ajuster les marges (augmente la marge gauche)
par(mar = c(5, 7, 4, 2))  # bas, gauche, haut, droite

# Créer le boxplot (sans ylab automatique)
boxplot(concentration ~ temperature,
        data = df,
        log = "y",
        main = "Ratio Videodrop 2/Videodrop 1",
        xlab = "",
        ylab = "",
        col = colors,
        border = "gray30",
        las = 1)

# 🔹 Ajouter le titre de l'axe vertical avec un décalage contrôlé
mtext("Ratio Videodrop 2/Videodrop 1",
      side = 2,
      line = 5)  # augmenter si besoin

# Tests de Wilcoxon entre groupes
wilcox_15_20 <- wilcox.test(data15, data20, exact = FALSE)
wilcox_15_26 <- wilcox.test(data15, data26, exact = FALSE)
wilcox_20_26 <- wilcox.test(data20, data26, exact = FALSE)

# Afficher les résultats
print(wilcox_15_20)
print(wilcox_15_26)
print(wilcox_20_26)