library(readxl)

data <- read_excel(
  "bioinformatics/boxplot_cyto_vs_videodrop/Bilan extraction ADN.xlsx",
  sheet = "Compare cyto videodrop"
)

# 15°C
wilcox_15 <- wilcox.test(
  data$`cyto 15`,
  data$`videodrop 15`,
  paired = FALSE,
  na.action = na.omit
)

# 20°C
wilcox_20 <- wilcox.test(
  data$`cyto 20`,
  data$`videodrop 20`,
  paired = FALSE,
  na.action = na.omit
)

# 26°C
wilcox_26 <- wilcox.test(
  data$`cyto 26`,
  data$`videodrop 26`,
  paired = FALSE,
  na.action = na.omit,
  exact=FALSE
)

# -----------------------
# RESULTS
# -----------------------

wilcox_15
wilcox_20
wilcox_26