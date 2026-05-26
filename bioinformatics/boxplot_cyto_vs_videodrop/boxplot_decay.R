library(readxl)

# Import data from Excel
data <- read_excel(
  "bioinformatics/boxplot_cyto_vs_videodrop/Bilan extraction ADN.xlsx",
  sheet = "Boxplot decay 20 vs 26"
)

# Extract each group
data20 <- na.omit(data$`data20`)
data26 <- na.omit(data$`data26`)

# Combine into a data frame
df <- data.frame(
  concentration = c(data20, data26),
  temperature = factor(
    c(rep("20°C", length(data20)),
      rep("26°C", length(data26)))
  )
)

# Define colors
colors <- c("yellow", "lightblue")

# Create boxplot
boxplot(concentration ~ temperature,
        data = df,
        main = "Decay",
        xlab = "",
        ylab = "Decay rates (day-1)",
        col = colors,
        border = "gray30",
        las = 1)


# wilcox.test(size20, size26, exact = FALSE)