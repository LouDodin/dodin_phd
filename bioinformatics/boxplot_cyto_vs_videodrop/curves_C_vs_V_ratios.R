# source("curves.R")

library(readr)

df <- read_csv("bioinformatics/boxplot_cyto_vs_videodrop/Bilan extraction ADN - curves.csv")

col_names <- list(
  coevo_15_C = "coevo 15 C", coevo_15_V = "coevo 15 V",
  coevo_20_C = "coevo 20 C", coevo_20_V = "coevo 20 V",
  coevo_26_C = "coevo 26 C", coevo_26_V = "coevo 26 V",
  virus_15_C = "virus 15 C", virus_15_V = "virus 15 V",
  virus_20_C = "virus 20 C", virus_20_V = "virus 20 V",
  virus_26_C = "virus 26 C", virus_26_V = "virus 26 V"
)

data <- lapply(col_names, function(cn) {
  vals <- as.numeric(df[[cn]])
  vals[!is.na(vals)]
})

col_ratio <- "#542788"

temps <- c(15, 20, 26)
traitements <- c("coevo", "virus")

par(mfrow = c(2, 3), mar = c(4, 5, 3, 1), oma = c(0, 0, 2, 0))

for (trt in traitements) {
  for (temp in temps) {
    
    key_C <- paste0(trt, "_", temp, "_C")
    key_V <- paste0(trt, "_", temp, "_V")
    
    yC <- data[[key_C]]
    yV <- data[[key_V]]
    
    # Garder uniquement positions valides
    valid <- which(!is.na(yC) & !is.na(yV) & yV > 0)
    yC <- yC[valid]
    yV <- yV[valid]
    
    ratio <- yV/(yV+yC)
    gen <- seq_along(ratio)
    
    ymin <- min(ratio, na.rm = TRUE)
    ymax <- max(ratio, na.rm = TRUE)
    
    plot(gen, ratio,
         type = "b",
         pch = 16,
         col = col_ratio,
         ylim = c(0, ymax * 1.1),
         xlab = "Generations",
         ylab = "V/(V+C)",
         main = paste0(trt, " — ", temp, "°C"),
         axes = TRUE,
         frame.plot = FALSE)
    
    abline(h = 0.5, lty = 2)  # équilibre C = V
  }
}

mtext("V/(V+C)",
      outer = TRUE,
      cex = 1.2,
      font = 2)