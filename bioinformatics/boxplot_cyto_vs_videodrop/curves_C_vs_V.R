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

col_C <- "#2166ac"
col_V <- "#d6604d"

temps <- c(15, 20, 26)
traitements <- c("coevo", "virus")

par(mfrow = c(2, 3), mar = c(4, 5, 3, 1), oma = c(0, 0, 2, 0))

for (trt in traitements) {
  for (temp in temps) {
    
    key_C <- paste0(trt, "_", temp, "_C")
    key_V <- paste0(trt, "_", temp, "_V")
    
    yC <- data[[key_C]]
    yV <- data[[key_V]]
    
    # Supprimer valeurs <= 0 (log impossible)
    yC <- yC[yC > 0]
    yV <- yV[yV > 0]
    
    n_gen <- max(length(yC), length(yV))
    gen_C <- seq_along(yC)
    gen_V <- seq_along(yV)
    
    ymin <- min(c(yC, yV), na.rm = TRUE)
    ymax <- max(c(yC, yV), na.rm = TRUE)
    
    # Définir bornes log10 propres
    log_min <- floor(log10(ymin))
    log_max <- ceiling(log10(ymax))
    y_ticks <- 10^(log_min:log_max)
    
    plot(gen_C, yC,
         type = "b",
         pch = 16,
         col = col_C,
         log = "y",
         xlim = c(1, n_gen),
         ylim = c(10^log_min, 10^log_max),
         xlab = "Generations",
         ylab = "Concentration (parts/mL)",
         main = paste0(trt, " — ", temp, "°C"),
         axes = FALSE,
         frame.plot = FALSE)
    
    axis(1, at = 1:n_gen, labels = 1:n_gen)
    axis(2, at = y_ticks,
         labels = parse(text = paste0("10^", log_min:log_max)),
         las = 1)
    
    box()
    
    lines(gen_V, yV,
          type = "b",
          pch = 17,
          col = col_V)
    
    legend("topleft",
           legend = c("Cytometer", "Videodrop"),
           col = c(col_C, col_V),
           pch = c(16, 17),
           lty = 1,
           bty = "n",
           cex = 0.85)
  }
}