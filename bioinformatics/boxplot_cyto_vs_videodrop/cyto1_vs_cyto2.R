# -----------------------------
# Create data
# -----------------------------
cyto1 <- c(2.23E+08, 3.41E+08, 4.29E+08, 8.00E+08, 1.27E+09,
           3.37E+08, 7.94E+05, 4.96E+06, 1.14E+06, 3.89E+05, 4.93E+05)

cyto2 <- c(1.18E+08, 1.79E+08, 1.86E+08, 3.56E+08, 5.84E+08,
           1.76E+08, 5.52E+05, 2.89E+05, 2.84E+05, 2.06E+05, 3.40E+05)

df <- data.frame(cyto1, cyto2)

# -----------------------------
# Linear regression (forced through origin)
# -----------------------------
model <- lm(cyto1 ~ 0 + cyto2, data = df)

a <- coef(model)[1]
r2 <- summary(model)$r.squared

# -----------------------------
# Plot
# -----------------------------
plot(cyto2, cyto1,
     pch = 3,            # croix
     cex = 1.2,          # taille des croix
     lwd = 2, 
     col = "darkblue",
     xlab = "Cyto 2",
     ylab = "Cyto 1",
     main = "Linear regression: y = a*x")

# Create regression line manually
x_seq <- seq(min(cyto2), max(cyto2), length.out = 200)
y_seq <- a * x_seq

lines(x_seq, y_seq, col = "red", lwd = 2)

# -----------------------------
# Add equation and R2
# -----------------------------
eq <- paste0("y = ",
             format(a, scientific = TRUE, digits = 4),
             " x")

r2_text <- paste0("R² = ", round(r2, 4))

legend("topleft",
       legend = c(eq, r2_text),
       bty = "n")