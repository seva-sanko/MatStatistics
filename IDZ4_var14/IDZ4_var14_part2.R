# IDZ 4, Variant 14 (513201226)
# Part 2: Two-way ANOVA

alpha1 <- 0.01
h2    <- 1.40

# Table 2 data from the assignment image (Var 14)
# No:  1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17
Y <- c(21.01,17.88,19.87,20.91,18.23,13.81,16.51,16.99,19.62,19.88,19.16,21.05,13.61,15.73,15.72,13.87,21.18,
# No: 18    19    20    21    22    23    24    25    26    27    28    29    30    31    32    33    34
       20.55,23.00,19.47,16.13,14.77,16.34,17.25,18.14,16.48,19.56,18.17,14.08,13.35,15.11,13.93,20.14,17.55,
# No: 35    36    37    38    39    40
       17.25,18.43,13.69, 7.03,13.99,13.02)

A <- c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,
       3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,
       5,5,5,5,5,5)

B <- c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4,1,
       1,2,2,3,3,4,4,1,1,2,2,3,3,4,4,1,1,
       2,2,3,3,4,4)

n <- length(Y)
cat("n =", n, "| Check A:", length(A), "| Check B:", length(B), "\n")

dat <- data.frame(Y = Y, A = A, B = B)
dat$Af <- as.factor(dat$A)
dat$Bf <- as.factor(dat$B)

d1 <- length(unique(A))
d2 <- length(unique(B))

cat("====================================================\n")
cat("IDZ 4, Variant 14 | PART 2: Two-way ANOVA\n")
cat("====================================================\n\n")
cat("n =", n, "| alpha1 =", alpha1, "| h =", h2, "\n")
cat("Levels A:", sort(unique(A)), "| Levels B:", sort(unique(B)), "\n\n")

# ===== ZADANIE A =====
cat("---------- ZADANIE A: Model + MNK estimates ----------\n")

q <- lm(Y ~ Af * Bf, data = dat)
qs <- summary(q)
qc <- qs$coefficients

cat("Full model Y ~ A * B\n")
print(qc[, c(1,2,3,4)])
cat("\n")

dat$ow <- paste0(dat$A, ":", dat$B)
ow.levels <- levels(as.factor(dat$ow))
n.levels  <- length(ow.levels)

XM <- matrix(0, nrow = n.levels, ncol = n)
for (i in 1:n.levels) {
  XM[i, (dat$ow == ow.levels[i])] <- 1
}
YM <- as.matrix(dat$Y)
S  <- XM %*% t(XM)
S1 <- solve(S)
bhat <- S1 %*% XM %*% YM

SSE <- sum((YM - t(XM) %*% bhat)^2)
s2  <- SSE / (n - d1 * d2)
cat("Unbiased variance s^2 =", round(s2, 6), "\n")
cat("sigma^2 from lm =", round(qs$sigma^2, 6), "\n\n")

# ===== ZADANIE B =====
cat("---------- ZADANIE B: Interaction plot ----------\n")

cell_means <- tapply(dat$Y, list(dat$A, dat$B), mean)
cat("Cell means:\n")
print(round(cell_means, 4))
cat("\n")

A_levels <- sort(unique(A))
B_levels <- sort(unique(B))
colors   <- c("red", "blue", "green", "orange", "purple")

png("/home/claude/plot2_B_interaction.png", width = 800, height = 600)
ylim_r <- range(cell_means, na.rm = TRUE)
plot(A_levels, cell_means[, 1], type = "l", col = colors[1], lwd = 2,
     ylim = ylim_r,
     xlab = "A factor levels", ylab = "E(Y)",
     main = "Interaction plot: A levels at fixed B")
for (j in 2:d2) {
  lines(A_levels, cell_means[, j], col = colors[j], lwd = 2)
  points(A_levels, cell_means[, j], col = colors[j], pch = 19, cex = 1)
}
points(A_levels, cell_means[, 1], col = colors[1], pch = 19, cex = 1)
legend("topright", legend = paste0("B = ", B_levels),
       col = colors[1:d2], lwd = 2, pch = 19)
dev.off()
cat("Interaction plot saved.\n\n")

cat("ANOVA (full model):\n")
print(anova(q))

# ===== ZADANIE C =====
cat("\n---------- ZADANIE C: Residuals normality ----------\n")

res2 <- q$residuals
cat("Residuals summary:\n")
print(summary(res2))
cat("\n")

min_v  <- min(res2)
max_v  <- max(res2)
brks2  <- seq(from = min_v, to = ceiling(max_v / h2) * h2, by = h2)

png("/home/claude/plot2_C_hist.png", width = 800, height = 600)
hh2 <- hist(res2, breaks = brks2, probability = TRUE,
            main = "Residual histogram (two-factor model)",
            xlab = "Values", col = "lightgray", border = "white")
curve(dnorm(x, mean = mean(res2), sd = sd(res2)),
      add = TRUE, col = "red", lwd = 2)
dev.off()

cat("Histogram breaks:", round(hh2$breaks, 4), "\n")
cat("counts:", hh2$counts, "\n\n")

chi2_f2 <- function(sigma) {
  ep  <- diff(pnorm(hh2$breaks, mean = 0, sd = sigma))
  ex  <- n * ep
  idx <- ex > 0
  sum((hh2$counts[idx] - ex[idx])^2 / ex[idx])
}
opt2    <- optimize(chi2_f2, interval = c(0.1, 20))
chi2_s2 <- opt2$objective
k2      <- length(hh2$counts)
df2_chi <- k2 - 1 - 1
crit2   <- qchisq(1 - alpha1, df2_chi)
pval2   <- pchisq(chi2_s2, df2_chi, lower.tail = FALSE)

cat("Chi-square normality test:\n")
cat("  X^2 =", round(chi2_s2, 6), "\n")
cat("  df =", df2_chi, "\n")
cat("  Critical chi^2_{", df2_chi, ",", 1-alpha1, "} =", round(crit2, 6), "\n")
cat("  P-value:", round(pval2, 6), "\n")
cat("  Result:", ifelse(chi2_s2 > crit2, "REJECT H0 (not normal)", "Do NOT reject H0 (normal)"), "\n\n")

# ===== ZADANIE D =====
cat("---------- ZADANIE D: ANOVA table ----------\n")

q12 <- lm(Y ~ Af + Bf, data = dat)
q1  <- lm(Y ~ Bf, data = dat)
q2  <- lm(Y ~ Af, data = dat)
q0  <- lm(Y ~ 1,  data = dat)

rdf <- n - d1 * d2

make_row <- function(label, mod_red, mod_full) {
  a    <- anova(mod_red, mod_full)
  df_h <- a$Df[2]
  ssh  <- a$`Sum of Sq`[2]
  f_v  <- a$F[2]
  pv   <- a$`Pr(>F)`[2]
  xal  <- qf(1 - alpha1, df_h, rdf)
  mrss <- ssh / df_h
  data.frame(H = label, Res.Df = rdf, SS = round(ssh,4), Df = df_h,
             MS = round(mrss,4), F.stat = round(f_v,4),
             F.crit = round(xal,4), p.value = round(pv,6))
}

err_rss <- deviance(q)
AOV <- rbind(
  make_row("H12 (no A:B)", q12, q),
  make_row("H1  (no A)",   q1,  q),
  make_row("H2  (no B)",   q2,  q),
  make_row("H0  (const)",  q0,  q)
)
AOV <- rbind(AOV, data.frame(H="Err", Res.Df=rdf, SS=round(err_rss,4),
                              Df=NA, MS=round(err_rss/rdf,4),
                              F.stat=NA, F.crit=NA, p.value=NA))
cat("ANOVA table:\n")
print(AOV)
cat("\n")

# ===== ZADANIE E =====
cat("---------- ZADANIE E: AIC and BIC ----------\n")

model_names <- c("Full (A*B)", "Additive (A+B)", "Only B", "Only A", "Constant")
models_list <- list(q, q12, q1, q2, q0)

aic_vals <- sapply(models_list, AIC)
bic_vals <- sapply(models_list, BIC)

mc2 <- data.frame(Model = model_names,
                  AIC   = round(aic_vals, 4),
                  BIC   = round(bic_vals, 4))
print(mc2)
cat("\nBest by AIC:", model_names[which.min(aic_vals)], "\n")
cat("Best by BIC:", model_names[which.min(bic_vals)], "\n\n")

# ===== ZADANIE F =====
cat("---------- ZADANIE F: Summary ----------\n")
cat("1. s^2 =", round(s2, 4), "\n")
cat("2. Interaction p-value:", round(anova(q)[3,5], 6), "\n")
cat("3. Normality:", ifelse(chi2_s2 <= crit2, "confirmed", "rejected"), "\n")
cat("4. Best model AIC:", model_names[which.min(aic_vals)], "\n")
cat("5. Best model BIC:", model_names[which.min(bic_vals)], "\n")
