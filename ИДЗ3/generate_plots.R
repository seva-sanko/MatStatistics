# ============================================================
#  Generate all plots for the report
#  Output: PNG files in ./plots/
# ============================================================

dir.create("./plots", showWarnings = FALSE)

# ---- Task 1 data ----
x1 <- c(1,0,1,1,3,1,0,1,0,1,1,1,0,1,1,0,0,1,2,0,
         1,1,0,1,1,0,2,0,0,0,0,1,1,0,0,2,0,1,0,1,
         1,0,1,0,0,0,1,2,0,2)
n1 <- length(x1)
lambda0_1 <- 1.40; lambda_hat_1 <- mean(x1)  # 0.70

# ---- Task 2 data ----
x2 <- c(8.752,2.202,2.779,2.717,2.033,2.755,3.688,4.318,2.188,2.878,
         2.243,3.204,3.076,4.700,2.244,4.634,3.069,1.574,5.394,4.323,
         2.691,1.079,2.726,1.260,-3.224,3.050,3.714,4.449,2.096,3.271,
         3.834,1.621,3.466,3.764,2.349,-1.848,4.301,1.119,0.639,7.699,
         4.572,3.146,2.898,3.383,3.535,3.054,6.165,4.822,2.797,-0.056)
n2 <- length(x2)
h  <- 0.80

# ============================================================
# Plot 1: Task1 — Histogram of frequencies
# ============================================================
png("./plots/task1_histogram.png", width=700, height=500, res=100)
freq <- table(x1)
barplot(freq/n1,
        names.arg = names(freq),
        main = "Histogram of relative frequencies (Task 1)",
        xlab = "x", ylab = "Relative frequency",
        col = "steelblue", border = "white",
        ylim = c(0, 0.55))
grid(nx=NA, ny=NULL, lty=2, col="grey80")
dev.off()

# ============================================================
# Plot 2: Task1 — ECDF
# ============================================================
png("./plots/task1_ecdf.png", width=700, height=500, res=100)
plot(ecdf(x1),
     main = "Empirical CDF (Task 1)",
     xlab = "x", ylab = expression(F^"*"(x)),
     col = "darkred", lwd = 2.5,
     verticals = TRUE, do.points = FALSE)
grid()
dev.off()

# ============================================================
# Plot 3: Task1 — ECDF vs theoretical Poisson(lambda0)
# ============================================================
png("./plots/task1_ecdf_vs_theor.png", width=700, height=500, res=100)
xs <- 0:5
plot(ecdf(x1),
     main = expression(paste("ECDF vs Poisson(", lambda[0], "=1.4) and Poisson(", hat(lambda), "=0.7)")),
     xlab = "x", ylab = "F(x)",
     col = "black", lwd = 2,
     verticals = TRUE, do.points = FALSE,
     xlim = c(-0.5, 5))
# Theoretical CDF Poisson(lambda0)
lines(xs - 0.01, ppois(xs, lambda0_1), type="s", col="red", lwd=2, lty=2)
# Theoretical CDF Poisson(lambda_hat)
lines(xs + 0.01, ppois(xs, lambda_hat_1), type="s", col="blue", lwd=2, lty=3)
legend("bottomright",
       legend = c(expression(F^"*"(x)),
                  expression(paste("Poisson(", lambda[0], "=1.4)")),
                  expression(paste("Poisson(", hat(lambda), "=0.7)"))),
       col = c("black", "red", "blue"),
       lwd = 2, lty = c(1, 2, 3))
grid()
dev.off()

# ============================================================
# Plot 4: Task1 — Chi-squared residuals bar plot (simple H0)
# ============================================================
png("./plots/task1_chi2_simple.png", width=700, height=450, res=100)
obs  <- c(22, 22, 5, 1)
exp0 <- n1 * c(dpois(0:2, lambda0_1), 1 - ppois(2, lambda0_1))
grp  <- c("X=0", "X=1", "X=2", "X>=3")
barplot(rbind(obs, exp0), beside=TRUE,
        names.arg = grp,
        col = c("steelblue", "tomato"),
        main = expression(paste("Observed vs Expected: simple H0 Poisson(", lambda[0], "=1.4)")),
        xlab = "Group", ylab = "Count",
        legend.text = c("Observed", "Expected"),
        args.legend = list(x="topright"))
grid(nx=NA, ny=NULL, lty=2, col="grey80")
dev.off()

# ============================================================
# Plot 5: Task1 — Chi-squared residuals bar plot (complex H0)
# ============================================================
png("./plots/task1_chi2_complex.png", width=700, height=450, res=100)
expc <- n1 * c(dpois(0:2, lambda_hat_1), 1 - ppois(2, lambda_hat_1))
barplot(rbind(obs, expc), beside=TRUE,
        names.arg = grp,
        col = c("steelblue", "seagreen3"),
        main = expression(paste("Observed vs Expected: complex H0 Poisson(", hat(lambda), "=0.7)")),
        xlab = "Group", ylab = "Count",
        legend.text = c("Observed", "Expected"),
        args.legend = list(x="topright"))
grid(nx=NA, ny=NULL, lty=2, col="grey80")
dev.off()

# ============================================================
# Plot 6: Task2 — Histogram + frequency polygon
# ============================================================
png("./plots/task2_histogram.png", width=800, height=500, res=100)
xmin <- floor(min(x2) / h) * h
xmax <- ceiling(max(x2) / h) * h
brks <- seq(xmin, xmax, by=h)
h_obj <- hist(x2, breaks=brks, freq=FALSE,
              main = "Histogram and frequency polygon (Task 2)",
              xlab = "x", ylab = "Density",
              col = "steelblue", border = "white")
# Frequency polygon
lines(h_obj$mids, h_obj$density, type="o",
      col="darkred", lwd=2, pch=19, cex=0.7)
# Normal density overlay (MLE estimates)
curve(dnorm(x, mean=mean(x2), sd=sqrt(mean((x2-mean(x2))^2))),
      add=TRUE, col="orange", lwd=2, lty=2)
legend("topright",
       legend=c("Frequency polygon", expression(paste("N(", hat(a), ",", hat(sigma)^2, ")"))),
       col=c("darkred","orange"), lwd=2, lty=c(1,2), pch=c(19,NA))
dev.off()

# ============================================================
# Plot 7: Task2 — ECDF
# ============================================================
png("./plots/task2_ecdf.png", width=700, height=500, res=100)
plot(ecdf(x2),
     main = "Empirical CDF (Task 2)",
     xlab = "x", ylab = expression(F^"*"(x)),
     col = "darkred", lwd = 2, verticals=TRUE, do.points=FALSE)
grid()
dev.off()

# ============================================================
# Plot 8: Task2 — ECDF vs N(a0, sigma0^2) and N(ahat, sigmahat^2)
# ============================================================
png("./plots/task2_ecdf_vs_theor.png", width=800, height=500, res=100)
plot(ecdf(x2),
     main = "ECDF vs theoretical N distributions",
     xlab = "x", ylab = "F(x)",
     col = "black", lwd = 2, verticals=TRUE, do.points=FALSE,
     xlim = range(x2) + c(-1,1))
# N(a0, sigma0^2) — simple H0
curve(pnorm(x, mean=9, sd=2), add=TRUE, col="red", lwd=2, lty=2)
# N(ahat, sigmahat^2) — complex H0
curve(pnorm(x, mean=mean(x2), sd=sqrt(mean((x2-mean(x2))^2))),
      add=TRUE, col="blue", lwd=2, lty=3)
legend("topleft",
       legend=c(expression(F^"*"(x)),
                expression(paste("N(", a[0], "=9, ", sigma[0], "=2)")),
                expression(paste("N(", hat(a), ",", hat(sigma), ")"))),
       col=c("black","red","blue"), lwd=2, lty=c(1,2,3))
grid()
dev.off()

# ============================================================
# Plot 9: Task2 — QQ-plot vs Normal
# ============================================================
png("./plots/task2_qqplot.png", width=600, height=600, res=100)
qqnorm(x2, main="Q-Q plot vs Normal distribution",
       col="steelblue", pch=19, cex=0.7)
qqline(x2, col="red", lwd=2)
grid()
dev.off()

# ============================================================
# Plot 10: Task2 — Chi-squared residuals (complex H0 Normal)
# ============================================================
png("./plots/task2_chi2_complex.png", width=750, height=450, res=100)
a_hat <- mean(x2); s_hat <- sqrt(mean((x2-a_hat)^2))
k <- 8
probs <- seq(0, 1, length.out=k+1)
bk <- qnorm(probs, mean=a_hat, sd=s_hat)
bk[1] <- -Inf; bk[k+1] <- Inf
obs2 <- as.integer(table(cut(x2, breaks=bk, include.lowest=TRUE)))
exp2 <- rep(n2/k, k)
grp2 <- paste0("G", 1:k)
barplot(rbind(obs2, exp2), beside=TRUE,
        names.arg=grp2,
        col=c("steelblue","seagreen3"),
        main=expression(paste("Observed vs Expected: complex H0 N(", hat(a), ",", hat(sigma)^2, ")")),
        xlab="Group (equal probability)", ylab="Count",
        legend.text=c("Observed","Expected"),
        args.legend=list(x="topright"))
abline(h=n2/k, col="red", lty=2, lwd=1.5)
dev.off()

# ============================================================
# Plot 11: Task1 — NMP power function
# ============================================================
png("./plots/task1_nmk_power.png", width=700, height=450, res=100)
c_crit <- 57
lambdas <- seq(0.3, 2.5, by=0.05)
power <- ppois(c_crit, lambda=n1*lambdas)
plot(lambdas, power, type="l", col="steelblue", lwd=2.5,
     main="Power function of NMP test (Task 1)",
     xlab=expression(lambda), ylab=expression(beta(lambda)))
abline(v=lambda0_1, col="red", lty=2, lwd=1.5)
abline(v=lambda_hat_1, col="blue", lty=2, lwd=1.5)
abline(h=0.05, col="grey50", lty=3)
legend("bottomright",
       legend=c("Power",
                expression(paste(lambda[0],"=1.4")),
                expression(paste(lambda[1],"=0.7")),
                expression(paste(alpha,"=0.05"))),
       col=c("steelblue","red","blue","grey50"),
       lwd=c(2,1.5,1.5,1), lty=c(1,2,2,3))
grid()
dev.off()

cat("All", 11, "plots saved to ./plots/\n")
list.files("./plots/")

