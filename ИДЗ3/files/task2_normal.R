# ============================================================
#  Вариант 14  —  Задание 2: Нормальное распределение
#  Параметры: alpha2=0.10, c=2.60, d=4.60, h=0.80,
#             a0=9.00, sigma0=2.00, a1=3.00, sigma1=2.00
# ============================================================

# ---------- Исходные данные ----------
x <- c(8.752, 2.202, 2.779, 2.717, 2.033, 2.755, 3.688, 4.318, 2.188, 2.878,
       2.243, 3.204, 3.076, 4.700, 2.244, 4.634, 3.069, 1.574, 5.394, 4.323,
       2.691, 1.079, 2.726, 1.260, -3.224, 3.050, 3.714, 4.449, 2.096, 3.271,
       3.834, 1.621, 3.466, 3.764, 2.349, -1.848, 4.301, 1.119, 0.639, 7.699,
       4.572, 3.146, 2.898, 3.383, 3.535, 3.054, 6.165, 4.822, 2.797, -0.056)

n      <- length(x)
alpha2 <- 0.10
c_par  <- 2.60
d_par  <- 4.60
h      <- 0.80
a0     <- 9.00
sigma0 <- 2.00
a1     <- 3.00
sigma1 <- 2.00

cat("=== Исходные данные ===\n")
cat("n =", n, "\n")
cat("Наблюдения:\n"); print(x); cat("\n")

# ============================================================
# (a) Вариационный ряд, ECDF, гистограмма, полигон частот (шаг h)
# ============================================================
cat("=== (a) Вариационный ряд и графики ===\n")
var_row <- sort(x)
cat("Вариационный ряд:\n"); print(var_row); cat("\n")

# Границы интервалов с шагом h
x_min  <- floor(min(x) / h) * h
x_max  <- ceiling(max(x) / h) * h
breaks <- seq(x_min, x_max, by = h)

# Таблица интервальных частот
freq_cut <- cut(x, breaks = breaks, include.lowest = TRUE, right = FALSE)
freq_tbl <- table(freq_cut)
rel_freq <- freq_tbl / n
cum_freq <- cumsum(rel_freq)
midpoints <- (breaks[-length(breaks)] + breaks[-1]) / 2

cat("Интервальная таблица частот (шаг h =", h, "):\n")
print(data.frame(
  interval  = names(freq_tbl),
  midpoint  = midpoints,
  n_i       = as.integer(freq_tbl),
  rel_n_i   = round(as.numeric(rel_freq), 4),
  F_emp     = round(as.numeric(cum_freq), 4)
))
cat("\n")

par(mfrow = c(1, 3))

# Гистограмма
hist(x, breaks = breaks, freq = FALSE,
     main = "Гистограмма частот",
     xlab = "x", ylab = "Плотность",
     col = "steelblue", border = "white")

# Полигон частот
plot(midpoints, as.numeric(rel_freq) / h,
     type = "b", pch = 19, col = "darkgreen", lwd = 2,
     main = "Полигон частот",
     xlab = "x", ylab = "Плотность")
grid()

# ECDF
plot(ecdf(x),
     main = "Эмпирическая функция распределения",
     xlab = "x", ylab = "F*(x)",
     col = "darkred", lwd = 2)
grid()

par(mfrow = c(1, 1))

# ============================================================
# (b) Выборочные характеристики
# ============================================================
cat("=== (b) Выборочные характеристики ===\n")

x_mean <- mean(x)
x_var  <- var(x)                      # несмещённая s²
x_var_b <- mean((x - x_mean)^2)       # смещённая D
x_med  <- median(x)

m2 <- mean((x - x_mean)^2)
m3 <- mean((x - x_mean)^3)
m4 <- mean((x - x_mean)^4)

asymm <- m3 / m2^(3/2)
kurt  <- m4 / m2^2 - 3

p_cd <- mean(x >= c_par & x <= d_par)

cat("(i)   Выборочное среднее:       x̄ =", round(x_mean, 4), "\n")
cat("(ii)  Дисперсия (несмещ. s²):   s² =", round(x_var, 4), "\n")
cat("      Дисперсия (смещ. D):        D =", round(x_var_b, 4), "\n")
cat("(iii) Медиана:                   Me =", round(x_med, 4), "\n")
cat("(iv)  Коэффициент асимметрии:    As =", round(asymm, 4), "\n")
cat("(v)   Эксцесс:                   Ex =", round(kurt, 4), "\n")
cat("(vi)  P(X ∈ [", c_par, ",", d_par, "]):  P* =", round(p_cd, 4), "\n\n")

# ============================================================
# (c) Оценки параметров нормального распределения (a, σ²)
#     МПП и метод моментов; смещение оценок
# ============================================================
cat("=== (c) Оценки параметров N(a, σ²) ===\n")

# МПП: â = x̄,  σ̂²_mle = (1/n)Σ(xi - x̄)² (смещённая)
a_hat_mle     <- x_mean
sigma2_hat_mle <- x_var_b    # смещённая

# Метод моментов: те же формулы
a_hat_mm      <- x_mean
sigma2_hat_mm <- x_var_b

# Несмещённая оценка дисперсии
sigma2_unbiased <- x_var     # s² = 1/(n-1)*Σ(xi-x̄)²

cat("МПП:\n")
cat("  â_mle     =", round(a_hat_mle, 4), "  (несмещённая, E[â]=a)\n")
cat("  σ̂²_mle   =", round(sigma2_hat_mle, 4), "  (СМЕЩЁННАЯ: E[σ̂²_mle] = σ²·(n-1)/n)\n")
cat("  Смещение σ̂²_mle: bias = -σ²/n ≈", round(-sigma2_unbiased / n, 4), "\n\n")

cat("Метод моментов:\n")
cat("  â_mm      =", round(a_hat_mm, 4), "  (совпадает с МПП)\n")
cat("  σ̂²_mm    =", round(sigma2_hat_mm, 4), "  (совпадает с МПП, тоже смещённая)\n\n")

cat("Несмещённая оценка дисперсии: s² =", round(sigma2_unbiased, 4), "\n\n")

# ============================================================
# (d) Доверительные интервалы для a и σ²  (уровень 1 - alpha2)
# ============================================================
cat("=== (d) Доверительные интервалы уровня", 1 - alpha2, "===\n")

# ДИ для a: x̄ ± t_{α/2, n-1} * s/√n
t_crit <- qt(1 - alpha2 / 2, df = n - 1)
s      <- sqrt(sigma2_unbiased)
ci_a_low  <- x_mean - t_crit * s / sqrt(n)
ci_a_high <- x_mean + t_crit * s / sqrt(n)

cat("ДИ для a (t-распределение, df =", n-1, "):\n")
cat("  t_{1-α/2} =", round(t_crit, 4), "\n")
cat("  [", round(ci_a_low, 4), ";", round(ci_a_high, 4), "]\n\n")

# ДИ для σ²: χ²-распределение
chi2_lo <- qchisq(alpha2 / 2, df = n - 1)
chi2_hi <- qchisq(1 - alpha2 / 2, df = n - 1)
ci_s2_low  <- (n - 1) * sigma2_unbiased / chi2_hi
ci_s2_high <- (n - 1) * sigma2_unbiased / chi2_lo

cat("ДИ для σ² (χ²-распределение, df =", n-1, "):\n")
cat("  χ²_{α/2}  =", round(chi2_lo, 4), "\n")
cat("  χ²_{1-α/2}=", round(chi2_hi, 4), "\n")
cat("  [", round(ci_s2_low, 4), ";", round(ci_s2_high, 4), "]\n\n")

# ============================================================
# (e) Критерий Колмогорова: H0: N(a0, σ0²)
# ============================================================
cat("=== (e) Критерий Колмогорова: H0: N(a0 =", a0, ", σ0 =", sigma0, ") ===\n")

# D_n = max|F*(x) - F0(x)|
sorted_x <- sort(x)
F_emp_hi <- (1:n) / n          # F*(x_i+)
F_emp_lo <- (0:(n-1)) / n      # F*(x_i-)
F0       <- pnorm(sorted_x, mean = a0, sd = sigma0)

D_plus  <- max(F_emp_hi - F0)
D_minus <- max(F0 - F_emp_lo)
D_n     <- max(D_plus, D_minus)

# Статистика Колмогорова
KS_stat <- sqrt(n) * D_n

# Критическое значение (таблица Колмогорова, уровень alpha2=0.10)
# P(K > k_alpha) = alpha: k_{0.10} = 1.2238
ks_crit <- 1.2238   # стандартное табличное значение для alpha=0.10

# p-value через ks.test
ks_result <- ks.test(x, "pnorm", mean = a0, sd = sigma0)

cat("D_n =", round(D_n, 4), "\n")
cat("√n · D_n =", round(KS_stat, 4), "\n")
cat("Критическое значение (α =", alpha2, "): k_α =", ks_crit, "\n")
cat("p-значение =", round(ks_result$p.value, 4), "\n")
if (KS_stat > ks_crit) {
  cat("Вывод: ОТВЕРГАЕМ H0 на уровне", alpha2, "\n")
} else {
  cat("Вывод: НЕТ оснований отвергнуть H0 на уровне", alpha2, "\n")
}
cat("Наибольший уровень значимости (p-value) =", round(ks_result$p.value, 4), "\n\n")

# ============================================================
# (f) Критерий χ²: простая гипотеза H0: N(a0, σ0²)
# ============================================================
cat("=== (f) χ²-критерий: простая гипотеза H0: N(a0 =", a0, ", σ0² =", sigma0^2, ") ===\n")

# Равновероятные интервалы (каждый ~5 наблюдений => k=10 групп)
k_f <- 10
probs_f <- seq(0, 1, length.out = k_f + 1)
breaks_f <- qnorm(probs_f, mean = a0, sd = sigma0)
breaks_f[1]   <- -Inf
breaks_f[k_f+1] <- Inf

obs_f <- as.integer(table(cut(x, breaks = breaks_f, include.lowest = TRUE)))
exp_f <- rep(n / k_f, k_f)

chi2_f_stat <- sum((obs_f - exp_f)^2 / exp_f)
df_f_stat   <- k_f - 1
chi2_cr_f   <- qchisq(1 - alpha2, df = df_f_stat)
p_val_f     <- 1 - pchisq(chi2_f_stat, df = df_f_stat)

cat("Число групп k =", k_f, ", ожидаемое в каждой =", n/k_f, "\n")
cat(sprintf("%-5s  %-6s  %-8s\n", "Группа", "Набл.", "Ожид."))
for (i in 1:k_f) cat(sprintf("%-5d  %-6d  %-8.3f\n", i, obs_f[i], exp_f[i]))
cat("\nχ² =", round(chi2_f_stat, 4), " df =", df_f_stat,
    " крит. =", round(chi2_cr_f, 4), " p =", round(p_val_f, 4), "\n")
cat("Вывод:", ifelse(chi2_f_stat > chi2_cr_f, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"),
    "на уровне", alpha2, "\n")
cat("Наибольший уровень значимости =", round(p_val_f, 4), "\n\n")

# ============================================================
# (g) Критерий χ²: сложная гипотеза H0: N(a, σ²) (оба параметра оцениваются)
# ============================================================
cat("=== (g) χ²-критерий: сложная гипотеза H0: N(â, σ̂²) ===\n")

a_est     <- x_mean
sigma_est <- sqrt(x_var_b)    # смещённая оценка σ (МПП)

k_g       <- 8    # число интервалов
probs_g   <- seq(0, 1, length.out = k_g + 1)
breaks_g  <- qnorm(probs_g, mean = a_est, sd = sigma_est)
breaks_g[1]     <- -Inf
breaks_g[k_g+1] <- Inf

obs_g <- as.integer(table(cut(x, breaks = breaks_g, include.lowest = TRUE)))
exp_g <- rep(n / k_g, k_g)

chi2_g_stat <- sum((obs_g - exp_g)^2 / exp_g)
df_g_stat   <- k_g - 1 - 2    # -2 оцениваемых параметра
chi2_cr_g   <- qchisq(1 - alpha2, df = df_g_stat)
p_val_g     <- 1 - pchisq(chi2_g_stat, df = df_g_stat)

cat("â =", round(a_est, 4), "  σ̂ =", round(sigma_est, 4), "\n")
cat("Число групп k =", k_g, ", df =", df_g_stat, "\n")
cat(sprintf("%-5s  %-6s  %-8s\n", "Группа", "Набл.", "Ожид."))
for (i in 1:k_g) cat(sprintf("%-5d  %-6d  %-8.3f\n", i, obs_g[i], exp_g[i]))
cat("\nχ² =", round(chi2_g_stat, 4), " df =", df_g_stat,
    " крит. =", round(chi2_cr_g, 4), " p =", round(p_val_g, 4), "\n")
cat("Вывод:", ifelse(chi2_g_stat > chi2_cr_g, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"),
    "на уровне", alpha2, "\n")
cat("Наибольший уровень значимости =", round(p_val_g, 4), "\n\n")

# ============================================================
# (h) НМК Неймана-Пирсона:
#     H0: (a,σ²)=(a0,σ0²)  vs  H1: (a,σ²)=(a1,σ1²)
# ============================================================
cat("=== (h) НМК Неймана-Пирсона ===\n")
cat("H0: N(a0=", a0, ", σ0=", sigma0, ")  vs  H1: N(a1=", a1, ", σ1=", sigma1, ")\n")

# Отношение правдоподобий: Λ = L(H1)/L(H0)
# log Λ = Σ[ log f1(xi) - log f0(xi) ]
log_L0 <- sum(dnorm(x, mean = a0, sd = sigma0, log = TRUE))
log_L1 <- sum(dnorm(x, mean = a1, sd = sigma1, log = TRUE))
log_Lambda <- log_L1 - log_L0

cat("log L(H0) =", round(log_L0, 4), "\n")
cat("log L(H1) =", round(log_L1, 4), "\n")
cat("log Λ = log L(H1) - log L(H0) =", round(log_Lambda, 4), "\n\n")

# Критическая область: log Λ > c  (отвергаем H0 в пользу H1)
# Находим c через распределение log Λ при H0 (bootstrap)
set.seed(42)
B <- 10000
log_Lambda_boot <- replicate(B, {
  x_b    <- rnorm(n, mean = a0, sd = sigma0)
  sum(dnorm(x_b, mean = a1, sd = sigma1, log = TRUE)) -
    sum(dnorm(x_b, mean = a0, sd = sigma0, log = TRUE))
})
c_crit_h <- quantile(log_Lambda_boot, 1 - alpha2)
real_alpha_h <- mean(log_Lambda_boot > c_crit_h)

cat("Критическое значение c (bootstrap, α =", alpha2, ") =", round(c_crit_h, 4), "\n")
cat("Реальный уровень значимости ≈", round(real_alpha_h, 4), "\n")
cat("log Λ_набл =", round(log_Lambda, 4), ":",
    ifelse(log_Lambda > c_crit_h, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# Мощность
log_Lambda_h1 <- replicate(B, {
  x_b <- rnorm(n, mean = a1, sd = sigma1)
  sum(dnorm(x_b, mean = a1, sd = sigma1, log = TRUE)) -
    sum(dnorm(x_b, mean = a0, sd = sigma0, log = TRUE))
})
power_h <- mean(log_Lambda_h1 > c_crit_h)
cat("Мощность критерия ≈", round(power_h, 4), "\n\n")

# Что если поменять H0 и H1 местами?
cat("--- Если поменять H0 и H1 ---\n")
cat("H0: N(a1 =", a1, ", σ1 =", sigma1, ")  vs  H1: N(a0 =", a0, ", σ0 =", sigma0, ")\n")
log_Lambda_swap <- -log_Lambda   # меняем знак
log_Lambda_boot_swap <- replicate(B, {
  x_b <- rnorm(n, mean = a1, sd = sigma1)
  sum(dnorm(x_b, mean = a0, sd = sigma0, log = TRUE)) -
    sum(dnorm(x_b, mean = a1, sd = sigma1, log = TRUE))
})
c_crit_swap  <- quantile(log_Lambda_boot_swap, 1 - alpha2)
real_alpha_swap <- mean(log_Lambda_boot_swap > c_crit_swap)
cat("c' =", round(c_crit_swap, 4), "  α* ≈", round(real_alpha_swap, 4), "\n")
cat("log Λ'_набл =", round(log_Lambda_swap, 4), ":",
    ifelse(log_Lambda_swap > c_crit_swap, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# ============================================================
# (i) Замена на распределение Лапласа f(x) = 1/(√2 σ) exp(-√2/σ |x-a|)
# ============================================================
cat("=== (i) Распределение Лапласа f(x) = 1/(√2·σ)·exp(-√2/σ·|x-a|) ===\n")

dlaplace <- function(x, mu, sigma) {
  (1 / (sqrt(2) * sigma)) * exp(-sqrt(2) / sigma * abs(x - mu))
}
plaplace <- function(x, mu, sigma) {
  ifelse(x < mu,
         0.5 * exp( sqrt(2) / sigma * (x - mu)),
         1 - 0.5 * exp(-sqrt(2) / sigma * (x - mu)))
}

# МПП для Лапласа: â = медиана, σ̂ = √2 * mean|xi - â|
a_lap   <- median(x)
sigma_lap <- sqrt(2) * mean(abs(x - a_lap))

cat("(c) МПП для Лапласа:\n")
cat("  â_lap (медиана) =", round(a_lap, 4), "\n")
cat("  σ̂_lap           =", round(sigma_lap, 4), "\n")
cat("  Смещение â: E[медиана] = a (несмещённая)\n")
cat("  Смещение σ̂: асимптотически несмещённая\n\n")

# ДИ для a Лапласа через CLT: Var(медиана) = σ²/(2n) * (1/f²(a))
# f(a) = 1/(√2 σ), => 1/f²(a) = 2σ², Var(медиана) = σ²/n
se_lap_a <- sigma_lap / sqrt(n)
z_cr <- qnorm(1 - alpha2 / 2)
ci_lap_a_lo <- a_lap - z_cr * se_lap_a
ci_lap_a_hi <- a_lap + z_cr * se_lap_a
cat("(d) ДИ для a_Лапласа уровня", 1 - alpha2, ":\n")
cat("  [", round(ci_lap_a_lo, 4), ";", round(ci_lap_a_hi, 4), "]\n\n")

# Колмогоров: H0: Laplace(a0, sigma0)
cat("(e) Критерий Колмогорова: H0: Laplace(a0=", a0, ", σ0=", sigma0, ")\n")
F0_lap   <- plaplace(sorted_x, a0, sigma0)
D_lap    <- max(pmax(abs(F_emp_hi - F0_lap), abs(F0_lap - F_emp_lo)))
KS_lap   <- sqrt(n) * D_lap
cat("  √n·D_n =", round(KS_lap, 4), "  крит. =", ks_crit, "\n")
cat("  Вывод:", ifelse(KS_lap > ks_crit, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# χ² простая: H0: Laplace(a0, sigma0)
cat("(f) χ²: простая H0: Laplace(a0, σ0)\n")
breaks_lap <- c(-Inf, quantile(x, probs = seq(0.1, 0.9, by = 0.1)), Inf)
obs_lap_s  <- as.integer(table(cut(x, breaks = breaks_lap, include.lowest = TRUE)))
p_lap_s    <- diff(plaplace(breaks_lap[-c(1, length(breaks_lap))], a0, sigma0))
p_lap_s    <- c(plaplace(breaks_lap[2], a0, sigma0),
                p_lap_s,
                1 - plaplace(breaks_lap[length(breaks_lap)-1], a0, sigma0))
exp_lap_s  <- n * p_lap_s
chi2_lap_s <- sum((obs_lap_s - exp_lap_s)^2 / exp_lap_s)
df_lap_s   <- length(obs_lap_s) - 1
cat("  χ² =", round(chi2_lap_s, 4), " df =", df_lap_s,
    " крит. =", round(qchisq(1 - alpha2, df_lap_s), 4),
    " p =", round(1 - pchisq(chi2_lap_s, df_lap_s), 4), "\n")
cat("  Вывод:", ifelse(chi2_lap_s > qchisq(1-alpha2, df_lap_s),
                       "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# χ² сложная: H0: Laplace(â, σ̂)
cat("(g) χ²: сложная H0: Laplace(â, σ̂)\n")
p_lap_c    <- diff(plaplace(breaks_lap[-c(1, length(breaks_lap))], a_lap, sigma_lap))
p_lap_c    <- c(plaplace(breaks_lap[2], a_lap, sigma_lap),
                p_lap_c,
                1 - plaplace(breaks_lap[length(breaks_lap)-1], a_lap, sigma_lap))
exp_lap_c  <- n * p_lap_c
chi2_lap_c <- sum((obs_lap_s - exp_lap_c)^2 / exp_lap_c)
df_lap_c   <- length(obs_lap_s) - 1 - 2
cat("  χ² =", round(chi2_lap_c, 4), " df =", df_lap_c,
    " крит. =", round(qchisq(1 - alpha2, df_lap_c), 4),
    " p =", round(1 - pchisq(chi2_lap_c, df_lap_c), 4), "\n")
cat("  Вывод:", ifelse(chi2_lap_c > qchisq(1-alpha2, df_lap_c),
                       "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# НМК: H0: Lap(a0,σ0) vs H1: Lap(a1,σ1)
cat("(h) НМК: H0: Laplace(a0,σ0) vs H1: Laplace(a1,σ1)\n")
log_L0_lap <- sum(log(dlaplace(x, a0, sigma0)))
log_L1_lap <- sum(log(dlaplace(x, a1, sigma1)))
log_Lam_lap <- log_L1_lap - log_L0_lap
cat("  log Λ =", round(log_Lam_lap, 4), "\n")

set.seed(123)
boot_lap <- replicate(B, {
  xb <- a0 + sigma0 / sqrt(2) * log(2 * runif(n)) * sample(c(-1,1), n, replace=TRUE)
  sum(log(dlaplace(xb, a1, sigma1))) - sum(log(dlaplace(xb, a0, sigma0)))
})
c_lap <- quantile(boot_lap, 1 - alpha2)
cat("  Критическое значение c_lap =", round(c_lap, 4), "\n")
cat("  Вывод:", ifelse(log_Lam_lap > c_lap, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

cat("=== РЕШЕНИЕ ЗАВЕРШЕНО ===\n")
