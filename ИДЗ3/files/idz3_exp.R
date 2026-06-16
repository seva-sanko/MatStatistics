# ============================================================
#  ИДЗ-3, Часть 2: Показательное распределение Exp(λ)
#  Параметры: alpha2=0.20, c=2.40, d=4.80, h=1.20,
#             lambda0=0.17, lambda1=0.33
# ============================================================

alpha2  <- 0.20
c_par   <- 2.40
d_par   <- 4.80
h       <- 1.20
lambda0 <- 0.17   # sigma0
lambda1 <- 0.33   # sigma1

x <- c(0.66, 2.88, 2.17, 9.70, 3.05, 0.64, 0.68, 0.96, 2.36, 0.49,
       6.69, 8.83, 6.37, 0.17, 1.21, 4.73, 1.52, 1.06, 7.64, 0.23,
       7.21, 4.37, 1.89, 0.13, 0.41, 1.36, 1.80, 0.50, 3.34, 1.34,
       1.04, 0.40, 2.48, 0.52, 4.57, 9.11, 0.99, 7.19, 3.15, 2.63,
       0.72, 0.19, 13.33, 4.69, 0.26, 5.01, 8.70, 2.07, 3.27, 0.05)

n  <- length(x)
xs <- sort(x)

cat("=== Исходные данные ===\n")
cat("n =", n, "\n")
cat("Данные:", x, "\n\n")

# ============================================================
# (1) Вариационный ряд, ECDF, гистограмма и полигон частот
# ============================================================
cat("=== (1) Вариационный ряд ===\n")
cat(xs, "\n\n")

# Границы интервалов с шагом h
m1   <- floor(min(x) / h) * h
brks <- seq(m1, max(x) + h, by = h)

hist_data  <- hist(x, breaks = brks, plot = FALSE)
midpoints  <- hist_data$mids
frequencies <- hist_data$counts

par(mfrow = c(1, 2))

# Гистограмма + полигон
hist(x, breaks = brks, freq = TRUE,
     main = "Гистограмма и полигон частот",
     xlab = "Значения", ylab = "Частота",
     col = "steelblue", border = "white")
lines(midpoints, frequencies, type = "o", col = "red", lwd = 2, pch = 19)

# ECDF
plot.ecdf(x, pch = NA, verticals = TRUE,
          main = "Эмпирическая функция распределения",
          ylab = "F*(x)", xlab = "x", lwd = 1.5)
grid()

par(mfrow = c(1, 1))

# ============================================================
# (2) Выборочные характеристики
# ============================================================
cat("=== (2) Выборочные характеристики ===\n")

m    <- mean(x)
# Смещённая дисперсия (момент)
s2   <- mean((x - m)^2)
# Медиана
med  <- median(xs)
# Асимметрия и эксцесс
mu3  <- mean((x - m)^3)
mu4  <- mean((x - m)^4)
asi  <- mu3 / s2^(3/2)
exc  <- mu4 / s2^2 - 3
# Вероятность
prb  <- mean(x >= c_par & x <= d_par)

cat("(i)   Выборочное среднее:     x̄ =", round(m, 4), "\n")
cat("(ii)  Дисперсия (смещ.):       D =", round(s2, 4), "\n")
cat("      Дисперсия (несмещ.):     s² =", round(var(x), 4), "\n")
cat("(iii) Медиана:                 Me =", round(med, 4), "\n")
cat("(iv)  Коэф. асимметрии:        As =", round(asi, 4), "\n")
cat("(v)   Эксцесс:                 Ex =", round(exc, 4), "\n")
cat("(vi)  P(X ∈ [", c_par, ";", d_par, "]):  P* =", round(prb, 4), "\n\n")

# ============================================================
# (3) Оценки параметра λ показательного распределения
#     Exp(λ): f(x) = λ·exp(-λx), E[X] = 1/λ
#     МПП: λ̂ = n / Σxᵢ = 1/x̄
#     ММ:  λ̂ = 1/x̄  (совпадает)
#     Смещение МПП: E[λ̂_МПП] = n/(n-1) · λ  => bias = λ/(n-1)
# ============================================================
cat("=== (3) Оценки параметра λ ===\n")

lambda_mle <- 1 / m          # МПП = ММ для показательного
lambda_mm  <- 1 / m
bias       <- lambda_mle / (n - 1)  # смещение λ̂_МПП

cat("МПП:            λ̂_МПП =", round(lambda_mle, 6), "\n")
cat("Метод моментов: λ̂_ММ  =", round(lambda_mm,  6), "\n")
cat("Смещение МПП:   bias  =  λ̂/(n-1) =", round(bias, 6),
    "(положительное — оценка завышает λ)\n\n")

# ============================================================
# (4) Асимптотический доверительный интервал для λ (уровень 1-α₂)
#     Var(λ̂) ≈ λ²/n  =>  SE = λ̂/√n
# ============================================================
cat("=== (4) Доверительный интервал для λ (уровень", 1-alpha2, ") ===\n")

z    <- qnorm(1 - alpha2 / 2)
se   <- lambda_mle / sqrt(n)
ci_lo <- lambda_mle - z * se
ci_hi <- lambda_mle + z * se

cat("z_{1-α₂/2} =", round(z, 6), "\n")
cat("SE(λ̂) =", round(se, 6), "\n")
cat("ДИ: [", round(ci_lo, 6), ";", round(ci_hi, 6), "]\n\n")

# ============================================================
# (5) Критерий Колмогорова: H0: Exp(λ0)  (простая гипотеза)
# ============================================================
cat("=== (5) Критерий Колмогорова: H0: Exp(λ0 =", lambda0, ") ===\n")

F_emp_hi <- (1:n) / n
F_emp_lo <- (0:(n-1)) / n
F0       <- pexp(xs, rate = lambda0)

ldif   <- abs(F0 - F_emp_lo)
rdif   <- abs(F0 - F_emp_hi)
maxdif <- pmax(ldif, rdif)

j <- which.max(maxdif)
cat("Точка максимального расхождения:\n")
KS.tab <- data.frame(i = 1:n, xi = xs, lecdf = F_emp_lo, recdf = F_emp_hi,
                     hcdf = round(F0,6), ldif = round(ldif,6),
                     rdif = round(rdif,6), maxdif = round(maxdif,6))
print(KS.tab[j, ])

Dn    <- max(maxdif)
D_stat <- Dn * sqrt(n)

# Критическое значение: D_α = C(α)/√n, для α=0.20: C(0.20)≈0.8
D_alpha <- 0.8 / sqrt(n)   # ≈ 0.1131

cat("\nDn =", round(Dn, 6), "\n")
cat("√n · Dn =", round(D_stat, 6), "\n")
cat("Критическое D_α ≈ 0.8/√n =", round(D_alpha, 6), "\n")

# p-value через ks.test
ks_res <- ks.test(xs, "pexp", rate = lambda0)
cat("p-value =", ks_res$p.value, "\n")
if (Dn > D_alpha) {
  cat("Вывод: ОТВЕРГАЕМ H0 (Dn > D_α)\n")
} else {
  cat("Вывод: НЕТ оснований отвергнуть H0 (Dn <= D_α)\n")
}
cat("(Наибольший уровень, при котором принимаем: p =",
    round(ks_res$p.value, 6), ")\n\n")

# График
plot(ecdf(x), main = "Эмпирическая и теоретическая ФР (λ₀ = 0.17)",
     xlab = "x", ylab = "F(x)", verticals = TRUE, col = "black", lwd = 1.5)
curve(pexp(x, rate = lambda0), add = TRUE, col = "red", lwd = 2)
legend("bottomright", legend = c("Эмпирическая F*(x)", "Теор. Exp(λ₀)"),
       col = c("black", "red"), lwd = 2)

# ============================================================
# (6а) χ²-критерий: ПРОСТАЯ гипотеза H0: Exp(λ0)
# ============================================================
cat("=== (6а) χ²-критерий: простая гипотеза H0: Exp(λ0 =", lambda0, ") ===\n")

brk   <- c(-Inf, 1.26, 2.46, 3.66, 6.06, 8.45, Inf)
l.brk <- length(brk)
lw    <- brk[1:(l.brk-1)]
up    <- brk[2:l.brk]

nu    <- hist(x, breaks = brk, plot = FALSE)$counts
pr    <- pexp(up, lambda0) - pexp(lw, lambda0)
npr   <- pr * n
res   <- (nu - npr) / sqrt(npr)
res2  <- res^2
X2    <- sum(res2)
ngr   <- length(nu)

csq.t <- data.frame(Group = 1:ngr, lower = lw, upper = up,
                    count = nu, prob = round(pr,4),
                    expr = round(npr,4), resid2 = round(res2,4))
print(csq.t)

df_simple <- ngr - 1
xal_simple <- qchisq(1 - alpha2, df_simple)
pv_simple  <- pchisq(X2, df_simple, lower.tail = FALSE)

cat("\nX² =", round(X2,4), " df =", df_simple,
    " крит. =", round(xal_simple,4), " p =", round(pv_simple,4), "\n")
cat("Вывод:", ifelse(X2 > xal_simple, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# ============================================================
# (6б) χ²-критерий: СЛОЖНАЯ гипотеза H0: Exp(λ̂)
#      Параметр λ оценивается минимизацией χ²
# ============================================================
cat("=== (6б) χ²-критерий: сложная гипотеза H0: Exp(λ̂) ===\n")

prob.exp <- function(par) pexp(up, par[1]) - pexp(lw, par[1])

csq.stat <- function(par) {
  ex <- n * prob.exp(par)
  ex[ex <= 0] <- 1e-10
  sum((nu - ex)^2 / ex)
}

qq    <- nlm(csq.stat, lambda_mle)
X2C   <- qq$minimum
lam.hat.csq <- qq$estimate[1]

df_complex <- ngr - 2
xal_complex <- qchisq(1 - alpha2, df_complex)
pv_complex  <- pchisq(X2C, df_complex, lower.tail = FALSE)

pr_c   <- prob.exp(lam.hat.csq)
npr_c  <- n * pr_c
res_c  <- (nu - npr_c) / sqrt(pmax(npr_c, 1e-10))
res2_c <- res_c^2

csq.t2 <- data.frame(Group = 1:ngr, lower = lw, upper = up,
                     count = nu, prob = round(pr_c,4),
                     expr = round(npr_c,4), resid2 = round(res2_c,4))
print(csq.t2)

cat("\nОценка λ̂ (мин. χ²) =", round(lam.hat.csq, 6), "\n")
cat("X²min =", round(X2C,4), " df =", df_complex,
    " крит. =", round(xal_complex,4), " p =", round(pv_complex,4), "\n")
cat("Вывод:", ifelse(X2C > xal_complex, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# ============================================================
# (7) НМК Неймана–Пирсона: H0: λ=λ0  vs  H1: λ=λ1
#     Exp(λ): T = Σxᵢ ~ Gamma(n, λ)
#     λ1=0.33 > λ0=0.17 => E[T|H1] < E[T|H0]
#     => отношение правдоподобий убывает в T
#     => критическая область ЛЕВОСТОРОННЯЯ: T ≤ c
# ============================================================
cat("=== (7) НМК: H0: λ=λ0=", lambda0, "  vs  H1: λ=λ1=", lambda1, " ===\n")
cat("(λ1 > λ0 => меньшее среднее => T=Σx малое => левосторонний критерий)\n\n")

T_obs <- sum(x)
cat("T_набл = Σxᵢ =", round(T_obs, 4), "\n")
cat("T|H0 ~ Gamma(n=50, rate=λ0) => E[T|H0] =", round(n/lambda0, 2), "\n\n")

# Критическое значение: P(T ≤ c | H0) ≤ α2
c_left <- qgamma(alpha2, shape = n, rate = lambda0)
real_alpha <- pgamma(c_left, shape = n, rate = lambda0)

cat("Критическое значение c =", round(c_left, 4), "\n")
cat("Реальный уровень значимости α* =", round(real_alpha, 4), "\n")
cat("T_набл =", round(T_obs,2), "vs c =", round(c_left,2), "\n")
cat("Вывод:", ifelse(T_obs <= c_left, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# Мощность
power <- pgamma(c_left, shape = n, rate = lambda1)
cat("Мощность = P(T ≤ c | H1) =", round(power, 6), "\n\n")

# p-value
p_val_nmk <- pgamma(T_obs, shape = n, rate = lambda0)
cat("p-value (левосторонний) =", round(p_val_nmk, 6), "\n\n")

# --- Если поменять H0 и H1 местами ---
cat("--- Меняем H0 и H1 местами ---\n")
cat("H0: λ=λ1=", lambda1, "  vs  H1: λ=λ0=", lambda0, "\n")
cat("(λ0 < λ1 => E[T|H1] > E[T|H0] => правосторонний критерий: T ≥ c')\n")

c_right     <- qgamma(1 - alpha2, shape = n, rate = lambda1)
real_alpha2 <- 1 - pgamma(c_right, shape = n, rate = lambda1)
power2      <- 1 - pgamma(c_right, shape = n, rate = lambda0)
p_val_nmk2  <- 1 - pgamma(T_obs, shape = n, rate = lambda1)

cat("c' =", round(c_right, 4), "\n")
cat("Реальный уровень α* =", round(real_alpha2, 4), "\n")
cat("T_набл =", round(T_obs,2), "vs c' =", round(c_right,2), "\n")
cat("Вывод:", ifelse(T_obs >= c_right, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n")
cat("Мощность =", round(power2, 6), "\n")
cat("p-value (правосторонний) =", round(p_val_nmk2, 6), "\n\n")

cat("=== РЕШЕНИЕ ЗАВЕРШЕНО ===\n")
