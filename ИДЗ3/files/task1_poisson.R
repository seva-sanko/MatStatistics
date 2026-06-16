# ============================================================
#  Вариант 14  —  Задание 1: Распределение Пуассона
#  Параметры: alpha1=0.05, a=0.00, b=1.20, lambda0=1.40, lambda1=0.70
# ============================================================

# ---------- Исходные данные ----------
x <- c(1,0,1,1,3,1,0,1,0,1,1,1,0,1,1,0,0,1,2,0,
       1,1,0,1,1,0,2,0,0,0,0,1,1,0,0,2,0,1,0,1,
       1,0,1,0,0,0,1,2,0,2)
n    <- length(x)
alpha1  <- 0.05
a_par   <- 0.00
b_par   <- 1.20
lambda0 <- 1.40
lambda1 <- 0.70

cat("=== Исходные данные ===\n")
cat("n =", n, "\n")
cat("Наблюдения:", x, "\n\n")

# ============================================================
# (a) Вариационный ряд, эмпирическая функция распределения,
#     гистограмма частот
# ============================================================
cat("=== (a) Вариационный ряд ===\n")
var_row <- sort(x)
cat("Вариационный ряд:\n", var_row, "\n\n")

# Таблица частот
freq_table <- table(x)
rel_freq   <- freq_table / n
cum_freq   <- cumsum(rel_freq)

cat("Таблица частот:\n")
print(data.frame(
  x_i       = as.integer(names(freq_table)),
  n_i       = as.integer(freq_table),
  rel_n_i   = round(as.numeric(rel_freq), 4),
  F_emp     = round(as.numeric(cum_freq), 4)
))
cat("\n")

# Эмпирическая функция распределения (график)
vals <- as.integer(names(freq_table))

par(mfrow = c(1, 2))

# Гистограмма частот
barplot(rel_freq,
        names.arg = names(rel_freq),
        main = "Гистограмма относительных частот",
        xlab = "x", ylab = "Частота",
        col  = "steelblue", border = "white")

# Эмпирическая функция распределения
plot(ecdf(x),
     main = "Эмпирическая функция распределения",
     xlab = "x", ylab = "F*(x)",
     col  = "darkred", lwd = 2)
grid()

par(mfrow = c(1, 1))

# ============================================================
# (b) Выборочные характеристики
# ============================================================
cat("=== (b) Выборочные характеристики ===\n")

x_mean <- mean(x)
x_var  <- var(x)           # несмещённая s^2
x_var_b <- mean((x - x_mean)^2)  # смещённая (момент)
x_med  <- median(x)

# Асимметрия (3-й центральный момент / s^3)
m3 <- mean((x - x_mean)^3)
m2 <- mean((x - x_mean)^2)
asymm <- m3 / m2^(3/2)

# Эксцесс (4-й момент / m2^2 - 3)
m4   <- mean((x - x_mean)^4)
kurt <- m4 / m2^2 - 3

# P(X in [a, b])
p_ab <- mean(x >= a_par & x <= b_par)

cat("(i)   Выборочное среднее (мат. ожидание): x̄ =", round(x_mean, 4), "\n")
cat("(ii)  Выборочная дисперсия (несмещ. s²):  s² =", round(x_var, 4), "\n")
cat("      Выборочная дисперсия (смещ. D):      D =", round(x_var_b, 4), "\n")
cat("(iii) Медиана:                             Me =", x_med, "\n")
cat("(iv)  Коэффициент асимметрии:              As =", round(asymm, 4), "\n")
cat("(v)   Эксцесс:                             Ex =", round(kurt, 4), "\n")
cat("(vi)  P(X ∈ [", a_par, ",", b_par, "]):       P* =", round(p_ab, 4), "\n\n")

# ============================================================
# (c) Оценки параметра λ: МПП и метод моментов; смещение
# ============================================================
cat("=== (c) Оценки λ ===\n")

# МПП для Пуассона: λ_mle = x̄
lambda_mle <- x_mean

# Метод моментов: E[X]=λ  =>  λ_mm = x̄  (совпадает)
lambda_mm  <- x_mean

# Смещение E[λ̂] - λ: теоретически обе оценки несмещённые
# Практическое (Bootstrap) смещение
set.seed(42)
B <- 5000
boot_mle <- replicate(B, mean(sample(x, n, replace = TRUE)))
bias_mle <- mean(boot_mle) - lambda_mle

cat("МПП:            λ̂_mle =", round(lambda_mle, 4), "\n")
cat("Метод моментов: λ̂_mm  =", round(lambda_mm,  4), "\n")
cat("(Обе оценки совпадают для Пуассона)\n")
cat("Bootstrap смещение МПП: bias ≈", round(bias_mle, 6),
    "(≈ 0 — оценка несмещённая)\n\n")

# ============================================================
# (d) Асимптотический доверительный интервал для λ (МПП)
# ============================================================
cat("=== (d) Доверительный интервал для λ ===\n")

z_alpha <- qnorm(1 - alpha1 / 2)   # z_{0.975}
se_lambda <- sqrt(lambda_mle / n)

ci_low  <- lambda_mle - z_alpha * se_lambda
ci_high <- lambda_mle + z_alpha * se_lambda

cat("z_{1-α/2} =", round(z_alpha, 4), "\n")
cat("SE(λ̂) =", round(se_lambda, 4), "\n")
cat("ДИ уровня", 1 - alpha1, ": [",
    round(ci_low, 4), ";", round(ci_high, 4), "]\n\n")

# ============================================================
# (e) Критерий χ² — простая гипотеза H0: Poisson(λ0)
# ============================================================
cat("=== (e) χ²-критерий: простая гипотеза H0: Poisson(λ0 =", lambda0, ") ===\n")

vals_all <- 0:max(x)
obs <- as.integer(table(factor(x, levels = vals_all)))
# Теоретические вероятности
p_theor <- dpois(vals_all, lambda = lambda0)

# Объединяем хвосты чтобы ожидаемые частоты >= 5
expected <- n * p_theor
# Группировка: объединяем 3 и выше
grp_obs <- c(obs[1], obs[2], obs[3], sum(obs[4:length(obs)]))
grp_exp <- c(expected[1], expected[2], expected[3], sum(expected[4:length(expected)]))
grp_lbl <- c("0", "1", "2", "≥3")

cat("Группы | Наблюд. | Ожидаем. (λ0)\n")
for (i in seq_along(grp_lbl)) {
  cat(sprintf("  %-5s |   %3d   |   %6.3f\n", grp_lbl[i], grp_obs[i], grp_exp[i]))
}

chi2_e <- sum((grp_obs - grp_exp)^2 / grp_exp)
df_e   <- length(grp_obs) - 1   # простая гипотеза: нет оцениваемых параметров
chi2_cr_e <- qchisq(1 - alpha1, df = df_e)
p_val_e   <- 1 - pchisq(chi2_e, df = df_e)

cat("\nχ² статистика =", round(chi2_e, 4), "\n")
cat("Степени свободы df =", df_e, "\n")
cat("Критическое значение χ²_{1-α}(df) =", round(chi2_cr_e, 4), "\n")
cat("p-значение =", round(p_val_e, 4), "\n")
if (chi2_e > chi2_cr_e) {
  cat("Вывод: ОТВЕРГАЕМ H0 на уровне", alpha1, "\n")
} else {
  cat("Вывод: НЕТ оснований отвергнуть H0 на уровне", alpha1, "\n")
}
cat("Наибольший уровень значимости (p-value) =", round(p_val_e, 4), "\n\n")

# ============================================================
# (f) Критерий χ² — сложная гипотеза H0: Poisson(λ) (λ оценивается)
# ============================================================
cat("=== (f) χ²-критерий: сложная гипотеза (λ оценивается) ===\n")

p_theor_f <- dpois(vals_all, lambda = lambda_mle)
expected_f <- n * p_theor_f

grp_obs_f <- c(obs[1], obs[2], obs[3], sum(obs[4:length(obs)]))
grp_exp_f <- c(expected_f[1], expected_f[2], expected_f[3],
               sum(expected_f[4:length(expected_f)]))

cat("Группы | Наблюд. | Ожидаем. (λ̂)\n")
for (i in seq_along(grp_lbl)) {
  cat(sprintf("  %-5s |   %3d   |   %6.3f\n", grp_lbl[i], grp_obs_f[i], grp_exp_f[i]))
}

chi2_f  <- sum((grp_obs_f - grp_exp_f)^2 / grp_exp_f)
df_f    <- length(grp_obs_f) - 1 - 1   # -1 за оценённый параметр λ
chi2_cr_f <- qchisq(1 - alpha1, df = df_f)
p_val_f   <- 1 - pchisq(chi2_f, df = df_f)

cat("\nχ² статистика =", round(chi2_f, 4), "\n")
cat("Степени свободы df =", df_f, "\n")
cat("Критическое значение χ²_{1-α}(df) =", round(chi2_cr_f, 4), "\n")
cat("p-значение =", round(p_val_f, 4), "\n")
if (chi2_f > chi2_cr_f) {
  cat("Вывод: ОТВЕРГАЕМ H0 на уровне", alpha1, "\n")
} else {
  cat("Вывод: НЕТ оснований отвергнуть H0 на уровне", alpha1, "\n")
}
cat("Наибольший уровень значимости (p-value) =", round(p_val_f, 4), "\n\n")

# ============================================================
# (g) Наиболее мощный критерий: H0: λ=λ0  vs  H1: λ=λ1
#     (λ1 < λ0 => левосторонний критерий по достаточной статистике T=sum(x))
# ============================================================
cat("=== (g) Наиболее мощный критерий (НМК) Неймана-Пирсона ===\n")
cat("H0: λ = λ0 =", lambda0, "  vs  H1: λ = λ1 =", lambda1, "\n")

T_obs <- sum(x)
cat("Достаточная статистика T = Σxi =", T_obs, "\n")
cat("T ~ Poisson(n*λ) при H0: n*λ0 =", n * lambda0, "\n\n")

# λ1 < λ0 => ОО: T <= c (левостороннее)
# Находим критическое значение c: P(T <= c | H0) <= alpha1
c_crit_g <- qpois(alpha1, lambda = n * lambda0)
real_alpha_g <- ppois(c_crit_g, lambda = n * lambda0)

cat("Критическое значение c =", c_crit_g, "\n")
cat("Реальный уровень значимости α* = P(T ≤ c | H0) =", round(real_alpha_g, 4), "\n")
cat("T_наблюд =", T_obs, " > c =", c_crit_g, "=> ",
    ifelse(T_obs <= c_crit_g, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# Мощность критерия
power_g <- ppois(c_crit_g, lambda = n * lambda1)
cat("Мощность = P(T ≤ c | H1: λ=λ1) =", round(power_g, 4), "\n\n")

# Если поменять местами H0 и H1:
cat("--- Если поменять H0 и H1 ---\n")
cat("H0: λ = λ1 =", lambda1, "  vs  H1: λ = λ0 =", lambda0, "\n")
# Теперь λ0 > λ1 => правосторонний критерий
c_crit_g2 <- qpois(1 - alpha1, lambda = n * lambda1)
real_alpha_g2 <- 1 - ppois(c_crit_g2 - 1, lambda = n * lambda1)
cat("Критическое значение c' =", c_crit_g2, "\n")
cat("Реальный уровень α* =", round(real_alpha_g2, 4), "\n")
cat("T_наблюд =", T_obs, ":",
    ifelse(T_obs > c_crit_g2, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# ============================================================
# (h) То же, но для геометрического распределения
#     P_λ(X=k) = λ^k / (λ+1)^(k+1),  k = 0,1,2,...
#     E[X] = λ,  оценка λ̂ = x̄  (метод моментов = МПП)
# ============================================================
cat("=== (h) Геометрическое распределение P_λ(X=k) = λ^k/(λ+1)^(k+1) ===\n")

dgeom_custom <- function(k, lam) lam^k / (lam + 1)^(k + 1)

# (c) Оценка λ
lambda_geom_mle <- x_mean   # E[X]=λ => МПП и ММ совпадают
cat("(c) λ̂_geom (МПП = ММ) =", round(lambda_geom_mle, 4), "\n")
cat("    Смещение: теоретически 0 (несмещённая оценка)\n\n")

# (d) ДИ для λ геометрического через CLT (дисперсия D=λ(λ+1))
se_geom <- sqrt(lambda_geom_mle * (lambda_geom_mle + 1) / n)
ci_geom_low  <- lambda_geom_mle - z_alpha * se_geom
ci_geom_high <- lambda_geom_mle + z_alpha * se_geom
cat("(d) ДИ для λ_geom уровня", 1 - alpha1, ": [",
    round(ci_geom_low, 4), ";", round(ci_geom_high, 4), "]\n\n")

# (e) χ²: простая гипотеза H0: Geom(λ0)
cat("(e) χ²-критерий: простая гипотеза H0: Geom(λ0 =", lambda0, ")\n")
p_theor_h0 <- dgeom_custom(vals_all, lambda0)
p_theor_h0[length(p_theor_h0)] <- 1 - sum(p_theor_h0[-length(p_theor_h0)])  # хвост
exp_h0 <- n * p_theor_h0
grp_obs_h <- c(obs[1], obs[2], obs[3], sum(obs[4:length(obs)]))
grp_exp_h0 <- c(exp_h0[1], exp_h0[2], exp_h0[3], sum(exp_h0[4:length(exp_h0)]))

chi2_h0 <- sum((grp_obs_h - grp_exp_h0)^2 / grp_exp_h0)
df_h0   <- length(grp_obs_h) - 1
chi2_cr_h0 <- qchisq(1 - alpha1, df = df_h0)
p_h0    <- 1 - pchisq(chi2_h0, df = df_h0)
cat("  χ² =", round(chi2_h0,4), " df =", df_h0,
    " крит. =", round(chi2_cr_h0,4), " p =", round(p_h0,4), "\n")
cat("  Вывод:", ifelse(chi2_h0 > chi2_cr_h0, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# (f) χ²: сложная гипотеза H0: Geom(λ̂)
cat("(f) χ²-критерий: сложная гипотеза H0: Geom(λ̂ =", round(lambda_geom_mle,4), ")\n")
p_theor_hf <- dgeom_custom(vals_all, lambda_geom_mle)
p_theor_hf[length(p_theor_hf)] <- 1 - sum(p_theor_hf[-length(p_theor_hf)])
exp_hf <- n * p_theor_hf
grp_exp_hf <- c(exp_hf[1], exp_hf[2], exp_hf[3], sum(exp_hf[4:length(exp_hf)]))

chi2_hf <- sum((grp_obs_h - grp_exp_hf)^2 / grp_exp_hf)
df_hf   <- length(grp_obs_h) - 1 - 1
chi2_cr_hf <- qchisq(1 - alpha1, df = df_hf)
p_hf    <- 1 - pchisq(chi2_hf, df = df_hf)
cat("  χ² =", round(chi2_hf,4), " df =", df_hf,
    " крит. =", round(chi2_cr_hf,4), " p =", round(p_hf,4), "\n")
cat("  Вывод:", ifelse(chi2_hf > chi2_cr_hf, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n\n")

# (g) НМК для геометрического: H0: λ=λ0 vs H1: λ=λ1
cat("(g) НМК для геометрического: H0: λ=λ0 =", lambda0, " vs H1: λ=λ1 =", lambda1, "\n")
# T = Σxi является достаточной статистикой
# λ1 < λ0 => левостороннее ОО: T <= c
# T|H0: сумма n независимых геом., E=nλ0, Var=nλ0(λ0+1)
# Используем CLT для нахождения c
mu_T_h0  <- n * lambda0
sd_T_h0  <- sqrt(n * lambda0 * (lambda0 + 1))
c_crit_geom <- floor(mu_T_h0 + qnorm(alpha1) * sd_T_h0)
real_alpha_geom <- pnorm((c_crit_geom - mu_T_h0) / sd_T_h0)
cat("  c (CLT) =", c_crit_geom, " реальный α* ≈", round(real_alpha_geom,4), "\n")
cat("  T_наблюд =", T_obs, ":",
    ifelse(T_obs <= c_crit_geom, "ОТВЕРГАЕМ H0", "НЕ ОТВЕРГАЕМ H0"), "\n")
mu_T_h1 <- n * lambda1
sd_T_h1 <- sqrt(n * lambda1 * (lambda1 + 1))
power_geom <- pnorm((c_crit_geom - mu_T_h1) / sd_T_h1)
cat("  Мощность ≈", round(power_geom, 4), "\n\n")

cat("=== РЕШЕНИЕ ЗАВЕРШЕНО ===\n")
