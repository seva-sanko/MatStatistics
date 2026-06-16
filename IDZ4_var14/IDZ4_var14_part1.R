# ИДЗ 4, Вариант 14 (513201226)
# Часть 1: Регрессионный анализ

# ===== ДАННЫЕ =====
alpha <- 0.10
h <- 1.80

Y <- c(11.83, 15.46, 18.12, 14.43, 11.32, 14.44, 14.26, 10.10, 21.75, 17.73,
       14.21, 8.80, 16.14, 9.24, 16.23, 14.48, 15.63, 15.02, 14.07, 15.80,
       17.78, 16.76, 11.90, 16.83, 12.48, 16.53, 10.09, 19.93, 12.73, 16.56,
       18.11, 17.59, 11.02, 15.09, 12.08, 19.75, 18.15, 14.99, 14.53, 19.69,
       13.29, 13.62, 20.22, 12.08, 15.13, 16.99, 18.30, 13.71, 14.84, 13.92)

X <- c(3, 2, 1, 3, 1, 3, 1, 3, 2, 2, 2, 3, 1, 0, 2, 0, 1, 0, 2, 0, 1, 1, 1,
       2, 1, 1, 3, 3, 2, 3, 3, 2, 2, 0, 3, 1, 3, 3, 3, 0, 3, 3, 2, 1, 0, 3,
       0, 3, 1, 2)

n <- length(Y)
dat <- data.frame(Y = Y, X = X)

cat("====================================================\n")
cat("ИДЗ 4, Вариант 14 | ЧАСТЬ 1: Регрессионный анализ\n")
cat("====================================================\n\n")
cat("n =", n, "| alpha =", alpha, "| h =", h, "\n\n")


# ===== ЗАДАНИЕ A: Линейная модель =====
cat("---------- ЗАДАНИЕ A: Линейная регрессия ----------\n")

x0 <- array(1, dim = n)
Xl <- matrix(c(x0, X), nrow = 2, ncol = n, byrow = TRUE)
Sl <- Xl %*% t(Xl)
S1l <- solve(Sl)
bhatl <- S1l %*% Xl %*% as.matrix(Y)

cat("МНК оценки линейной модели Y = beta1 + beta2*X + eps:\n")
cat("  beta1_hat =", round(bhatl[1], 7), "\n")
cat("  beta2_hat =", round(bhatl[2], 7), "\n\n")

# График линейной модели
png("/home/claude/plot_A_linear.png", width = 800, height = 600)
plot(X, Y, cex = 1.5, pch = 1,
     main = "Линейная модель",
     xlab = "x", ylab = "y",
     col = "black")
x1 <- c(min(X), max(X))
y1 <- bhatl[1] + bhatl[2] * x1
lines(x1, y1, col = "blue", lwd = 2)
legend("topright", legend = c("Данные", "Линейная регрессия"),
       col = c("black", "blue"), pch = c(1, NA), lty = c(NA, 1), lwd = 2)
dev.off()
cat("График линейной модели сохранен: plot_A_linear.png\n\n")

# Через lm для проверки
q1 <- lm(Y ~ X, data = dat)
q1s <- summary(q1)


# ===== ЗАДАНИЕ B: Полиномиальная модель =====
cat("---------- ЗАДАНИЕ B: Полиномиальная (квадратичная) модель ----------\n")

X2 <- matrix(c(x0, X, X^2), nrow = 3, ncol = n, byrow = TRUE)
S2 <- X2 %*% t(X2)
S21 <- solve(S2)
bhat2 <- S21 %*% X2 %*% as.matrix(Y)

cat("МНК оценки квадратичной модели Y = beta1 + beta2*X + beta3*X^2 + eps:\n")
cat("  beta1_hat =", round(bhat2[1], 7), "\n")
cat("  beta2_hat =", round(bhat2[2], 7), "\n")
cat("  beta3_hat =", round(bhat2[3], 7), "\n\n")

# График сравнения линейной и полиномиальной модели
png("/home/claude/plot_B_poly.png", width = 800, height = 600)
plot(X, Y, cex = 1.5, pch = 1,
     main = "Сравнение моделей",
     xlab = "x", ylab = "y",
     col = "black")
lines(x1, y1, col = "blue", lwd = 2)
x2_seq <- seq(min(X), max(X), length.out = 1001)
y2_seq <- bhat2[1] + bhat2[2] * x2_seq + bhat2[3] * x2_seq^2
lines(x2_seq, y2_seq, col = "red", lwd = 2)
legend("topright",
       legend = c("Данные", "Линейная", "Полиномиальная"),
       col = c("black", "blue", "red"),
       pch = c(1, NA, NA), lty = c(NA, 1, 1), lwd = 2)
dev.off()
cat("График полиномиальной модели сохранен: plot_B_poly.png\n\n")

# Через lm
q2 <- lm(Y ~ I(X) + I(X^2), data = dat)
q2s <- summary(q2)

B_cov <- q2s$cov.unscaled
V_cov <- q2s$sigma^2 * B_cov


# ===== ЗАДАНИЕ C: Гистограмма остатков, хи-квадрат =====
cat("---------- ЗАДАНИЕ C: Проверка нормальности остатков ----------\n")

res <- as.numeric(Y - t(X2) %*% bhat2)
res_sorted <- sort(res)
cat("Упорядоченные остатки:\n")
print(round(res_sorted, 7))
cat("\n")

# Гистограмма с шагом h
min_val <- min(res)
max_val <- max(res)
brks <- seq(from = min_val,
            to   = ceiling(max_val / h) * h,
            by   = h)

png("/home/claude/plot_C_hist.png", width = 800, height = 600)
hh <- hist(res, breaks = brks, probability = TRUE,
           main = "Гистограмма остатков для полиномиальной модели",
           xlab = "Values", col = "lightgray", border = "white")
curve(dnorm(x, mean = mean(res), sd = sd(res)),
      add = TRUE, col = "red", lwd = 2)
dev.off()

cat("Параметры гистограммы:\n")
cat("  breaks:", round(hh$breaks, 7), "\n")
cat("  counts:", hh$counts, "\n\n")

# Тест хи-квадрат (сложная гипотеза, оцениваем sigma)
brk2    <- hh$breaks
observed <- hh$counts
expected_probs <- diff(pnorm(brk2, mean = 0, sd = sd(res)))
expected <- n * expected_probs

# Нелинейная оптимизация: минимизируем хи-квадрат по sigma
chi2_func <- function(sigma) {
  ep <- diff(pnorm(brk2, mean = 0, sd = sigma))
  ex <- n * ep
  # избегаем деления на ноль
  idx <- ex > 0
  sum((observed[idx] - ex[idx])^2 / ex[idx])
}
opt <- optimize(chi2_func, interval = c(0.1, 20))
sigma_hat <- opt$minimum
chi2_stat <- opt$objective

k      <- length(observed)          # число интервалов
df_chi <- k - 1 - 1                 # сложная гипотеза: оцениваем 1 параметр (sigma)
critical_value <- qchisq(1 - alpha, df_chi)
p_value <- pchisq(chi2_stat, df_chi, lower.tail = FALSE)

cat("Тест хи-квадрат (нормальность остатков):\n")
cat("  Оптимальная sigma:", round(sigma_hat, 6), "\n")
cat("  Статистика X^2 =", round(chi2_stat, 6), "\n")
cat("  Число интервалов k =", k, "\n")
cat("  Степени свободы df =", df_chi, "(k-1-1)\n")
cat("  Критическое значение chi^2_{", df_chi, ", ", 1-alpha, "} =",
    round(critical_value, 6), "\n")
cat("  P-значение:", round(p_value, 6), "\n")
if (chi2_stat > critical_value) {
  cat("  Вывод: Отвергаем H0 — отклонение от нормальности значимо\n\n")
} else {
  cat("  Вывод: Не отвергаем H0 — данные совместимы с нормальностью\n\n")
}


# ===== ЗАДАНИЕ D: Доверительные интервалы =====
cat("---------- ЗАДАНИЕ D: Доверительные интервалы для beta2, beta3 ----------\n")

ci_q2 <- confint(q2, level = 1 - alpha)
cat("Доверительные интервалы (уровень доверия", 1 - alpha, "):\n")
cat("  beta1: [", round(ci_q2[1, 1], 7), ";", round(ci_q2[1, 2], 7), "]\n")
cat("  beta2: [", round(ci_q2[2, 1], 7), ";", round(ci_q2[2, 2], 7), "]\n")
cat("  beta3: [", round(ci_q2[3, 1], 7), ";", round(ci_q2[3, 2], 7), "]\n\n")

# Совместный (эллиптический) доверительный регион для beta2 и beta3
# через F-статистику (q=2 параметра)
# Интерпретация
cat("Интерпретация:\n")
for (i in 1:3) {
  nm <- c("beta1", "beta2", "beta3")[i]
  lo <- ci_q2[i, 1]; hi <- ci_q2[i, 2]
  if (lo > 0 || hi < 0) {
    cat(sprintf("  %s: [%.4f; %.4f] — не включает 0, параметр значим\n", nm, lo, hi))
  } else {
    cat(sprintf("  %s: [%.4f; %.4f] — включает 0, параметр не значим\n", nm, lo, hi))
  }
}
cat("\n")


# ===== ЗАДАНИЕ E: Проверка гипотез линейности и независимости =====
cat("---------- ЗАДАНИЕ E: Гипотезы линейности и независимости ----------\n")

q0 <- lm(Y ~ 1, data = dat)

# H0: зависимость линейна (beta3 = 0) — сравниваем q1 и q2
test_lin <- anova(q1, q2)
F_lin    <- test_lin$F[2]
df1_lin  <- test_lin$Df[2]
df2_lin  <- test_lin$Res.Df[2]
F_crit_lin <- qf(1 - alpha, df1_lin, df2_lin)
p_lin    <- test_lin$`Pr(>F)`[2]

cat("Гипотеза линейности H0: beta3 = 0\n")
cat("  F-статистика =", round(F_lin, 7), "\n")
cat("  Критическое значение F_{", df1_lin, ",", df2_lin, ",", 1-alpha, "} =",
    round(F_crit_lin, 7), "\n")
cat("  P-значение:", round(p_lin, 7), "\n")
if (F_lin > F_crit_lin) {
  cat("  Вывод: Отвергаем H0 — зависимость нелинейна\n\n")
} else {
  cat("  Вывод: Не отвергаем H0 — зависимость линейна\n\n")
}

# H0: независимость (beta2 = beta3 = 0) — сравниваем q0 и q2
test_ind <- anova(q0, q2)
F_ind    <- test_ind$F[2]
df1_ind  <- test_ind$Df[2]
df2_ind  <- test_ind$Res.Df[2]
F_crit_ind <- qf(1 - alpha, df1_ind, df2_ind)
p_ind    <- test_ind$`Pr(>F)`[2]

cat("Гипотеза независимости H0: beta2 = 0, beta3 = 0\n")
cat("  F-статистика =", round(F_ind, 7), "\n")
cat("  Критическое значение F_{", df1_ind, ",", df2_ind, ",", 1-alpha, "} =",
    round(F_crit_ind, 7), "\n")
cat("  P-значение:", round(p_ind, 7), "\n")
if (F_ind > F_crit_ind) {
  cat("  Вывод: Отвергаем H0 — Y зависит от X\n\n")
} else {
  cat("  Вывод: Не отвергаем H0 — Y не зависит от X\n\n")
}


# ===== ЗАДАНИЕ F: Выбор модели по AIC и BIC =====
cat("---------- ЗАДАНИЕ F: Сравнение моделей по AIC и BIC ----------\n")

aic_q0 <- AIC(q0)
aic_q1 <- AIC(q1)
aic_q2 <- AIC(q2)

bic_q0 <- BIC(q0)
bic_q1 <- BIC(q1)
bic_q2 <- BIC(q2)

model_comp <- data.frame(
  Model = c("(1) Линейная (q1)", "(2) Квадратичная (q2)", "(3) Константа (q0)"),
  AIC = round(c(aic_q1, aic_q2, aic_q0), 4),
  BIC = round(c(bic_q1, bic_q2, bic_q0), 4)
)

cat("Таблица сравнения моделей:\n")
print(model_comp)
cat("\n")

best_aic <- model_comp$Model[which.min(model_comp$AIC)]
best_bic <- model_comp$Model[which.min(model_comp$BIC)]
cat("Лучшая модель по AIC:", best_aic, "\n")
cat("Лучшая модель по BIC:", best_bic, "\n\n")


# ===== ЗАДАНИЕ G: Итоговый вывод =====
cat("---------- ЗАДАНИЕ G: Интерпретация результатов ----------\n")
cat("1. Линейная модель: Y =", round(bhatl[1], 4), "+",
    round(bhatl[2], 4), "* X\n")
cat("2. Квадратичная модель: Y =", round(bhat2[1], 4), "+",
    round(bhat2[2], 4), "* X +", round(bhat2[3], 4), "* X^2\n")
cat("3. Нормальность остатков:", ifelse(chi2_stat <= critical_value, "подтверждена", "отвергнута"), "\n")
cat("4. Зависимость линейна:", ifelse(F_lin <= F_crit_lin, "Да (H0 принята)", "Нет (H0 отвергнута)"), "\n")
cat("5. Y зависит от X:", ifelse(F_ind > F_crit_ind, "Да (H0 отвергнута)", "Нет (H0 принята)"), "\n")
cat("6. Наилучшая модель по AIC:", best_aic, "\n")
cat("7. Наилучшая модель по BIC:", best_bic, "\n")
