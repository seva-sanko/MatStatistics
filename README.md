# MatStatistics

Mathematical statistics assignments (variant 14) — hypothesis testing, distribution fitting, and regression analysis in R.

## IDZ3 — Hypothesis Testing & Distribution Fitting

Tasks on Poisson and Normal distributions: goodness-of-fit tests, empirical CDF, histograms, power analysis.

**Variant parameters:** α=0.05, λ₀=1.40, λ₁=0.70 (Poisson); a=0.00, b=1.20 (Normal)

```
ИДЗ3/
  files/
    task1_poisson.R      — Poisson distribution: chi-squared test, ECDF, NMK power
    task2_normal.R       — Normal distribution: chi-squared test, ECDF, Q-Q plot
    idz3_exp.R           — experiment runner
  generate_plots.R       — generates all plots into plots/
  report_with_plots.tex  — main LaTeX report
  report_var14.tex       — variant-specific report
  plots/                 — generated PNG plots
```

Tests performed:
- Chi-squared goodness-of-fit (simple and composite hypotheses)
- Kolmogorov-Smirnov via ECDF
- NMK power analysis at varying sample sizes

## IDZ4 — Regression Analysis

Linear and polynomial regression, ANOVA, interaction effects.

**Variant parameters:** variant 14, α=0.10

```
IDZ4_var14/
  IDZ4_var14_part1.R    — regression models (linear, polynomial)
  IDZ4_var14_part2.R    — ANOVA, interaction effects, residual analysis
  IDZ4_var14_report.tex — LaTeX report
  plot_A_linear.png     — linear regression fit
  plot_B_poly.png       — polynomial regression fit
  plot_C_hist.png       — residual histogram
  plot2_B_interaction.png
  plot2_C_hist.png
```

## Requirements

- R 4.x
- Packages: `ggplot2`, `stats`, `car`
