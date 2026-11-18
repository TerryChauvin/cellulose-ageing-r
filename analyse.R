# Full pipeline on the synthetic example data. Run from the repository root:
#   Rscript analyse.R
suppressPackageStartupMessages(library(ggplot2))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

df <- read.csv("data/exemple_vieillissement.csv")
df$dpv <- dpv_from_viscosity(df$viscosite_intrinseque_mL_g)

# Kinetics: Ekenstam k per material and temperature, Arrhenius per material
kin <- ekenstam_by_group(df)
arr <- do.call(rbind, unname(lapply(split(kin, kin$materiau), function(g)
  cbind(materiau = g$materiau[1], arrhenius_fit(g$k, g$temperature_C)))))
write.csv(kin, "outputs/ekenstam_k.csv", row.names = FALSE)
write.csv(arr, "outputs/arrhenius.csv", row.names = FALSE)

# Before / after at 80 C and between-material comparison
ref <- 80
vars <- c("dpv", "pH", "reserve_alcaline_mol_kg")
summary_tab <- do.call(rbind, lapply(vars, function(v) before_after(df, v, ref)))
write.csv(summary_tab, "outputs/summary_before_after.csv", row.names = FALSE)
print(summary_tab, digits = 3)
print(kin, digits = 3)
print(arr, digits = 3)
for (v in vars)
  cat(v, ": ANOVA on change between materials, p =",
      signif(compare_materials(df, v, ref)$anova_p, 3), "\n")

# Figures
p1 <- ggplot(df[df$temperature_C == ref, ], aes(jours, dpv, colour = materiau)) +
  geom_point(alpha = 0.6) + stat_summary(fun = mean, geom = "line") +
  labs(x = "Ageing time (days)", y = "DPv", title = "DPv during ageing at 80 C (synthetic data)")
ggsave("figures/dpv_vs_time.png", p1, width = 7, height = 4, dpi = 150)

df$ekenstam <- NA_real_
for (m in unique(df$materiau)) for (tc in unique(df$temperature_C)) {
  i <- df$materiau == m & df$temperature_C == tc
  df$ekenstam[i] <- 1 / df$dpv[i] - 1 / mean(df$dpv[i & df$jours == 0])
}
p2 <- ggplot(df, aes(jours, ekenstam, colour = materiau, shape = factor(temperature_C))) +
  geom_point() +
  geom_smooth(aes(group = interaction(materiau, temperature_C)),
              method = "lm", formula = y ~ 0 + x, se = FALSE, linewidth = 0.4) +
  labs(x = "Ageing time (days)", y = "1/DP - 1/DP0", shape = "T (C)",
       title = "Ekenstam plot (synthetic data)")
ggsave("figures/ekenstam.png", p2, width = 7, height = 4, dpi = 150)

p3 <- ggplot(kin, aes(1000 / (temperature_C + 273.15), k, colour = materiau)) +
  geom_point() + geom_line() + scale_y_log10() +
  labs(x = "1000 / T (1/K)", y = "k (1/day)", title = "Arrhenius plot (synthetic data)")
ggsave("figures/arrhenius.png", p3, width = 7, height = 4, dpi = 150)

d80 <- df[df$temperature_C == ref, ]
long <- rbind(
  data.frame(materiau = d80$materiau, jours = d80$jours, variable = "pH", value = d80$pH),
  data.frame(materiau = d80$materiau, jours = d80$jours, variable = "reserve_alcaline_mol_kg",
             value = d80$reserve_alcaline_mol_kg))
p4 <- ggplot(long, aes(factor(jours), value, fill = materiau)) + geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Ageing time (days)", y = NULL,
       title = "pH and alkaline reserve at 80 C (synthetic data)")
ggsave("figures/ph_alkaline_reserve.png", p4, width = 8, height = 4, dpi = 150)
