# Before / after ageing comparisons between packaging materials.

#' Before vs after ageing, per material, at one temperature.
#'
#' "Before" = jours == 0, "after" = longest ageing time. Paired t-test when both
#' groups have the same size (replicate i before is assumed to be the same sample
#' as replicate i after, in row order); Welch t-test otherwise.
#' @param var column name (e.g. "dpv", "pH", "reserve_alcaline_mol_kg")
before_after <- function(df, var, temp_c) {
  d <- df[df$temperature_C == temp_c, ]
  t_end <- max(d$jours)
  rows <- lapply(split(d, d$materiau), function(g) {
    b <- g[[var]][g$jours == 0]
    a <- g[[var]][g$jours == t_end]
    paired <- length(a) == length(b)
    tt <- t.test(a, b, paired = paired)
    data.frame(materiau = g$materiau[1], variable = var, days = t_end,
               mean_before = mean(b), mean_after = mean(a),
               diff = mean(a) - mean(b), ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
               p_value = tt$p.value, test = if (paired) "paired t" else "Welch t")
  })
  do.call(rbind, unname(rows))
}

#' One-way ANOVA + Tukey HSD on the change (after - before) between materials.
#' The change is computed per replicate, pairing by row order (see before_after).
compare_materials <- function(df, var, temp_c) {
  d <- df[df$temperature_C == temp_c, ]
  t_end <- max(d$jours)
  delta <- do.call(rbind, unname(lapply(split(d, d$materiau), function(g) {
    data.frame(materiau = g$materiau[1],
               delta = g[[var]][g$jours == t_end] - g[[var]][g$jours == 0])
  })))
  fit <- aov(delta ~ materiau, data = delta)
  list(anova_p = summary(fit)[[1]][["Pr(>F)"]][1],
       tukey = as.data.frame(TukeyHSD(fit)$materiau))
}
