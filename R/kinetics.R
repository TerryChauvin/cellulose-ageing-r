# Degradation kinetics of cellulose: viscometric DP, Ekenstam, chain scissions, Arrhenius.

#' Viscometric degree of polymerisation from intrinsic viscosity.
#'
#' Mark-Houwink form used by ISO 5351 (cellulose in cupriethylenediamine, CED):
#' DP^a = K * [eta], with [eta] in mL/g, K = 0.75 and a = 0.905.
#' Constants as reported by Henniges et al. (2011, Cellulose 18) for the ISO 5351
#' procedure; check them against the edition of the standard you apply.
#' Older practice (e.g. TAPPI T230 / Evans & Wallis) uses other constants: pass K, a.
#'
#' @param eta intrinsic viscosity [eta] in mL/g (numeric vector)
#' @return DPv (numeric vector)
dpv_from_viscosity <- function(eta, K = 0.75, a = 0.905) {
  if (any(eta < 0, na.rm = TRUE)) stop("eta must be >= 0")
  (K * eta)^(1 / a)
}

#' Number of chain scissions per initial chain: DP0 / DP - 1.
chain_scissions <- function(dp, dp0) dp0 / dp - 1

#' Ekenstam pseudo-first-order fit: 1/DP(t) - 1/DP0 = k * t (regression through the origin).
#'
#' @param days ageing time (days), values > 0 are used
#' @param dp DPv measured at those times
#' @param dp0 DPv before ageing (scalar, e.g. mean of the t = 0 samples)
#' @return one-row data.frame: k (1/day), 95 % CI, n
ekenstam_fit <- function(days, dp, dp0, level = 0.95) {
  keep <- days > 0 & is.finite(dp)
  if (sum(keep) < 3) stop("need at least 3 points with days > 0")
  y <- 1 / dp[keep] - 1 / dp0
  t <- days[keep]
  fit <- lm(y ~ 0 + t)
  ci <- confint(fit, level = level)[1, ]
  data.frame(k = unname(coef(fit)), k_low = ci[[1]], k_high = ci[[2]], n = sum(keep))
}

#' Ekenstam fit for each (material, temperature) of a data.frame that has
#' columns materiau, temperature_C, jours, dpv. DP0 is the mean DPv at jours == 0
#' for the same material and temperature.
ekenstam_by_group <- function(df) {
  out <- lapply(split(df, list(df$materiau, df$temperature_C), drop = TRUE), function(g) {
    dp0 <- mean(g$dpv[g$jours == 0])
    cbind(materiau = g$materiau[1], temperature_C = g$temperature_C[1],
          dp0 = dp0, ekenstam_fit(g$jours, g$dpv, dp0))
  })
  do.call(rbind, unname(out))
}

#' Arrhenius fit: ln k = ln A - Ea / (R T). Needs >= 3 temperatures for a CI.
#'
#' @param k rate constants, same unit for all
#' @param temperature_C temperatures in degrees Celsius
#' With exactly 3 temperatures the fit has 1 residual degree of freedom, so the 95 % CI uses
#' t(0.975, 1) = 12.7 and is very wide (expected, not a bug): use more temperatures to tighten it.
#' @return one-row data.frame: Ea (kJ/mol) with 95 % CI, ln A
arrhenius_fit <- function(k, temperature_C, level = 0.95) {
  if (length(unique(temperature_C)) < 3) stop("need >= 3 distinct temperatures")
  R_gas <- 8.314462618 # J/(mol K)
  inv_T <- 1 / (temperature_C + 273.15)
  fit <- lm(log(k) ~ inv_T)
  ci <- confint(fit, level = level)[2, ]
  data.frame(Ea_kJ_mol = -coef(fit)[[2]] * R_gas / 1000,
             Ea_low = -ci[[2]] * R_gas / 1000, Ea_high = -ci[[1]] * R_gas / 1000,
             lnA = coef(fit)[[1]])
}
