fit_pk_models <- function(sim_data, k = 4) {
  list(
    gam_Cmax = gam(Cmax ~ s(eGFR, k = k), data = sim_data),
    gam_Tmax = gam(Tmax ~ s(eGFR, k = k), data = sim_data),
    gam_Cmax_norm = gam(Cmax_norm ~ s(eGFR, k = k), data = sim_data)
  )
}

make_prediction_grid <- function(models, egfr_seq = seq(10, 120, by = 1)) {
  grid <- tibble(eGFR = egfr_seq)

  pred_Cmax <- predict(models$gam_Cmax, grid, se.fit = TRUE)
  pred_Tmax <- predict(models$gam_Tmax, grid, se.fit = TRUE)
  pred_Cmax_norm <- predict(models$gam_Cmax_norm, grid, se.fit = TRUE)

  cmax_fit <- as.numeric(pred_Cmax$fit)
  cmax_se <- as.numeric(pred_Cmax$se.fit)
  tmax_fit <- as.numeric(pred_Tmax$fit)
  tmax_se <- as.numeric(pred_Tmax$se.fit)
  cmax_norm_fit <- as.numeric(pred_Cmax_norm$fit)
  cmax_norm_se <- as.numeric(pred_Cmax_norm$se.fit)

  grid %>%
    mutate(
      Cmax_fit = cmax_fit,
      Cmax_lwr = cmax_fit - 1.96 * cmax_se,
      Cmax_upr = cmax_fit + 1.96 * cmax_se,
      Tmax_fit = tmax_fit,
      Tmax_lwr = tmax_fit - 1.96 * tmax_se,
      Tmax_upr = tmax_fit + 1.96 * tmax_se,
      Cmax_norm_fit = cmax_norm_fit,
      Cmax_norm_lwr = cmax_norm_fit - 1.96 * cmax_norm_se,
      Cmax_norm_upr = cmax_norm_fit + 1.96 * cmax_norm_se
    )
}

models <- fit_pk_models(sim_data = sim_data, k = 4)
prediction_grid <- make_prediction_grid(models = models)

dir.create(here("data", "derived"), recursive = TRUE, showWarnings = FALSE)
write_csv(prediction_grid, here("data", "derived", "prediction_grid.csv"))
