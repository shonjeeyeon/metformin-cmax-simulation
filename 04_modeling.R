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

  grid %>%
    mutate(
      Cmax_fit = pred_Cmax$fit,
      Cmax_lwr = Cmax_fit - 1.96 * pred_Cmax$se.fit,
      Cmax_upr = Cmax_fit + 1.96 * pred_Cmax$se.fit,
      Tmax_fit = pred_Tmax$fit,
      Tmax_lwr = Tmax_fit - 1.96 * pred_Tmax$se.fit,
      Tmax_upr = Tmax_fit + 1.96 * pred_Tmax$se.fit,
      Cmax_norm_fit = pred_Cmax_norm$fit,
      Cmax_norm_lwr = Cmax_norm_fit - 1.96 * pred_Cmax_norm$se.fit,
      Cmax_norm_upr = Cmax_norm_fit + 1.96 * pred_Cmax_norm$se.fit
    )
}

models <- fit_pk_models(sim_data = sim_data, k = 4)
prediction_grid <- make_prediction_grid(models = models)

write_csv(prediction_grid, here("data", "derived", "prediction_grid.csv"))
