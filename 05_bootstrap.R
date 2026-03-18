bootstrap_band <- function(pk_summary, grid, param = c("Cmax", "Cmax_norm", "Tmax"), B = 300, n_sim = 1000, k = 4) {
  param <- match.arg(param)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = nrow(grid))

  for (b in seq_len(B)) {
    df_b <- pk_summary %>%
      mutate(
        eGFR_mean_b = rnorm(n(), eGFR_mean, eGFR_sd / sqrt(n)),
        Cmax_mean_b = rnorm(n(), Cmax_mean, Cmax_sd / sqrt(n)),
        Tmax_mean_b = rnorm(n(), Tmax_mean, Tmax_sd / sqrt(n))
      )

    sim_b <- df_b %>%
      group_by(group) %>%
      do({
        this <- .
        p_cmax <- ln_params(this$Cmax_mean_b, this$Cmax_sd)
        p_tmax <- ln_params(this$Tmax_mean_b, this$Tmax_sd)

        tibble(
          eGFR = rnorm(n_sim, this$eGFR_mean_b, this$eGFR_sd),
          Cmax = rlnorm(n_sim, p_cmax$meanlog, p_cmax$sdlog),
          Tmax = rlnorm(n_sim, p_tmax$meanlog, p_tmax$sdlog),
          dose_mg = this$dose_mg
        )
      }) %>%
      ungroup() %>%
      mutate(
        eGFR = pmin(pmax(eGFR, 5), 130),
        Cmax_norm = Cmax / dose_mg * 850
      )

    formula <- switch(
      param,
      Cmax = Cmax ~ s(eGFR, k = k),
      Cmax_norm = Cmax_norm ~ s(eGFR, k = k),
      Tmax = Tmax ~ s(eGFR, k = k)
    )

    gam_b <- gam(formula, data = sim_b)
    boot_mat[b, ] <- predict(gam_b, newdata = grid)
  }

  tibble(
    eGFR = grid$eGFR,
    fit = apply(boot_mat, 2, mean),
    lwr = apply(boot_mat, 2, quantile, 0.025),
    upr = apply(boot_mat, 2, quantile, 0.975)
  )
}

boot_Cmax <- bootstrap_band(pk_summary, prediction_grid, param = "Cmax", B = 300, n_sim = Nsim, k = 4)
boot_Cmax_norm <- bootstrap_band(pk_summary, prediction_grid, param = "Cmax_norm", B = 300, n_sim = Nsim, k = 4)
boot_Tmax <- bootstrap_band(pk_summary, prediction_grid, param = "Tmax", B = 300, n_sim = Nsim, k = 4)

bootstrap_curves <- bind_rows(
  boot_Cmax %>% mutate(parameter = "Cmax"),
  boot_Cmax_norm %>% mutate(parameter = "Cmax_norm"),
  boot_Tmax %>% mutate(parameter = "Tmax")
)

write_csv(bootstrap_curves, here("data", "derived", "bootstrap_curves.csv"))
