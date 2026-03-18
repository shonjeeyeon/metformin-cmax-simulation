ln_params <- function(m, s) {
  s <- ifelse(s <= 0, 1e-6, s)
  v <- s^2
  phi <- sqrt(log(1 + v / m^2))
  mu <- log(m) - 0.5 * phi^2
  list(meanlog = mu, sdlog = phi)
}

simulate_pk_data <- function(pk_summary, n_sim = 1000) {
  pk_summary %>%
    group_by(group) %>%
    do({
      this <- .
      p_cmax <- ln_params(this$Cmax_mean, this$Cmax_sd)
      p_tmax <- ln_params(this$Tmax_mean, this$Tmax_sd)

      tibble(
        group   = this$group,
        regimen = this$regimen,
        dose_mg = this$dose_mg,
        eGFR    = rnorm(n_sim, this$eGFR_mean, this$eGFR_sd),
        Cmax    = rlnorm(n_sim, p_cmax$meanlog, p_cmax$sdlog),
        Tmax    = rlnorm(n_sim, p_tmax$meanlog, p_tmax$sdlog)
      )
    }) %>%
    ungroup() %>%
    mutate(
      eGFR = pmin(pmax(eGFR, 5), 130),
      Cmax_norm = Cmax / dose_mg * 850
    )
}

Nsim <- 1000
sim_data <- simulate_pk_data(pk_summary = pk_summary, n_sim = Nsim)

dir.create(here("data", "derived"), recursive = TRUE, showWarnings = FALSE)
write_csv(sim_data, here("data", "derived", "sim_data.csv"))

