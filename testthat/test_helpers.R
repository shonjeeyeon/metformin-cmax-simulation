source(here::here("R", "01_packages.R"))
source(here::here("R", "03_simulation.R"))

stopifnot(is.list(ln_params(1, 0.2)))
stopifnot(all(c("meanlog", "sdlog") %in% names(ln_params(1, 0.2))))
