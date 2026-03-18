Metformin Cmax and Tmax Simulation
================

# The Original Pharmacokinetics Table

Source:
<https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=56d13a1c-b289-4528-b23c-60f5427b4552>

# 1. Libraries and Dataset

The project now loads packages and constructs the summary dataset from
reusable scripts.

``` r
pk_summary
```

    ## # A tibble: 9 × 10
    ##   group    dose_mg regimen eGFR_mean eGFR_sd Cmax_mean Cmax_sd Tmax_mean Tmax_sd
    ##   <chr>      <dbl> <chr>       <dbl>   <dbl>     <dbl>   <dbl>     <dbl>   <dbl>
    ## 1 Healthy…     500 single      102.     14.2      1.03    0.33      2.75    0.81
    ## 2 Healthy…     850 single      102.     14.2      1.6     0.38      2.64    0.82
    ## 3 Healthy…     850 multip…     102.     14.2      2.01    0.42      1.79    0.94
    ## 4 T2D 850…     850 single       96.3    17.2      1.48    0.5       3.32    1.08
    ## 5 T2D 850…     850 multip…      96.3    17.2      1.9     0.62      2.01    1.22
    ## 6 Elderly…     850 single       78.3    14.9      2.45    0.7       2.71    1.05
    ## 7 CKD mil…     850 single       75      10        1.86    0.52      3.2     0.45
    ## 8 CKD mod…     850 single       45      10        4.12    1.83      3.75    0.5 
    ## 9 CKD sev…     850 single       20       6        3.93    0.92      4.01    1.1 
    ## # ℹ 1 more variable: n <dbl>

# 2. Simulate the Population

Simulated individual-level data are generated from the summary table
using normal assumptions for eGFR and lognormal assumptions for Cmax and
Tmax.

``` r
dplyr::glimpse(sim_data)
```

    ## Rows: 9,000
    ## Columns: 7
    ## $ group     <chr> "CKD mild (G2) 850mg", "CKD mild (G2) 850mg", "CKD mild (G2)…
    ## $ regimen   <chr> "single", "single", "single", "single", "single", "single", …
    ## $ dose_mg   <dbl> 850, 850, 850, 850, 850, 850, 850, 850, 850, 850, 850, 850, …
    ## $ eGFR      <dbl> 69.39524, 72.69823, 90.58708, 75.70508, 76.29288, 92.15065, …
    ## $ Cmax      <dbl> 1.3631176, 1.3467054, 1.7824991, 1.7275249, 0.8901151, 2.383…
    ## $ Tmax      <dbl> 2.949888, 3.275649, 2.937536, 3.758330, 3.246988, 2.907404, …
    ## $ Cmax_norm <dbl> 1.3631176, 1.3467054, 1.7824991, 1.7275249, 0.8901151, 2.383…

# 3. Fit the Parameters to Spline (`k = 4`)

Generalized additive models are used to smooth each PK endpoint over
eGFR.

``` r
dplyr::glimpse(prediction_grid)
```

    ## Rows: 111
    ## Columns: 10
    ## $ eGFR          <dbl> 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, …
    ## $ Cmax_fit      <dbl> 4.211826, 4.194010, 4.176183, 4.158337, 4.140462, 4.1225…
    ## $ Cmax_lwr      <dbl> 4.131524, 4.116568, 4.101555, 4.086472, 4.071301, 4.0560…
    ## $ Cmax_upr      <dbl> 4.292128, 4.271452, 4.250811, 4.230203, 4.209623, 4.1890…
    ## $ Tmax_fit      <dbl> 4.133506, 4.120597, 4.107683, 4.094759, 4.081822, 4.0688…
    ## $ Tmax_lwr      <dbl> 4.046930, 4.037093, 4.027202, 4.017246, 4.007214, 3.9970…
    ## $ Tmax_upr      <dbl> 4.220082, 4.204101, 4.188164, 4.172273, 4.156430, 4.1406…
    ## $ Cmax_norm_fit <dbl> 4.234914, 4.215123, 4.195322, 4.175505, 4.155664, 4.1357…
    ## $ Cmax_norm_lwr <dbl> 4.155968, 4.138987, 4.121951, 4.104850, 4.087666, 4.0703…
    ## $ Cmax_norm_upr <dbl> 4.313861, 4.291259, 4.268693, 4.246161, 4.223661, 4.2011…

# 4. Build Plots for Cmax and Tmax

<img src="metformin_cmax_tmax_simulation_files/figure-gfm/joint-plot-1.png" width="80%" style="display: block; margin: auto auto auto 0;" />

# 5. Bootstrapped Confidence Bands

``` r
head(bootstrap_curves)
```

    ## # A tibble: 6 × 5
    ##    eGFR   fit   lwr   upr parameter
    ##   <dbl> <dbl> <dbl> <dbl> <chr>    
    ## 1    10  4.20  3.26  5.09 Cmax     
    ## 2    11  4.18  3.24  5.05 Cmax     
    ## 3    12  4.17  3.24  5.01 Cmax     
    ## 4    13  4.15  3.26  4.96 Cmax     
    ## 5    14  4.13  3.28  4.93 Cmax     
    ## 6    15  4.11  3.27  4.89 Cmax

# 6. Plot the Confidence Bands

<img src="metformin_cmax_tmax_simulation_files/figure-gfm/boot-cmax-1.png" width="80%" style="display: block; margin: auto auto auto 0;" />

<img src="metformin_cmax_tmax_simulation_files/figure-gfm/boot-cmax-norm-1.png" width="80%" style="display: block; margin: auto auto auto 0;" />

<img src="metformin_cmax_tmax_simulation_files/figure-gfm/boot-tmax-1.png" width="80%" style="display: block; margin: auto auto auto 0;" />

# 7. Compare the Three Curves

<img src="metformin_cmax_tmax_simulation_files/figure-gfm/combined-boot-1.png" width="80%" style="display: block; margin: auto auto auto 0;" />
