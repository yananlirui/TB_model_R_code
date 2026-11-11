################################################################################
# Code accompanying Li et al., "Evaluating the impact of expanding tuberculosis
# diagnostic and treatment availability of general hospitals in China".
#
# Purpose: Project the health and economic effects of expanding general-hospital
# TB services in China during 2026-2035, using a dynamic transmission model
# with explicit care-seeking, referral, treatment and recovery pathways.
#
# Scope: Pulmonary TB across all ages; extrapulmonary TB is not modelled.
# Incidence counts disease episodes; notifications count flows into treatment.
# The workflow includes five calibrations and three complete analyses.
# Ordinary differential equations are solved using deSolve.
#
# Run: Rscript TB_model_main.R, or source the saved script in a fresh R session.
# Inputs are embedded below. Results, tables and figures are saved in output/.
# The latest saved RDS for the selected seed is reused when results exist.
# A fresh calibration can take many hours and uses parallel worker processes.
# See README.md for settings, dependencies and instructions for a fresh run.
#
# S0. Environment, run settings, parameters and optional resistance verification
# S1. Population, model equations, calibration, CEA and sensitivity functions
# S2. Calibration and analysis execution, diagnostics and result checkpoints
# S3. Excel tables, calibration diagnostics and publication figures
################################################################################
# ==============================================================================
# S0. Computational environment, model parameters, and analysis settings
# ==============================================================================
# ------------------------------------------------------------------------------
# 0.1 Computational environment and required packages
# ------------------------------------------------------------------------------
# Verified environment: R 4.6.1 (2026-06-24).
rm(list = ls())           # Remove all objects from the working environment
start_time <- Sys.time()  # Record script start time for performance monitoring

# Install and load required packages
# Package versions:
#   coda 0.19.4.1 | cowplot 1.2.0   | deSolve 1.42      | dplyr 1.2.1
#   FME 1.3.6.4   | ggplot2 4.0.3   | grid 4.6.1        | gridExtra 2.3.1
#   gtools 3.9.5  | lubridate 1.9.5 | parallel 4.6.1    | patchwork 1.3.2
#   pbapply 1.7.5 | pracma 2.4.6    | rstudioapi 0.19.0 | triangle 1.1.0
#   writexl 2.0.1
install_and_load <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = NA)
    library(package_name, character.only = TRUE)
  }
}

# Install and load required packages (sorted alphabetically)
required_packages <- c(
  "coda",        # MCMC diagnostics and posterior sample manipulation
  "cowplot",     # Combine and arrange ggplot figures
  "deSolve",     # Solve ordinary differential equations
  "dplyr",       # Data manipulation
  "FME",         # Inverse modelling, sensitivity, and Monte Carlo analysis
  "ggplot2",     # Publication-quality graphics
  "grid",        # Grid graphics and layouts
  "gridExtra",   # Additional grid plotting utilities
  "gtools",      # Miscellaneous R utilities
  "lubridate",   # Date and time manipulation
  "parallel",    # Parallel computing
  "patchwork",   # Combine ggplot objects
  "pbapply",     # Progress bars for apply functions
  "pracma",      # Numerical methods and integration
  "rstudioapi",  # Locate scripts in RStudio
  "triangle",    # Triangular distribution sampling
  "writexl"      # Export Excel workbooks
)
invisible(lapply(required_packages, install_and_load))

# To verify the computational environment on the current machine,
# change verify_environment to TRUE and run the block below.
verify_environment <- FALSE
if (verify_environment) {
  cat(R.version.string, "\n\n")
  print(vapply(
    required_packages,
    function(pkg) as.character(utils::packageVersion(pkg)),
    character(1)
  ))
}

# ------------------------------------------------------------------------------
# 0.2 Locate this script and define data and output directories
# ------------------------------------------------------------------------------

get_script_dir <- function() {
  # Run via source().
  for (frame in rev(sys.frames())) {
    if (is.character(frame$ofile) && length(frame$ofile) == 1L &&
      !is.na(frame$ofile) && nzchar(frame$ofile)) {
      return(dirname(normalizePath(
        frame$ofile, winslash = "/", mustWork = TRUE
      )))
    }
  }

  # Run via Rscript.
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) > 0L) {
    return(dirname(normalizePath(
      sub("^--file=", "", script_arg[1]),
      winslash = "/", mustWork = TRUE
    )))
  }

  # Run directly in RStudio.
  if (rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getSourceEditorContext()$path
    if (is.character(script_path) && length(script_path) == 1L &&
      !is.na(script_path) && nzchar(script_path)) {
      return(dirname(normalizePath(
        script_path, winslash = "/", mustWork = TRUE
      )))
    }
  }

  stop("Save this script and run it using RStudio, source(), or Rscript.")
}

file_dir <- get_script_dir()

# Save figures and tables in the output folder alongside this script.
save_dir <- file.path(file_dir, "output")
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(save_dir)) {
  stop("Cannot create the output directory: ", save_dir)
}

# ------------------------------------------------------------------------------
# 0.3 Analysis settings and result loading
# ------------------------------------------------------------------------------
# Set the seed, study periods and calibration controls.
settings <- list(
  seed       = 2026L,
  years      = 2000:2035,
  pred_years = 2015:2035,
  eval_years = 2026:2035,
  calib      = list(nchains = 4L, niter = 250000L)
)
igra_scenarios <- c(0, 0.5, 1)

set.seed(settings$seed)

# Read the latest results for this seed; exclude compiled backup files.
results_file <- list.files(
  save_dir, paste0("^TB_results_seed", settings$seed, "_[0-9]{8}[.]rds$"),
  full.names = TRUE
)
results_file <- if (length(results_file)) {
  tail(sort(results_file), 1L)
} else ""
model_results <- if (nzchar(results_file)) readRDS(results_file) else list()
calibration_results <- model_results$calibration$results

# Rerun everything unless all four results exist for this seed.
run_all <- !identical(as.integer(model_results$seed), settings$seed) ||
  is.null(calibration_results[["0"]]) ||
  is.null(calibration_results[["0.5"]]) ||
  is.null(calibration_results[["1"]]) || is.null(model_results$analysis)
tb_analysis_results <- if (run_all) NULL else model_results$analysis
if (run_all) message("A full run is required and may take many hours.")

# ------------------------------------------------------------------------------
# 0.4 Model parameters
# ------------------------------------------------------------------------------
# ---- 0.4.1 Base-case parameter values ----
base_par <- list(
  # Natural history and treatment; rates per year
  epsilon    = 1 / 1.2,    # Relative infectiousness: SA versus CA
  c_f        = 0.0000014,  # Fast-latent infection self-clearance
  c_s        = 0.033,      # Slow-latent infection self-clearance
  delta1     = 0.0859,     # Transition rate: FL to SA
  delta3     = 0.000622,   # Reactivation rate: SL to SA
  delta4     = 0.69,       # Transition rate: SA to CA
  delta5     = 0.58,       # Transition rate: CA to SA
  gamma_sen  = 2,          # Treatment transition rate, RS-TB
  gamma_res  = 0.5,        # Treatment transition rate, RR-TB
  e1         = 0.75,       # Reinfection protection in SL and R
  e2         = 0.5,        # E protection relative to SL and R
  p_igra_E   = 0,          # IGRA positivity after clearance
  or_sp      = 0.106,      # Smear-positive OR: SA relative to CA
  tao_sp     = 0.231,      # Natural cure of untreated SP-TB
  tao_sn     = 0.130,      # Natural cure of untreated SN-TB
  rho_low    = 0.032,      # Short-term relapse: R_low to SA
  rho_high   = 0.14,       # Short-term relapse: R_high to SA
  rho_long   = 0.0015,     # Long-term relapse: R_long to SA
  Kappa      = 0.5,        # R_low/R_high to R_long transition
  mu_sp      = 0.389,      # Untreated smear-positive TB mortality
  mu_sn      = 0.025,      # Untreated smear-negative TB mortality
  mu_ta      = 0.01449,    # TB mortality during treatment

  # Care-seeking pathway
  p_pc       = 0.385,      # Initial care-seeking at primary care
  p_gh       = 0.534,      # Initial care-seeking at GH
  p_ref      = 0.743,      # PA referral to GH rather than TA1
  p_lf1      = 1 - 0.132,  # Loss to follow-up during PA referral
  p_lf2      = 1 - 0.25,   # Loss to follow-up during GA referral

  # CEA parameters
  cost_dx_rs = 486.4,      # RS-TB diagnosis per patient
  cost_dx_rr = 1072.0,     # RR-TB diagnosis per patient
  cost_tx_rs = 1676.1,     # RS-TB treatment per patient
  cost_tx_rr = 26586.9,    # RR-TB treatment per patient
  mult_dx    = 3.779,      # Evaluated per diagnosed PTB
  mult_tx    = 1.107,      # Treated per true PTB case
  cost_equip = 571685.8,   # Equipment per hospital
  cost_train = 106430,     # Training per hospital
  prog_share = 0.5,        # Programme share of setup
  equip_life = 10,         # Equipment life, years
  dw_tb      = 0.333,      # Active-TB disability weight
  yll_tb     = 26.05,      # YLLs per TB death
  disc_rate  = 0.03        # Annual discount rate
)

# ---- 0.4.2 Parameter ranges for one-way sensitivity analysis ----
# Reported 95% intervals, plausible bounds, or assumed scenario ranges.
sen_prop <- c(0.8, 1.2)

par_ranges <- list(
  # TB natural history and treatment
  epsilon    = c(1 / 2.3,    1),
  c_f        = c(0.00000029, 0.00000261),
  c_s        = c(0.0254,     0.0467),
  delta1     = c(0.0696,     0.111),
  delta3     = c(0.000135,   0.00113),
  delta4     = c(0.54,       0.88),
  delta5     = c(0.45,       0.72),
  gamma_sen  = sen_prop * base_par$gamma_sen,
  gamma_res  = sen_prop * base_par$gamma_res,
  e1         = c(0.60,       0.85),
  e2         = sen_prop * base_par$e2,
  p_igra_E   = c(0,          1),
  or_sp      = c(0.0369,     0.2717),
  tao_sp     = c(0.177,      0.288),
  tao_sn     = c(0.073,      0.209),
  rho_low    = c(0.016,      0.048),
  rho_high   = c(0.07,       0.21),
  rho_long   = sen_prop * base_par$rho_long,
  Kappa      = sen_prop * base_par$Kappa,
  mu_sp      = c(0.335,      0.449),
  mu_sn      = c(0.017,      0.035),
  mu_ta      = c(0.00947,    0.02816),

  # Care-seeking pathway
  p_pc       = qbeta(c(0.025, 0.975), 114, 182),
  p_gh       = qbeta(c(0.025, 0.975), 158, 138),
  p_ref      = qbeta(c(0.025, 0.975), 85,  30),
  p_lf1      = c(0.8 * base_par$p_lf1, 1),
  p_lf2      = sen_prop * base_par$p_lf2,

  # CEA parameters
  cost_dx_rs = sen_prop * base_par$cost_dx_rs,
  cost_dx_rr = sen_prop * base_par$cost_dx_rr,
  cost_tx_rs = sen_prop * base_par$cost_tx_rs,
  cost_tx_rr = sen_prop * base_par$cost_tx_rr,
  mult_dx    = 1 / rev(qbeta(c(0.025, 0.975), 79913, 222095)),
  mult_tx    = 1 / rev(qbeta(c(0.025, 0.975), 1592,  171)),
  cost_equip = c(521525.8,   621845.9),
  cost_train = c(24650,      188210),
  prog_share = c(0.073,      0.76),
  equip_life = c(5,          20),
  dw_tb      = c(0.226,      0.450),
  yll_tb     = c(16.5,       41.2),
  disc_rate  = c(0,          0.06)
)

# ---- 0.4.3 Parameter sampling for probabilistic sensitivity analysis ----
sample_known_par <- function(base_par, par_ranges,
                             n = 10000L, seed = NULL) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n < 1 || n != floor(n) || n > .Machine$integer.max) {
    stop("n must be a positive integer within R's integer range.")
  }
  n <- as.integer(n)
  if (!is.null(seed)) set.seed(seed)

  tri_sample <- function(n, mode, range) {
    if (range[1] == range[2]) {
      return(rep(range[1], n))
    }
    triangle::rtriangle(n, a = range[1], b = range[2], c = mode)
  }

  # Joint sampling preserves dependence among PA, GH, and direct care.
  multinomial_samples <- gtools::rdirichlet(n, c(114, 158, 24))

  samples <- list(
    # TB natural history and treatment
    epsilon    = rbeta(n, 3.99, 0.80),
    c_f        = tri_sample(n, base_par$c_f,        par_ranges$c_f),
    c_s        = tri_sample(n, base_par$c_s,        par_ranges$c_s),
    delta1     = tri_sample(n, base_par$delta1,     par_ranges$delta1),
    delta3     = tri_sample(n, base_par$delta3,     par_ranges$delta3),
    delta4     = rlnorm(n, log(base_par$delta4), 0.125),
    delta5     = rlnorm(n, log(base_par$delta5), 0.120),
    gamma_sen  = tri_sample(n, base_par$gamma_sen,  par_ranges$gamma_sen),
    gamma_res  = tri_sample(n, base_par$gamma_res,  par_ranges$gamma_res),
    e1         = rbeta(n, 27.47, 9.16),
    e2         = tri_sample(n, base_par$e2,         par_ranges$e2),
    p_igra_E   = runif(n, par_ranges$p_igra_E[1], par_ranges$p_igra_E[2]),
    or_sp      = rlnorm(n, log(base_par$or_sp), 0.456),
    tao_sp     = rlnorm(n, log(base_par$tao_sp), 0.124),
    tao_sn     = rlnorm(n, log(base_par$tao_sn), 0.269),
    rho_low    = tri_sample(n, base_par$rho_low,    par_ranges$rho_low),
    rho_high   = tri_sample(n, base_par$rho_high,   par_ranges$rho_high),
    rho_long   = tri_sample(n, base_par$rho_long,   par_ranges$rho_long),
    Kappa      = tri_sample(n, base_par$Kappa,      par_ranges$Kappa),
    mu_sp      = rlnorm(n, log(base_par$mu_sp), 0.0747),
    mu_sn      = rlnorm(n, log(base_par$mu_sn), 0.1842),
    mu_ta      = tri_sample(n, base_par$mu_ta,      par_ranges$mu_ta),

    # Care-seeking pathway
    p_pc       = multinomial_samples[, 1],
    p_gh       = multinomial_samples[, 2],
    p_ref      = rbeta(n, 85, 30),
    p_lf1      = tri_sample(n, base_par$p_lf1,      par_ranges$p_lf1),
    p_lf2      = tri_sample(n, base_par$p_lf2,      par_ranges$p_lf2),

    # CEA parameters
    cost_dx_rs = tri_sample(n, base_par$cost_dx_rs, par_ranges$cost_dx_rs),
    cost_dx_rr = tri_sample(n, base_par$cost_dx_rr, par_ranges$cost_dx_rr),
    cost_tx_rs = tri_sample(n, base_par$cost_tx_rs, par_ranges$cost_tx_rs),
    cost_tx_rr = tri_sample(n, base_par$cost_tx_rr, par_ranges$cost_tx_rr),
    mult_dx    = 1 / rbeta(n, 79913, 222095),
    mult_tx    = 1 / rbeta(n, 1592, 171),
    cost_equip = tri_sample(n, base_par$cost_equip, par_ranges$cost_equip),
    cost_train = tri_sample(n, base_par$cost_train, par_ranges$cost_train),
    prog_share = tri_sample(n, base_par$prog_share, par_ranges$prog_share),
    equip_life = tri_sample(n, base_par$equip_life, par_ranges$equip_life),
    dw_tb      = rbeta(n, 22.2, 44.4),
    yll_tb     = rlnorm(n, 3.259954, 0.2334108),
    disc_rate  = runif(n, min = 0, max = 0.06)
  )

  as.data.frame(samples)
}

# ---- 0.4.4 Rifampicin resistance trend: Supplementary Figure 2 ----
# Three-parameter exponential trend fitted with R nlminb.
resistance_prop <- function(year) {
  3.722206 * exp(-0.000594817 * (year - 2000)) - 3.583639
}

# Set verify_resistance to TRUE to refit and redraw Supplementary Fig. 2.
verify_resistance <- FALSE
if (verify_resistance) {
  # Study start year, end year and RR-TB proportion; Supplementary Table 4.
  resistance_data <- as.data.frame(rbind(
    c(2000,   2001,   0.1658), c(2010,   2011,   0.0893),
    c(2015.7, 2016.7, 0.1142), c(2016,   2017,   0.0839),
    c(2017,   2018,   0.0625), c(2018,   2019,   0.0795),
    c(2019,   2020,   0.0593), c(2020,   2021,   0.0422),
    c(2021,   2022,   0.0458), c(2022,   2023,   0.0410),
    c(2004,   2005,   0.0560), c(2005.2, 2007.5, 0.0617),
    c(2010,   2011.5, 0.0857), c(2016,   2019,   0.0593),
    c(2017,   2020,   0.0566), c(2004.7, 2005.7, 0.1452),
    c(2014.7, 2017.7, 0.1010), c(2015,   2016,   0.1031),
    c(2016,   2017,   0.2414), c(2017,   2018,   0.2831),
    c(2018,   2019,   0.1662), c(2019,   2020,   0.1516),
    c(2020,   2021,   0.1342), c(2021,   2022,   0.1037),
    c(2022,   2023,   0.1125), c(2000,   2005,   0.1055),
    c(2006.7, 2010.4, 0.1395), c(2009,   2010.2, 0.0659),
    c(2011,   2015,   0.0807), c(2012,   2013,   0.0627),
    c(2013,   2014,   0.0782), c(2014,   2015,   0.0675),
    c(2015,   2016,   0.0678), c(2016,   2017,   0.0881),
    c(2017,   2018,   0.1013), c(2018,   2019,   0.0927),
    c(2019,   2020,   0.0866), c(2020,   2021,   0.0754),
    c(2016,   2017,   0.0478), c(2017,   2018,   0.0428),
    c(2018,   2019,   0.0799), c(2019,   2020,   0.0479),
    c(2020,   2021,   0.0981), c(2010.5, 2011.5, 0.0788),
    c(2012,   2016,   0.0831), c(2013.8, 2014.2, 0.1060),
    c(2015,   2016,   0.0948), c(2014,   2015,   0.1776),
    c(2020,   2021,   0.0849), c(2004,   2005,   0.1694),
    c(2017,   2021,   0.1569), c(2001.5, 2002.5, 0.1547),
    c(2010.2, 2013.2, 0.1081), c(2012.6, 2015.7, 0.1495),
    c(2014,   2015,   0.0758), c(2015,   2016,   0.0916),
    c(2016,   2017,   0.0858), c(2017,   2018,   0.0802),
    c(2018,   2019,   0.0354), c(2008.3, 2009,   0.1941),
    c(2009.8, 2014.5, 0.1379), c(2011,   2011.8, 0.1503),
    c(2013,   2014,   0.0530), c(2019,   2020,   0.0580),
    c(2020,   2021,   0.0530), c(2008.8, 2011.4, 0.1676),
    c(2012,   2012.5, 0.1886), c(2013,   2014,   0.2420),
    c(2014,   2015,   0.1857), c(2015,   2016,   0.1489),
    c(2016,   2017,   0.1887), c(2017,   2018,   0.1633),
    c(2018,   2019,   0.1053), c(2019,   2020,   0.1040),
    c(2013.2, 2014.3, 0.0992), c(2013,   2015,   0.0424),
    c(2018,   2019,   0.1818), c(2009.2, 2010.5, 0.0778),
    c(2017,   2019,   0.0984), c(2004.8, 2012,   0.0437),
    c(2004.1, 2005.1, 0.0695), c(2011,   2014,   0.1328),
    c(2010.5, 2013,   0.1066), c(2008,   2009.5, 0.0793),
    c(2006.4, 2007,   0.5161), c(2004.2, 2005.4, 0.0803),
    c(2009,   2012,   0.0871), c(2001,   2002,   0.1065),
    c(2002,   2003,   0.1165), c(2003,   2004,   0.1071),
    c(2004,   2005,   0.0974), c(2004.3, 2005,   0.0954),
    c(2008.5, 2009.6, 0.0800), c(2013.5, 2014.5, 0.0641),
    c(2018,   2019,   0.0347)
  ))
  colnames(resistance_data) <- c("year_start", "year_end", "proportion")

  # Equal-weight fit at study midpoints.
  resistance_data$year <- with(
    resistance_data, (year_start + year_end) / 2 - 2000
  )
  resistance_loss <- function(par) {
    fitted <- par[1] * exp(par[2] * resistance_data$year) + par[3]
    if (any(!is.finite(fitted)) || any(abs(fitted) > 1e100)) {
      return(1e200)
    }
    sum((fitted - resistance_data$proportion)^2)
  }
  resistance_fit <- nlminb(
    start     = c(a = 3, b = -0.1, c = -4),
    objective = resistance_loss,
    control   = list(iter.max = 5000, eval.max = 10000)
  )
  if (resistance_fit$convergence != 0) stop(resistance_fit$message)
  resistance_coef <- resistance_fit$par

  # Reproduce the fitted curve and compare with the formula above.
  resistance_curve <- data.frame(year = seq(2000, 2035, by = 0.1))
  resistance_curve$fitted <- with(resistance_curve,
    resistance_coef[1] * exp(resistance_coef[2] * (year - 2000)) +
      resistance_coef[3]
  )
  cat(sprintf(
    "a = %.6f, b = %.9f, c = %.6f\n",
    resistance_coef[1], resistance_coef[2], resistance_coef[3]
  ))
  cat("Maximum difference from stored formula:", max(abs(
    resistance_curve$fitted - resistance_prop(resistance_curve$year)
  )), "\n")
  # Show the original figure range; keep all records in the fit above.
  resistance_plot <- ggplot(resistance_data) +
    geom_segment(aes(
      x = year_start, xend = year_end,
      y = proportion * 100, yend = proportion * 100
    ), colour = "grey45", linewidth = 0.25) +
    geom_point(
      aes(x = year_start, y = proportion * 100),
      colour = "grey25", size = 0.35
    ) +
    geom_point(
      aes(x = year_end, y = proportion * 100),
      colour = "#A96464", size = 0.35
    ) +
    geom_point(aes(
      x = (year_start + year_end) / 2, y = proportion * 100
    ), shape = 1, colour = "#B68B8B", size = 1.4, stroke = 0.25) +
    geom_line(
      data = resistance_curve, aes(x = year, y = fitted * 100),
      colour = "#A8DDA0", linewidth = 0.5
    ) +
    scale_x_continuous(breaks = seq(2000, 2035, by = 5), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 30, by = 5), expand = c(0, 0)) +
    coord_cartesian(xlim = c(2000, 2035), ylim = c(0, 30)) +
    labs(x = "Year", y = "Proportion (%)") +
    theme_classic(base_size = 12, base_family = "Times New Roman") +
    theme(
      axis.line          = element_line(colour = "grey55", linewidth = 0.3),
      axis.ticks         = element_line(colour = "grey55", linewidth = 0.3),
      axis.ticks.length  = grid::unit(-2, "mm"),
      axis.text          = element_text(colour = "black"),
      axis.title         = element_text(size = 13),
      panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3),
      plot.margin        = margin(8, 10, 5, 5)
    )

  ggsave(
    file.path(save_dir, "Supplementary_Fig2_rifampicin_resistance.png"),
    resistance_plot, width = 10, height = 5.25, dpi = 300
  )
}

# ==============================================================================
# S1. Model and analysis functions
# ==============================================================================
# ------------------------------------------------------------------------------
# 1.1 Initialize the model population
# ------------------------------------------------------------------------------

initialize_population <- function(params) {
  # Read initial-condition parameters by name.
  cp         <- as.list(params)
  pre_e      <- cp$Pre_E
  pre_ltbi   <- cp$Pre_LTBI
  prop_fl    <- cp$Prop_FL
  pre_atb    <- cp$Pre_ATB
  prop_sa    <- cp$Prop_SA
  prop_pa    <- cp$Prop_PA
  prop_ga    <- cp$Prop_GA
  prop_ta    <- cp$Prop_TA
  population <- 100000

  # Check that the initial compartment allocation is valid.
  proportions <- c(
    pre_e, pre_ltbi, prop_fl, pre_atb,
    prop_sa, prop_pa, prop_ga, prop_ta
  )

  if (any(!is.finite(proportions)) ||
        any(proportions < 0 | proportions > 1) ||
        pre_ltbi + pre_e + pre_atb > 1 ||
        prop_sa + prop_pa + prop_ga + prop_ta > 1) {
    stop("Invalid initial population proportions.")
  }

  y <- c(
    S  = population * (1 - pre_ltbi - pre_e - pre_atb),
    E  = population * pre_e,
    FL = population * pre_ltbi * prop_fl,
    SL = population * pre_ltbi * (1 - prop_fl),
    SA = population * pre_atb * prop_sa,
    CA = population * pre_atb *
      (1 - prop_sa - prop_pa - prop_ga - prop_ta),
    PA     = population * pre_atb * prop_pa,
    GA     = population * pre_atb * prop_ga,
    TA1    = population * pre_atb * prop_ta,
    TA2    = 0,
    R_low  = 0,
    R_high = 0,
    R_long = 0,
    D      = 0
  )

  y
}

# ------------------------------------------------------------------------------
# 1.2 TB transmission model and model outcomes
# ------------------------------------------------------------------------------

# Prepare fixed inputs once; calibration and analysis share the same equations.
tb_model <- function(calib_par, base_par, scenarios = 0,
                     assumptions = "basecase", calibration_only = FALSE) {
  # Demography and treatment success, annual observations for 2000-2023.
  demographic_values <- list(
    birth   = c(
      14.03, 13.38, 12.86, 12.41, 12.29, 12.40, 12.09, 12.10,
      12.14, 11.95, 11.90, 13.27, 14.57, 13.03, 13.83, 11.99,
      13.57, 12.64, 10.86, 10.41, 8.52,  7.52,  6.77,  6.39
    ),
    death   = c(
      6.45,  6.43,  6.41,  6.40,  6.42,  6.51,  6.81,  6.93,
      7.06,  7.08,  7.11,  7.13,  7.13,  7.13,  7.12,  7.07,
      7.04,  7.06,  7.08,  7.09,  7.07,  7.18,  7.37,  7.87
    ),
    success = c(
      93,    95,    92,    93,    92,    92,    92,    93,
      93,    94,    95,    95,    95,    95,    94,    94,
      93,    93,    94,    94,    95,    94,    95,    95
    )
  )
  demographic_inputs <- lapply(demographic_values, function(values) {
    stats::approxfun(2000:2023, values, rule = 2)
  })
  if (!assumptions %in% c("basecase", "worst", "best")) {
    stop("Assumptions must be 'basecase', 'worst', or 'best'.")
  }
  if (length(scenarios) != 1L || !scenarios %in% 0:3) {
    stop("Scenarios must be 0, 1, 2, or 3.")
  }
  if (assumptions == "worst") {
    cum_hosp  <- cumsum(c(9.8, 23.7, 46.6, 19.8))
    cum_visit <- cumsum(c(64.7, 23.4, 8.3, 3.5))
  } else if (assumptions == "best") {
    cum_hosp  <- cumsum(c(19.8, 46.6, 23.7, 9.8))
    cum_visit <- cumsum(c(3.5, 8.3, 23.4, 64.7))
  }
  visit_map <- if (assumptions %in% c("worst", "best")) {
    stats::approxfun(c(0, cum_hosp), c(0, cum_visit), rule = 2)
  } else {
    NULL
  }
  # Time-varying demographics and hospital coverage.
  hospital_coverage <- function(t) {
    if (t < 2010) {
      result <- c(0, 0)
    } else if (t < 2015) {
      result <- c(
        min(51 / 13681 + 0.018588921 * (t - 2010), 2921 / 20133),
        min(51 / 2921 + 0.111879493 * (t - 2010), 1)
      )
    } else {
      result <- c(
        min(1685 / 17430 + 0.002877808 * (t - 2015), 2921 / 20133),
        min(1685 / 2921 + 0.037726806 * (t - 2015), 1)
      )
    }
    if (scenarios %in% 1:3 && t > 2025) {
      start_coverage     <- min(
        1685 / 17430 + 0.002877808 * (2025 - 2015), 2921 / 20133
      )
      target_coverage    <- c(0.5, 0.8, 1)[scenarios]
      expansion_progress <- pmin(1, (t - 2025) / (2030 - 2025))
      result[1] <- start_coverage +
        expansion_progress * (target_coverage - start_coverage)
    }
    hosp_cover <- result[1]
    if (!is.null(visit_map)) {
      result[1] <- visit_map(result[1] * 100) / 100
    }
    c(result, hosp_cover)
  }
  # Prepare fixed natural-history, care and cost parameters.
  epsilon    <- base_par$epsilon
  c_f        <- base_par$c_f
  c_s        <- base_par$c_s
  delta1     <- base_par$delta1
  delta3     <- base_par$delta3
  delta4     <- base_par$delta4
  delta5     <- base_par$delta5
  gamma_sen  <- base_par$gamma_sen
  gamma_res  <- base_par$gamma_res
  e1         <- base_par$e1
  e2         <- base_par$e2
  p_igra_e   <- base_par$p_igra_E
  or_sp      <- base_par$or_sp
  tao_sp     <- base_par$tao_sp
  tao_sn     <- base_par$tao_sn
  rho_low    <- base_par$rho_low
  rho_high   <- base_par$rho_high
  rho_long   <- base_par$rho_long
  kappa      <- base_par$Kappa
  mu_sp      <- base_par$mu_sp
  mu_sn      <- base_par$mu_sn
  mu_ta      <- base_par$mu_ta
  p_pc       <- base_par$p_pc
  p_gh       <- base_par$p_gh
  p_ref      <- base_par$p_ref
  p_lf1      <- base_par$p_lf1
  p_lf2      <- base_par$p_lf2
  cost_dx_rs <- base_par$cost_dx_rs
  cost_dx_rr <- base_par$cost_dx_rr
  cost_tx_rs <- base_par$cost_tx_rs
  cost_tx_rr <- base_par$cost_tx_rr
  mult_dx    <- base_par$mult_dx
  mult_tx    <- base_par$mult_tx
  cp         <- as.list(calib_par)
  beta       <- cp$Beta1
  delta2     <- cp$delta2
  theta_2023 <- cp$theta_2023
  eta        <- cp$eta
  omiga1     <- cp$omiga1
  omiga2     <- cp$omiga2
  ur         <- cp$Ur
  a          <- cp$a
  b          <- cp$b
  c          <- cp$c
  # Shared state equations; only the returned outputs differ by mode.
  function(t, y, parms) {
    # Read the current compartment sizes.
    state  <- unname(y)
    s      <- state[1]
    e      <- state[2]
    fl     <- state[3]
    sl     <- state[4]
    sa     <- state[5]
    ca     <- state[6]
    pa     <- state[7]
    ga     <- state[8]
    ta1    <- state[9]
    ta2    <- state[10]
    r_low  <- state[11]
    r_high <- state[12]
    r_long <- state[13]
    r      <- r_low + r_high + r_long
    n_pop  <- sum(state[1:13])
    ca_all <- ca + pa + ga + ta1 + ta2
    atb    <- sa + ca_all
    # Time-varying demographics, resistance and hospital coverage.
    resistance_share <- resistance_prop(t)
    smear_share      <- a * exp(-b * (t - 2000)) + c
    b_t              <- demographic_inputs$birth(t) / 1000
    mu_t             <- demographic_inputs$death(t) / 1000
    alpha_t          <- demographic_inputs$success(t) / 100
    coverage         <- hospital_coverage(t)
    gh_cover         <- coverage[1]
    hosp_cover       <- coverage[3]
    gh_share         <- coverage[2]
    gamma            <- 1 / (
      (1 - resistance_share) / gamma_sen + resistance_share / gamma_res
    )
    # Smear positivity and state-specific natural history.
    sa_share <- sa / max(atb, .Machine$double.eps)
    d        <- 1 - or_sp
    smear_a  <- 1 - sa_share * d
    smear_b  <- (1 - sa_share) * d
    smear_h  <- smear_a + smear_share * d
    smear_ca <- 2 * smear_share /
      (smear_h + sqrt(smear_h^2 - 4 * smear_b * smear_share))
    smear_sa <- or_sp * smear_ca / (1 - smear_ca + or_sp * smear_ca)
    mu_sa    <- mu_sp * smear_sa + mu_sn * (1 - smear_sa)
    mu_ca    <- mu_sp * smear_ca + mu_sn * (1 - smear_ca)
    tao_sa   <- tao_sp * smear_sa + tao_sn * (1 - smear_sa)
    tao_ca   <- tao_sp * smear_ca + tao_sn * (1 - smear_ca)
    # Transmission, recurrence and notification flows.
    lambda        <- beta * (epsilon * sa + ca + pa + ga) / n_pop
    reinf_rate    <- (1 - e1) * lambda * delta1 / (delta1 + delta2)
    rel_flow      <- rho_low * r_low + rho_high * r_high + rho_long * r_long
    reinf_flow    <- reinf_rate * r
    notif_passive <- eta * (1 - p_gh - p_pc + p_gh * gh_cover) * ca +
      omiga1 * (1 - p_lf1) * (1 - p_ref + p_ref * gh_cover) * pa +
      omiga2 * (1 - p_lf2) * ga
    theta_t     <- theta_2023 * pmin(1, pmax(0, (t - 2000) / (2023 - 2000)))
    acf_sa      <- theta_t * sa
    acf_ca      <- theta_t * ca
    acf_total   <- acf_sa + acf_ca
    acf_ta1     <- acf_total * (1 - gh_share)
    acf_ta2     <- acf_total * gh_share
    notif_total <- notif_passive + acf_total
    # State equations, including the three recovery compartments.
    d_s  <- b_t * n_pop - (lambda + mu_t) * s
    d_e  <- c_f * fl + c_s * sl - ((1 - e1 * e2) * lambda + mu_t) * e
    d_fl <- lambda * s + (1 - e1 * e2) * lambda * e + (1 - e1) *
      lambda * sl - (delta1 + delta2 + c_f + mu_t) * fl
    d_sl <- delta2 * fl - ((1 - e1) * lambda + delta3 + c_s + mu_t) * sl
    d_sa <- delta1 * fl + delta3 * sl + delta5 * ca + rel_flow +
      reinf_flow - (tao_sa + delta4 + mu_sa + mu_t) * sa - acf_sa
    d_ca <- delta4 * sa + omiga1 * p_lf1 * pa + omiga2 * p_lf2 * ga +
      (1 - alpha_t) * gamma * (ta1 + ta2) -
      (eta + tao_ca + delta5 + mu_ca + mu_t) * ca - acf_ca
    d_pa <- eta * p_pc * ca - (omiga1 + mu_ca + mu_t) * pa
    d_ga <- eta * p_gh * (1 - gh_cover) * ca + omiga1 * (1 - p_lf1) *
      p_ref * (1 - gh_cover) * pa - (omiga2 + mu_ca + mu_t) * ga
    d_ta1 <- eta * (1 - p_gh - p_pc) * ca +
      omiga1 * (1 - p_lf1) * (1 - p_ref) * pa +
      omiga2 * (1 - p_lf2) * (1 - gh_share) * ga +
      acf_ta1 - (gamma + mu_ta + mu_t) * ta1
    d_ta2 <- eta * p_gh * gh_cover * ca + omiga1 * (1 - p_lf1) * p_ref *
      gh_cover * pa + omiga2 * (1 - p_lf2) * gh_share *
      ga + acf_ta2 - (gamma + mu_ta + mu_t) * ta2
    d_r_low <- alpha_t * gamma * (ta1 + ta2) -
      (rho_low + kappa + reinf_rate + mu_t) * r_low
    d_r_high <- tao_sa * sa + tao_ca * ca -
      (rho_high + kappa + reinf_rate + mu_t) * r_high
    d_r_long <- kappa * (r_low + r_high) - (rho_long + reinf_rate + mu_t) *
      r_long
    tb_deaths   <- mu_sa * sa + mu_ca * (ca + pa + ga) + mu_ta * (ta1 + ta2)
    d_d         <- tb_deaths + mu_t * n_pop
    derivatives <- c(
      d_s, d_e, d_fl, d_sl, d_sa, d_ca, d_pa, d_ga, d_ta1, d_ta2,
      d_r_low, d_r_high, d_r_long, d_d
    )

    # Observation rates, incidence components and population counts.
    pop_scale    <- 12658.3
    population   <- n_pop * pop_scale
    notif_rate   <- notif_total / n_pop * 1e+05
    mort_rate    <- tb_deaths / n_pop * 1e+05
    mort_report  <- mort_rate * (1 - ur)
    igra_prev    <- (fl + sl + p_igra_e * e) / n_pop * 100
    tb_prev      <- atb / n_pop * 1e+05
    smear_prop   <- smear_share
    sa_prop      <- sa_share * 100
    acf_prop     <- acf_total / max(notif_total, .Machine$double.eps)
    inci_fast    <- delta1 * fl / n_pop * 1e+05
    inci_slow    <- delta3 * sl / n_pop * 1e+05
    rel_low  <- rho_low * r_low / n_pop * 1e+05
    rel_high <- rho_high * r_high / n_pop * 1e+05
    rel_long <- rho_long * r_long / n_pop * 1e+05
    reinfect     <- reinf_flow / n_pop * 1e+05
    inci_rate    <- inci_fast + inci_slow + rel_low +
      rel_high + rel_long + reinfect
    inci_num  <- inci_rate / 1e+05 * population
    death_num <- tb_deaths * pop_scale
    notif_num <- notif_total * pop_scale
    # Expected system delay among successfully treated patients.
    w0        <- (1 - p_gh - p_pc) + p_gh * gh_cover
    w1        <- p_pc * (1 - p_lf1) * (1 - p_ref + p_ref * gh_cover)
    w12       <- p_pc * (1 - p_lf1) * p_ref * (1 - gh_cover) * (1 - p_lf2)
    w2        <- p_gh * (1 - gh_cover) * (1 - p_lf2)
    p_success <- w0 + w1 + w12 + w2
    delay_days <- if (p_success > 0) {
      365 * (w1 / omiga1 + w12 * (1 / omiga1 + 1 / omiga2) + w2 / omiga2) /
        p_success
    } else {
      NA_real_
    }
    # Calibration needs only the target indicators and incidence components.
    if (calibration_only) {
      outcomes <- c(
        # Disease rates, counts and prevalence.
        inci_rate   = inci_rate,    # TB disease incidence /100,000
        notif_rate  = notif_rate,   # TB notifications /100,000
        mort_report = mort_report,  # Reported TB deaths /100,000
        tb_prev     = tb_prev,      # Active TB prevalence /100,000
        igra_prev   = igra_prev,    # IGRA-positive prevalence, %

        # Sources of incident TB disease, per 100,000.
        inci_fast   = inci_fast,    # Fast-latent progression
        inci_slow   = inci_slow,    # Slow-latent reactivation
        rel_low     = rel_low,      # Relapse from R_low
        rel_high    = rel_high,     # Relapse from R_high
        rel_long    = rel_long,     # Relapse from R_long
        reinfect    = reinfect,     # Disease after reinfection in R

        # Population and disease characteristics.
        r_total     = r,            # Recovered population, model scale
        sa_prop     = sa_prop,      # Subclinical share of active TB, %
        smear_prop  = smear_prop,   # Smear-positive share of active TB

        # Hospital coverage and care pathways.
        acf_prop    = acf_prop,     # ACF share of notifications
        delay_days  = delay_days    # Mean health system delay, days
      )
      return(list(derivatives, outcomes))
    }

    # Full analysis adds costs, health losses and care-pathway outcomes.
    cost_dx <- cost_dx_rs * (1 - resistance_share) +
      cost_dx_rr * resistance_share
    cost_tx <- cost_tx_rs * (1 - resistance_share) +
      cost_tx_rr * resistance_share
    cost_diagnosis <- notif_num * mult_dx * cost_dx
    recover_num    <- alpha_t * gamma * (ta1 + ta2) * pop_scale
    uncured_num    <- (1 - alpha_t) * gamma * (ta1 + ta2) * pop_scale
    cost_treatment <- (recover_num + 0.5 * uncured_num) * mult_tx * cost_tx
    cost_dt        <- cost_diagnosis + cost_treatment
    yld_num        <- ca_all * base_par$dw_tb * pop_scale
    yll_num        <- death_num * base_par$yll_tb
    care_num       <- eta * ca * pop_scale
    direct_tx      <- eta * (1 - p_gh - p_pc + p_gh * gh_cover) * ca * pop_scale
    ref_num      <- (omiga1 * pa + omiga2 * ga) * pop_scale
    ltfu_num       <- (omiga1 * pa * p_lf1 + omiga2 * ga * p_lf2) * pop_scale
    ref_need       <- eta * (p_gh + p_pc - p_gh * gh_cover) * ca * pop_scale
    outcomes <- c(
      # Disease rates, counts and prevalence.
      inci_rate   = inci_rate,    # TB disease incidence /100,000
      inci_num    = inci_num,     # Incident TB episodes /year
      notif_rate  = notif_rate,   # TB notifications /100,000
      notif_num   = notif_num,    # TB notifications /year
      mort_rate   = mort_rate,    # TB deaths /100,000
      death_num   = death_num,    # TB deaths /year
      mort_report = mort_report,  # Reported TB deaths /100,000
      tb_prev     = tb_prev,      # Active TB prevalence /100,000
      igra_prev   = igra_prev,    # IGRA-positive prevalence, %

      # Sources of incident TB disease, per 100,000.
      inci_fast   = inci_fast,    # Fast-latent progression
      inci_slow   = inci_slow,    # Slow-latent reactivation
      rel_low     = rel_low,      # Relapse from R_low
      rel_high    = rel_high,     # Relapse from R_high
      rel_long    = rel_long,     # Relapse from R_long
      reinfect    = reinfect,     # Disease after reinfection in R

      # Population and disease characteristics.
      population  = population,   # Living population, national scale
      r_total     = r,            # Recovered population, model scale
      sa_prop     = sa_prop,      # Subclinical share of active TB, %
      smear_prop  = smear_prop,   # Smear-positive share of active TB
      smear_sa    = smear_sa,     # Smear-positive share in SA
      smear_ca    = smear_ca,     # Smear-positive share in clinical TB
      mu_sa       = mu_sa,        # TB mortality rate in SA, /year
      mu_ca       = mu_ca,        # TB mortality rate in clinical TB, /year

      # Hospital coverage and care pathways.
      gh_cover    = gh_cover,     # Patient-weighted GH coverage
      hosp_cover  = hosp_cover,   # Actual GH fraction providing TB services
      gh_share    = gh_share,     # GH share of TB service facilities
      acf_prop    = acf_prop,     # ACF share of notifications
      delay_days  = delay_days,   # Mean health system delay, days
      care_num    = care_num,     # Care-seeking events /year
      direct_tx   = direct_tx,    # Direct treatment entries /year
      ref_num     = ref_num,      # Referral exits from PA/GA /year
      ltfu_num    = ltfu_num,     # Referral losses to follow-up /year
      ref_need    = ref_need,     # Patients needing referral /year

      # Costs and health losses.
      cost_dt     = cost_dt,      # Diagnosis and treatment cost, USD/year
      yld_num     = yld_num,      # YLDs accrued /year
      yll_num     = yll_num       # YLLs accrued /year
    )
    list(derivatives, outcomes)
  }
}

# ------------------------------------------------------------------------------
# 1.3 Cost-effectiveness analysis
# ------------------------------------------------------------------------------

run_cea <- function(calib_par, base_par, return_models = FALSE,
                    assumptions = "basecase") {
  # Annual evaluation years and continuous-time integration boundaries.
  eval_years <- settings$eval_years
  int_times  <- seq(min(eval_years), max(eval_years) + 1L, by = 1)

  # Include the end boundary of the final evaluation year.
  times <- sort(unique(c(settings$years, int_times)))

  y0 <- initialize_population(calib_par)

  run_model <- function(scenario_id) {
    as.data.frame(deSolve::ode(
      y     = y0,
      times = times,
      func  = tb_model(
        calib_par, base_par, scenarios = scenario_id,
        assumptions = assumptions
      )
    ))
  }

  models    <- lapply(0:3, run_model)
  eval_rows <- match(eval_years, times)
  int_rows  <- match(int_times, times)

  disc_years <- 1 / (1 + base_par$disc_rate)^(eval_years - 2024)
  disc_int   <- 1 / (1 + base_par$disc_rate)^(int_times - 2024)

  # Four absolute outcomes, followed by three differences from status quo.
  with_diffs <- function(values, change) {
    changes <- if (change == "increase") {
      values[-1] - values[1]
    } else {
      values[1] - values[-1]
    }
    c(values, changes)
  }

  # Integrate annual flows over the evaluation period.
  total <- function(metric, discount = FALSE,
                    change = c("increase", "reduction")) {
    change <- match.arg(change)

    weights <- if (discount) disc_int else rep(1, length(int_times))

    values <- sapply(models, function(model) {
      pracma::trapz(int_times, model[int_rows, metric] * weights)
    })

    with_diffs(values, change)
  }

  # Annual rates and other point-in-time outcomes.
  at_year <- function(year, metric, scale = 1,
                      change = c("increase", "reduction")) {
    change <- match.arg(change)
    row    <- match(year, times)

    if (is.na(row)) {
      stop("Requested year is not available in the model output.")
    }

    values <- sapply(models, function(model) {
      model[row, metric] * scale
    })

    with_diffs(values, change)
  }

  # Change over time, followed by differences from status quo.
  abs_change <- function(y2015, y2035) {
    changes <- y2035[1:4] - y2015[1:4]

    c(changes, changes[-1] - changes[1])
  }

  # Positive percentages indicate reductions over time.
  pct_red <- function(y2015, y2035) {
    reductions <- (y2015[1:4] - y2035[1:4]) / y2015[1:4] * 100

    c(reductions, reductions[-1] - reductions[1])
  }

  # General-hospital expansion costs.
  gh_cost <- function(scenario, baseline) {
    # China had 20,133 general hospitals in 2020.
    cover_diff <- scenario[eval_rows, "hosp_cover"] -
      baseline[eval_rows, "hosp_cover"]
    extra_gh   <- pmax(0, cover_diff * 20133)

    new_gh     <- c(extra_gh[1], pmax(0, diff(extra_gh)))
    equip_cost <- base_par$cost_equip / base_par$equip_life
    # Annualised equipment and programme costs continue while services operate.
    # Training is charged once, when each additional hospital starts services.
    annual_cost <- (
      cumsum(new_gh) * equip_cost + new_gh * base_par$cost_train
    ) / (1 - base_par$prog_share)

    sum(annual_cost * disc_years)
  }

  gh_totals <- c(0, sapply(models[-1], gh_cost, baseline = models[[1]]))

  # Epidemiological and care-engagement outcomes.
  inci_2015    <- at_year(2015, "inci_rate", change = "reduction")
  inci_2035    <- at_year(2035, "inci_rate", change = "reduction")
  deaths_2015  <- at_year(2015, "death_num", change = "reduction")
  deaths_2035  <- at_year(2035, "death_num", change = "reduction")
  inci_red     <- pct_red(inci_2015, inci_2035)
  death_red    <- pct_red(deaths_2015, deaths_2035)
  inci_total   <- total("inci_num", change = "reduction")
  deaths_total <- total("death_num", change = "reduction")
  notif_total  <- total("notif_num", change = "increase")
  direct_tx    <- total("direct_tx", change = "increase")
  ref_num      <- total("ref_num", change = "reduction")
  ltfu_total   <- total("ltfu_num", change = "reduction")

  # Active case-finding proportions are expressed as percentages.
  acf_2015   <- at_year(2015, "acf_prop", scale = 100, change = "reduction")
  acf_2035   <- at_year(2035, "acf_prop", scale = 100, change = "reduction")
  acf_change <- abs_change(acf_2015, acf_2035)
  delay_2015 <- at_year(2015, "delay_days", change = "reduction")
  delay_2035 <- at_year(2035, "delay_days", change = "reduction")
  delay_red  <- pct_red(delay_2015, delay_2035)

  # Costs, DALYs, and cost-effectiveness outcomes.
  cost_dt <- total("cost_dt", discount = TRUE, change = "increase")
  ylds    <- total("yld_num", discount = TRUE, change = "reduction")
  ylls    <- total("yll_num", discount = TRUE, change = "reduction")
  dalys   <- ylds + ylls
  cost_gh <- c(gh_totals, gh_totals[-1])
  costs   <- cost_dt + cost_gh

  # Three times China's GDP per capita in 2024 (US$13,445).
  wtp <- 13445 * 3

  cost_per_case  <- costs / inci_total
  cost_per_notif <- costs / notif_total
  cost_per_death <- costs / deaths_total
  # Absolute columns: DALY burden; incremental columns: DALYs averted.
  nmb  <- c(-wtp * dalys[1:4] - costs[1:4], wtp * dalys[5:7] - costs[5:7])
  icer <- costs / dalys

  seq_icer <- c(
    NA_real_, NA_real_, NA_real_, NA_real_,
    icer[5],
    (costs[3] - costs[2]) / (dalys[2] - dalys[3]),
    (costs[4] - costs[3]) / (dalys[3] - dalys[4])
  )

  # Legacy "infection" row keys refer to incident TB disease, not infection.
  # Keep these keys compatible with saved RDS results; exports use "case".
  cea <- rbind(
    Incidence_2015                       = inci_2015,
    Deaths_2015                          = deaths_2015,
    Incidence_2035                       = inci_2035,
    Incidence_reduction_2015_2035        = inci_red,
    Deaths_2035                          = deaths_2035,
    Death_reduction_2015_2035            = death_red,
    Cumulative_infections_averted        = inci_total,
    Cumulative_deaths_averted            = deaths_total,
    ACF_proportion_2015                  = acf_2015,
    Mean_system_delay_2015               = delay_2015,
    ACF_proportion_2035                  = acf_2035,
    ACF_change_2015_2035                 = acf_change,
    Mean_system_delay_2035               = delay_2035,
    D_sys_reduction_2015_2035            = delay_red,
    Reduced_referrals                    = ref_num,
    Reduced_loss_to_follow_up            = ltfu_total,
    Direct_treatment_difference          = direct_tx,
    Total_discounted_cost                = costs,
    Discounted_cost_DT                   = cost_dt,
    Discounted_expansion_cost            = cost_gh,
    Total_DALYs                          = dalys,
    DALYs_from_YLDs                      = ylds,
    DALYs_from_YLLs                      = ylls,
    Cost_per_infection_averted           = cost_per_case,
    Cost_per_additional_case             = cost_per_notif,
    Cost_per_death_averted               = cost_per_death,
    Net_monetary_benefit                 = nmb,
    ICER_USD_per_DALY_averted            = icer,
    Sequential_ICER_USD_per_DALY_averted = seq_icer
  )

  cea <- cea[, c(1, 2, 5, 3, 6, 4, 7)]

  colnames(cea) <- c(
    "status_quo",
    "scenario_50",
    "incremental_50",
    "scenario_80",
    "incremental_80",
    "scenario_100",
    "incremental_100"
  )

  if (return_models) {
    list(CEA_baseline = cea, models = models)
  } else {
    cea
  }
}

# ------------------------------------------------------------------------------
# 1.4 Model calibration
# ------------------------------------------------------------------------------

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)

# Worker budget shared by calibration and subsequent analyses.
detected_cores <- parallel::detectCores(logical = TRUE)
parallel_workers <- if (is.na(detected_cores)) {
  1L
} else {
  max(1L, detected_cores - 2L)
}

run_calibration <- function(p_igra_e, base_par, calib_settings,
                            assumptions = "basecase") {
  # Fixed IGRA assumption and MCMC controls for this calibration.
  calibration_start    <- Sys.time()
  base_par_one         <- base_par
  base_par_one$p_igra_E <- p_igra_e
  nchains <- calib_settings$nchains
  niter   <- calib_settings$niter
  burnin  <- niter / 2
  thin    <- 50L
  covprop <- 0.05
  seed    <- settings$seed

  # Uniform parameter bounds and initial proposal covariance.
  bounds <- rbind(
    Beta1      = c(0.01,   1.2),    # Transmission coefficient, beta
    delta2     = c(0.2,    0.9),    # Transition rate from FL to SL
    theta_2023 = c(0.0001, 0.2),    # Case-finding rate in SA/CA, 2023
    eta        = c(1,      4),      # Total care-seeking rate
    omiga1     = c(1,      3),      # Total referral rate: PA to GA/TA1/TA2
    omiga2     = c(1,      3),      # Total referral rate: GA to TA1/TA2
    Ur         = c(0.1,    1),      # Proportion of TB deaths not reported
    a          = c(0.01,   0.3),    # Initial excess over long-term asymptote
    b          = c(0.01,   1),      # Smear-positive exponential decline rate
    c          = c(0.01,   0.2),    # Long-term smear-positive asymptote
    Pre_E      = c(0.01,   0.3),    # Eliminated infection prevalence, 2000
    Pre_LTBI   = c(0.1,    0.4),    # Latent infection prevalence, 2000
    Prop_FL    = c(0.01,   0.3),    # FL proportion among LTBI, 2000
    Pre_ATB    = c(0.007,  0.015),  # ATB prevalence including SA and CA, 2000
    Prop_SA    = c(0.2,    0.5),    # SA proportion among ATB, 2000
    Prop_PA    = c(0.01,   0.15),   # PA proportion among ATB, 2000
    Prop_GA    = c(0.01,   0.15),   # GA proportion among ATB, 2000
    Prop_TA    = c(0.01,   0.15)    # TA1 proportion among ATB, 2000
  )
  colnames(bounds) <- c("lower", "upper")
  bounds   <- as.data.frame(bounds)
  parnames <- rownames(bounds)
  lower    <- bounds$lower
  upper    <- bounds$upper
  jump_cov <- diag(((upper - lower) * covprop)^2) * (2.4^2 / length(lower))

  # Notifications per 100,000, 2006-2020; assumed SD is 5% of each value.
  notif_obs <- c(
    106.64, 105.12, 99.69,  89.14,  81.94,
    75.97,  73.76,  69.43,  68.08,  65.10,
    61.67,  60.74,  59.38,  55.57,  47.58
  )
  # Reported TB mortality per 100,000, 2004-2010; lower/upper 95% limits.
  death_obs <- c(7.03, 7.13, 5.43, 5.04, 4.04, 4.33, 4.38)
  death_lo  <- c(6.84, 6.93, 5.26, 4.88, 3.90, 4.19, 4.23)
  death_hi  <- c(7.23, 7.32, 5.60, 5.20, 4.18, 4.48, 4.52)
  death_sd  <- (death_hi - death_lo) / 3.92
  acf_obs   <- 0.03    # ACF share of notifications in 2023
  acf_sd    <- 0.0015  # Assumed SD for the ACF target

  likelihood_env <- list2env(list(
    base_par_one          = base_par_one,
    assumptions           = assumptions,
    initialize_population = initialize_population,
    tb_model              = tb_model,
    resistance_prop       = resistance_prop,
    notif_obs             = notif_obs,
    death_obs             = death_obs,
    death_sd              = death_sd,
    acf_obs               = acf_obs,
    acf_sd                = acf_sd
  ), parent = .GlobalEnv)
  log_lik <- function(params) {
    out <- deSolve::ode(
      y     = initialize_population(params),
      times = 2000:2023,
      func  = tb_model(params, base_par_one, assumptions = assumptions,
                       calibration_only = TRUE)
    )
    required_columns <- c(
      "time", "notif_rate", "mort_report", "igra_prev", "tb_prev",
      "sa_prop", "delay_days", "smear_prop", "acf_prop"
    )
    if (!all(required_columns %in% colnames(out))) {
      stop("tb_model is missing required calibration outputs.")
    }
    years <- out[, "time"]
    if (nrow(out) != 24L || anyNA(years) || !all(years == 2000:2023)) {
      return(-Inf)
    }

    # Map model outputs to observation years and survey populations.
    notif_pred <- out[match(2006:2020, years), "notif_rate"]
    death_pred <- out[match(2004:2010, years), "mort_report"]
    delay_pred <- out[match(2019, years), "delay_days"]
    igra_p     <- 1.0357 * out[match(2013, years), "igra_prev"] / 100
    tb_p       <- 1.1246 * out[match(2010, years), "tb_prev"] / 1e5
    sa_p       <- out[match(2010, years), "sa_prop"] / 100
    smear_p    <- out[match(c(2000, 2010), years), "smear_prop"]
    acf_p      <- out[match(2023, years), "acf_prop"]

    # Reject invalid ODE predictions before evaluating their likelihood.
    preds <- c(notif_pred, death_pred, delay_pred)
    probs <- c(igra_p, tb_p, sa_p, smear_p, acf_p)
    if (any(!is.finite(c(preds, probs))) || any(preds < 0) ||
          any(probs < 0 | probs > 1)) {
      return(-Inf)
    }

    # Normal likelihoods for rates/delay; binomial for survey counts.
    components <- c(
      # L1: notifications, not disease incidence; 2006-2020.
      L1 = sum(dnorm(notif_obs, notif_pred, notif_obs * 0.05, log = TRUE)),
      # L2: reported TB mortality, 2004-2010.
      L2 = sum(dnorm(death_obs, death_pred, death_sd, log = TRUE)),
      # L3: IGRA positivity at age 5+, 2013.
      L3 = dbinom(3955, size = 21022, prob = igra_p, log = TRUE),
      # L4: active TB prevalence at age 15+, 2010.
      L4 = dbinom(1310, size = 252940, prob = tb_p, log = TRUE),
      # L5: subclinical share of active TB, 2010.
      L5 = dbinom(561, size = 1301, prob = sa_p, log = TRUE),
      # L6: mean health system delay in days, 2019.
      L6 = dnorm(99.88, delay_pred, sd = 21.66, log = TRUE),
      # L7: smear-positive share of active TB, 2000 and 2010.
      L7 = sum(dbinom(c(447, 181), c(1340, 1310), smear_p, log = TRUE)),
      # L8: ACF share of notifications, 2023.
      L8 = dnorm(acf_obs, acf_p, sd = acf_sd, log = TRUE)
    )
    sum(components)
  }
  environment(log_lik) <- likelihood_env
  cl <- parallel::makeCluster(min(nchains * 5L, parallel_workers), outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  invisible(parallel::clusterEvalQ(cl, library(FME)))
  parallel::clusterSetRNGStream(cl, iseed = seed)
  parallel::clusterExport(cl, c(
    "log_lik", "parnames", "lower", "upper", "niter", "burnin",
    "jump_cov"
  ), envir = environment())
  parallel::clusterExport(cl, c(
    "base_par_one", "initialize_population", "tb_model", "resistance_prop"
  ), envir = likelihood_env)
  # Twenty pre-fit jobs share the worker budget; MCMC uses four chains.
  run_prefit <- function(job) {
    assign(
      ".Random.seed", job$stream, # nolint: object_name_linter.
      envir = .GlobalEnv
    )
    started <- Sys.time()
    span    <- upper - lower
    decode  <- function(z) setNames(lower + z * span, parnames)
    valid   <- function(p) {
      p["Pre_ATB"] + p["Pre_LTBI"] + p["Pre_E"] <= 1 &&
        sum(p[c("Prop_SA", "Prop_PA", "Prop_GA", "Prop_TA")]) <= 1
    }
    objective <- function(z) {
      p <- decode(z)
      if (!valid(p)) {
        return(1e100)
      }
      v <- -2 * log_lik(p)
      if (is.finite(v)) v else 1e100
    }
    for (attempt in seq_len(1000L)) {
      z0 <- runif(length(lower))
      v0 <- objective(z0)
      if (v0 < 1e100) break
    }
    if (v0 >= 1e100) stop("No valid prefit start")
    fit <- stats::optim(z0, objective,
      method = "L-BFGS-B", lower = rep(0, length(lower)),
      upper  = rep(1, length(lower)), control = list(maxit = 500)
    )
    z       <- if (is.finite(fit$value) && fit$value < v0) fit$par else z0
    value   <- if (is.finite(fit$value)) min(v0, fit$value) else v0
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    cat(sprintf(
      "[Chain %d] Pre-fit %d/5 | Objective %.3f | Code %d | %.1f min\n",
      job$chain, job$start, value, fit$convergence, seconds / 60
    ))
    flush.console()
    list(
      chain   = job$chain,
      start   = job$start,
      p       = decode(z),
      value   = value,
      seconds = seconds,
      code    = fit$convergence
    )
  }
  environment(run_prefit) <- .GlobalEnv
  # Independent streams keep pre-fits and chains reproducible in parallel.
  original_kind <- RNGkind()
  original_seed <- get(".Random.seed", envir = .GlobalEnv)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  stream      <- get(".Random.seed", envir = .GlobalEnv)
  prefit_jobs <- list()
  for (chain in seq_len(nchains)) {
    for (start in 1:5) {
      prefit_jobs[[length(prefit_jobs) + 1L]] <- list(
        chain = chain, start = start, stream = stream
      )
      stream <- parallel::nextRNGStream(stream)
    }
  }
  mcmc_streams <- list()
  for (chain in seq_len(nchains)) {
    mcmc_streams[[chain]] <- stream
    stream <- parallel::nextRNGStream(stream)
  }
  do.call(RNGkind, as.list(original_kind))
  assign(
    ".Random.seed", original_seed, # nolint: object_name_linter.
    envir = .GlobalEnv
  )
  prefit_wall_start <- Sys.time()
  message(
    "Starting ", length(prefit_jobs), " pre-fit optimisations with ",
    length(cl), " workers."
  )
  fits                <- parallel::parLapplyLB(cl, prefit_jobs, run_prefit)
  prefit_wall_seconds <- as.numeric(difftime(
    Sys.time(), prefit_wall_start, units = "secs"
  ))
  chain_starts        <- lapply(seq_len(nchains), function(id) {
    candidates <- fits[vapply(fits, function(f) f$chain == id, logical(1))]
    candidates[[which.min(vapply(candidates, "[[", numeric(1), "value"))]]
  })
  # Switch from parallel pre-fits to one worker per MCMC chain.
  parallel::stopCluster(cl)
  cl <- parallel::makeCluster(min(nchains, parallel_workers), outfile = "")
  invisible(parallel::clusterEvalQ(cl, library(FME)))
  parallel::clusterExport(cl, c(
    "log_lik", "parnames", "lower", "upper", "niter", "burnin",
    "jump_cov", "chain_starts", "mcmc_streams"
  ), envir = environment())
  parallel::clusterExport(cl, c(
    "base_par_one", "initialize_population", "tb_model", "resistance_prop"
  ), envir = likelihood_env)
  run_chain <- function(chain_id) {
    assign(
      ".Random.seed", mcmc_streams[[chain_id]], # nolint: object_name_linter.
      envir = .GlobalEnv
    )
    valid_allocation <- function(p) {
      p["Pre_ATB"] + p["Pre_LTBI"] + p["Pre_E"] <= 1 &&
        sum(p[c("Prop_SA", "Prop_PA", "Prop_GA", "Prop_TA")]) <= 1
    }
    p <- chain_starts[[chain_id]]$p
    cat(sprintf(
      "[Chain %d] Starting MCMC from pre-fit objective %.3f\n",
      chain_id, chain_starts[[chain_id]]$value
    ))
    evaluations <- 0L
    chain_start <- Sys.time()
    objective   <- function(p) {
      value <- if (valid_allocation(p)) {
        -2 * log_lik(p)
      } else {
        Inf
      }
      evaluations <<- evaluations + 1L # nolint: assignment_linter.
      if (evaluations %% 10000L == 0L) {
        elapsed <- as.numeric(difftime(Sys.time(), chain_start, units = "mins"))
        cat(sprintf(
          "[Chain %d] %s objective evaluations | Elapsed: %.1f min\n",
          chain_id, format(evaluations, big.mark = ","), elapsed
        ))
        flush.console()
      }
      value
    }
    result <- FME::modMCMC(
      f            = objective,
      p            = p,
      lower        = lower,
      upper        = upper,
      niter        = niter,
      burninlength = burnin,
      updatecov    = 100,
      ntrydr       = 2,
      jump         = jump_cov,
      verbose      = 10000
    )
    cat(sprintf(
      "[Chain %d] MCMC finished in %.1f minutes.\n", chain_id,
      as.numeric(difftime(Sys.time(), chain_start, units = "mins"))
    ))
    result$performance <- list(
      prefit_seconds = chain_starts[[chain_id]]$seconds,
      mcmc_seconds   = as.numeric(difftime(
        Sys.time(), chain_start, units = "secs"
      ))
    )
    result
  }
  environment(run_chain) <- .GlobalEnv
  message("Running ", nchains, " MCMC chains. This may take a long time.")
  # Summarize retained draws; thin only the posterior prediction samples.
  mcmc_chains  <- parallel::parLapply(cl, seq_len(nchains), run_chain)
  mcmc_list <- coda::mcmc.list(lapply(mcmc_chains, function(ch) {
    coda::mcmc(ch$pars)
  }))
  chains_pars <- do.call(rbind, lapply(mcmc_chains, function(ch) ch$pars))
  a1_param_sum <- data.frame(
    Mean     = colMeans(chains_pars),
    par_q025 = apply(chains_pars, 2, quantile, probs = 0.025),
    par_q975 = apply(chains_pars, 2, quantile, probs = 0.975)
  )
  diagnostic_error <- function(e) list(error = conditionMessage(e))
  a2_gelman_diag <- tryCatch(
    coda::gelman.diag(
      mcmc_list, autoburnin = FALSE, multivariate = FALSE
    )$psrf,
    error = diagnostic_error
  )
  a3_effective_size <- tryCatch(
    coda::effectiveSize(mcmc_list), error = diagnostic_error
  )
  chain_posts       <- do.call(rbind, lapply(mcmc_chains, function(ch) {
    ch$pars[seq.int(1L, nrow(ch$pars), by = thin), , drop = FALSE]
  }))
  message("Calculating autocorrelation data...")
  autocorrelation_data <- do.call(rbind, lapply(
    seq_along(mcmc_chains), function(j) {
      draws <- as.matrix(mcmc_chains[[j]]$pars)
      do.call(rbind, lapply(seq_len(ncol(draws)), function(i) {
        x <- draws[, i]
        if (length(x) < 2L || any(!is.finite(x)) || stats::sd(x) == 0) {
          warning("ACF unavailable for ", colnames(draws)[i], " in Chain ", j)
          return(data.frame(
            Lag = 0, ACF = NA_real_, Parameter = colnames(draws)[i],
            Chain = paste("Chain", j)
          ))
        }
        ac <- stats::acf(x, lag.max = min(1000L, length(x) - 1L), plot = FALSE)
        data.frame(
          Lag       = as.numeric(ac$lag),
          ACF       = as.numeric(ac$acf),
          Parameter = colnames(draws)[i],
          Chain     = paste("Chain", j)
        )
      }))
    }
  ))
  # Posterior checks include all targets and the incidence breakdown.
  posterior_start <- Sys.time()
  plot_years      <- 2000:2025
  target_names <- c(
    "L1: TB notifications",
    "L2: Reported TB mortality",
    "L3: IGRA positivity, age 5+",
    "L4: Active TB prevalence, age 15+",
    "L5: Subclinical TB proportion",
    "L6: Mean health system delay",
    "L7: Smear-positive TB proportion",
    "L8: Active case-finding proportion"
  )
  message("Preparing calibration plot data through 2025...")
  predict_trajectory <- function(i) {
    p   <- chain_posts[i, ]
    out <- as.data.frame(deSolve::ode(
      y     = initialize_population(p),
      times = plot_years,
      func  = tb_model(p, base_par_one, assumptions = assumptions,
                       calibration_only = TRUE)
    ))
    required_columns <- c(
      "time", "notif_rate", "mort_report", "igra_prev", "tb_prev",
      "sa_prop", "delay_days", "smear_prop", "acf_prop"
    )
    if (!all(required_columns %in% names(out)) ||
          nrow(out) != length(plot_years) || anyNA(out$time) ||
          !all(out$time == plot_years)) {
      stop("Incomplete or incompatible ODE output for posterior sample ", i)
    }
    values <- cbind(
      out$notif_rate, out$mort_report, out$igra_prev * 1.0357,
      out$tb_prev * 1.1246, out$sa_prop, out$delay_days,
      out$smear_prop * 100, out$acf_prop * 100
    )
    if (any(!is.finite(values)) || any(values < 0)) {
      stop("Invalid predictions for posterior sample ", i)
    }
    attr(values, "incidence") <- out$inci_rate
    attr(values, "components") <- as.matrix(out[, c(
      "inci_fast", "inci_slow", "rel_low", "rel_high", "rel_long", "reinfect"
    )])
    values
  }
  environment(predict_trajectory) <- .GlobalEnv
  prediction_workers <- min(parallel_workers, nrow(chain_posts))
  if (length(cl) != prediction_workers) {
    parallel::stopCluster(cl)
    cl <- parallel::makeCluster(prediction_workers)
  }
  parallel::clusterExport(cl, c(
    "chain_posts", "plot_years", "base_par_one", "assumptions",
    "initialize_population",
    "tb_model", "resistance_prop"
  ), envir = environment())
  trajectory_list <- pbapply::pblapply(
    seq_len(nrow(chain_posts)), predict_trajectory, cl = cl
  )
  trajectory_array <- array(unlist(trajectory_list, use.names = FALSE), dim = c(
    length(plot_years),
    length(target_names), nrow(chain_posts)
  ))
  pred_mean       <- apply(trajectory_array, c(1, 2), mean)
  pred_low        <- apply(trajectory_array, c(1, 2), quantile, probs = 0.025)
  pred_high       <- apply(trajectory_array, c(1, 2), quantile, probs = 0.975)
  calibration_plot_data <- do.call(rbind, lapply(
    seq_along(target_names), function(i) {
      data.frame(
        Year      = plot_years,
        Mean      = pred_mean[, i],
        Low       = pred_low[, i],
        High      = pred_high[, i],
        Indicator = target_names[i]
      )
    }
  ))
  calib_par_mean <- setNames(a1_param_sum$Mean, rownames(a1_param_sum))
  incidence_draws <- vapply(
    trajectory_list, attr, numeric(length(plot_years)), which = "incidence"
  )
  incidence_prediction <- data.frame(
    Year      = plot_years,
    Mean      = rowMeans(incidence_draws),
    Low_2.5   = apply(incidence_draws, 1, quantile, probs = 0.025),
    High_97.5 = apply(incidence_draws, 1, quantile, probs = 0.975)
  )
  incidence_components <- data.frame(
    Year = plot_years,
    Reduce("+", lapply(trajectory_list, attr, which = "components")) /
      length(trajectory_list)
  )
  component_total <- rowSums(incidence_components[, -1, drop = FALSE])
  stopifnot(max(abs(component_total - incidence_prediction$Mean)) < 1e-8)
  performance <- list(
    prefit_wall_seconds = prefit_wall_seconds,
    posterior_seconds = as.numeric(difftime(
      Sys.time(), posterior_start, units = "secs"
    )),
    total_calibration_seconds = as.numeric(difftime(
      Sys.time(), calibration_start, units = "secs"
    )),
    chains = lapply(mcmc_chains, "[[", "performance")
  )
  cea_baseline <- run_cea(calib_par_mean, base_par_one,
                          assumptions = assumptions)
  icer         <- cea_baseline[
    "ICER_USD_per_DALY_averted",
    c("incremental_50", "incremental_80", "incremental_100")
  ]

  result <- list(
    p_igra_E = p_igra_e,
    assumptions = assumptions,
    Known_par = base_par_one,
    calibration_settings = utils::modifyList(
      calib_settings,
      list(seed = seed, thin = thin, covprop = covprop, burnin = burnin)
    ),
    parameter_bounds = data.frame(
      lower = lower, upper = upper, row.names = parnames
    ),
    A1_param_sum              = a1_param_sum,
    A2_gelman_diag            = a2_gelman_diag,
    A3_effective_size         = a3_effective_size,
    chain_posts               = chain_posts,
    mcmc_list                 = mcmc_list,
    mcmc_chains               = mcmc_chains,
    calibration_plot_data     = calibration_plot_data,
    incidence_prediction      = incidence_prediction,
    autocorrelation_plot_data = autocorrelation_data,
    CEA_baseline              = cea_baseline,
    ICER                      = icer,
    incidence_components      = incidence_components,
    performance               = performance
  )
  result
}
# Calibrate IGRA assumptions concurrently within one shared worker budget.
# Each MCMC job needs nchains workers plus one coordinating R process.
calibrate_scenarios <- function(scenarios, base_par, calib_settings,
                                worker_budget = parallel_workers) {
  if (!length(scenarios)) return(list())
  concurrent <- min(length(scenarios), max(1L, floor(
    worker_budget / (calib_settings$nchains + 1L)
  )))
  if (concurrent == 1L) {
    return(setNames(lapply(scenarios, function(value) {
      run_calibration(value, base_par, calib_settings)
    }), as.character(scenarios)))
  }
  workers_each <- max(calib_settings$nchains,
                      floor(worker_budget / concurrent) - 1L)
  message(
    "Calibrating ", length(scenarios), " IGRA assumptions with ", concurrent,
    " concurrent jobs and ", workers_each, " workers per job."
  )
  cl <- parallel::makeCluster(concurrent, outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c(
    "run_calibration",
    "initialize_population", "tb_model", "resistance_prop",
    "run_cea", "settings"
  ), envir = environment())
  parallel::clusterCall(cl, function(budget) {
    parallel_workers <<- budget
    set.seed(settings$seed)
    NULL
  }, workers_each)
  results <- parallel::parLapply(
    cl, scenarios, function(value, known, controls) {
      run_calibration(value, known, controls)
    }, known = base_par, controls = calib_settings
  )
  setNames(results, as.character(scenarios))
}

# ------------------------------------------------------------------------------
# 1.5 One-way sensitivity analysis
# ------------------------------------------------------------------------------
run_usa <- function(a1_param_sum, base_par, par_ranges, num = 7,
                    assumptions = "basecase") {
  if (num < 3L || num %% 2 == 0L) {
    stop("num must be an odd integer of at least 3.")
  }

  calib_par_mean <- setNames(
    a1_param_sum$Mean,
    rownames(a1_param_sum)
  )

  calib_par_ranges <- setNames(
    lapply(seq_len(nrow(a1_param_sum)), function(i) {
      unname(unlist(
        a1_param_sum[i, c("par_q025", "par_q975")],
        use.names = FALSE
      ))
    }),
    rownames(a1_param_sum)
  )

  scen_cols <- c(
    "status_quo", "scenario_50",
    "scenario_80", "scenario_100"
  )

  incr_cols <- c(
    "incremental_50",
    "incremental_80",
    "incremental_100"
  )

  # Generate equally spaced low-to-base and base-to-high parameter values.

  make_values <- function(lower, base, upper) {
    n_left  <- (num + 1L) / 2L
    n_right <- num - n_left + 1L

    c(
      seq(lower, base, length.out = n_left),
      seq(base, upper, length.out = n_right)[-1]
    )
  }

  # Extract the outcomes required for the tornado plots.

  extract_results <- function(cea_result) {
    c(
      setNames(
        cea_result["Total_DALYs", scen_cols],
        paste0("DALY", 0:3)
      ),
      setNames(
        cea_result["Total_discounted_cost", scen_cols],
        paste0("COST", 0:3)
      ),
      setNames(
        cea_result["ICER_USD_per_DALY_averted", incr_cols],
        paste0("ICER", 1:3)
      )
    )
  }

  run_parameter <- function(parameter, values, parameter_type) {
    results <- lapply(values, function(value) {
      if (parameter_type == "Known") {
        base_par_one <- base_par
        base_par_one[[parameter]] <- value

        run_cea(calib_par = calib_par_mean, base_par = base_par_one,
                assumptions = assumptions)
      } else {
        calib_par_one <- calib_par_mean
        calib_par_one[[parameter]] <- value

        run_cea(calib_par = calib_par_one, base_par = base_par,
                assumptions = assumptions)
      }
    })

    outcome_matrix <- do.call(
      rbind,
      lapply(results, extract_results)
    )

    data.frame(
      parameter      = parameter,
      parameter_type = parameter_type,
      par_value      = values,
      value_n        = seq_along(values),
      outcome_matrix,
      row.names   = NULL,
      check.names = FALSE
    )
  }

  usa_results <- list()

  # Vary known parameters over their specified sensitivity ranges.

  known_parameters <- intersect(
    names(base_par),
    names(par_ranges)
  )

  for (parameter in known_parameters) {
    values <- make_values(
      par_ranges[[parameter]][1],
      base_par[[parameter]],
      par_ranges[[parameter]][2]
    )

    usa_results[[parameter]] <- run_parameter(
      parameter, values,
      parameter_type = "Known"
    )
  }

  # Vary calibrated parameters over their posterior 95% credible intervals.

  for (parameter in names(calib_par_mean)) {
    values <- make_values(
      calib_par_ranges[[parameter]][1],
      calib_par_mean[[parameter]],
      calib_par_ranges[[parameter]][2]
    )

    usa_results[[parameter]] <- run_parameter(
      parameter, values,
      parameter_type = "Calibrated"
    )
  }

  usa_results
}

# ------------------------------------------------------------------------------
# 1.6 Probabilistic sensitivity analysis
# ------------------------------------------------------------------------------

run_psa <- function(chain_posts, known_par_dist, plf2 = base_par$p_lf2,
                    assumptions = "basecase") {
  num_cores <- min(parallel_workers, nrow(chain_posts))

  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  pred_years <- settings$pred_years

  pred_vars <- c(
    "inci_rate", "inci_num",
    "mort_rate", "death_num",
    "acf_prop", "delay_days",
    "ref_num", "ltfu_num", "direct_tx"
  )

  scen_cols <- c(
    "status_quo", "scenario_50",
    "scenario_80", "scenario_100"
  )

  incr_cols <- c(
    "incremental_50",
    "incremental_80",
    "incremental_100"
  )

  simulate_one <- function(i) {
    calib_par_one <- chain_posts[i, ]
    base_par_one <- known_par_dist[i, ]
    base_par_one$p_lf2 <- plf2

    cea_output <- run_cea(
      calib_par = calib_par_one,
      base_par = base_par_one,
      return_models = TRUE,
      assumptions = assumptions
    )

    cea_result <- cea_output$CEA_baseline

    psa_result <- c(
      setNames(
        cea_result["ICER_USD_per_DALY_averted", incr_cols],
        paste0("ICER", 1:3)
      ),
      setNames(
        cea_result["Net_monetary_benefit", incr_cols],
        paste0("INMB", 1:3)
      ),
      setNames(
        cea_result["Total_discounted_cost", incr_cols],
        paste0("cost_incr", 1:3)
      ),
      setNames(
        cea_result["Total_DALYs", incr_cols],
        paste0("effect_incr", 1:3)
      ),
      setNames(
        cea_result["Total_DALYs", scen_cols],
        paste0("DALY", 0:3)
      ),
      setNames(
        cea_result["Total_discounted_cost", scen_cols],
        paste0("COST", 0:3)
      )
    )

    pred_result <- array(
      unlist(lapply(cea_output$models, function(model) {
        as.matrix(
          model[
            match(pred_years, model$time),
            pred_vars,
            drop = FALSE
          ]
        )
      }), use.names = FALSE),
      dim = c(
        length(pred_years),
        length(pred_vars),
        length(scen_cols)
      ),
      dimnames = list(
        Year     = pred_years,
        Outcome  = pred_vars,
        Scenario = scen_cols
      )
    )

    list(
      cea_result  = cea_result,
      psa_result  = psa_result,
      pred_result = pred_result
    )
  }

  environment(simulate_one) <- .GlobalEnv

  parallel::clusterExport(
    cl,
    c(
      "chain_posts", "known_par_dist", "plf2", "assumptions",
      "pred_years", "pred_vars", "scen_cols", "incr_cols",
      "initialize_population", "tb_model", "resistance_prop", "run_cea",
      "settings"
    ),
    envir = environment()
  )

  simulation_results <- parallel::parLapply(
    cl, seq_len(nrow(chain_posts)), simulate_one
  )

  # ----------------------------------------------------------------------------
  # Summarize model-based cost-effectiveness outcomes
  # ----------------------------------------------------------------------------

  cea_results  <- lapply(simulation_results, `[[`, "cea_result")
  cea_template <- cea_results[[1]]

  cea_array <- array(
    unlist(cea_results, use.names = FALSE),
    dim = c(
      nrow(cea_template),
      ncol(cea_template),
      length(cea_results)
    )
  )

  cea_mean <- apply(cea_array, c(1, 2), mean, na.rm = TRUE)
  cea_low  <- apply(
    cea_array, c(1, 2),
    stats::quantile,
    probs = 0.025, na.rm = TRUE
  )
  cea_high <- apply(
    cea_array, c(1, 2),
    stats::quantile,
    probs = 0.975, na.rm = TRUE
  )

  cea_final <- data.frame(
    Outcome = rownames(cea_template),
    setNames(
      as.data.frame(cea_mean),
      paste0("mean_", colnames(cea_template))
    ),
    setNames(
      as.data.frame(cea_low),
      paste0("low_2.5_", colnames(cea_template))
    ),
    setNames(
      as.data.frame(cea_high),
      paste0("high_97.5_", colnames(cea_template))
    ),
    check.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # Calculate cost-effectiveness acceptability curve data
  # ----------------------------------------------------------------------------

  psa_results <- data.frame(
    simulation = seq_len(nrow(chain_posts)),
    do.call(rbind, lapply(simulation_results, `[[`, "psa_result")),
    check.names = FALSE
  )

  wtp_thresholds <- seq(0, 50000, by = 100)

  dalys <- as.matrix(psa_results[paste0("DALY", 0:3)])
  costs <- as.matrix(psa_results[paste0("COST", 0:3)])

  ceac_probabilities <- t(sapply(wtp_thresholds, function(wtp) {
    best_strategy <- max.col(
      -dalys * wtp - costs,
      ties.method = "first"
    )

    tabulate(best_strategy, nbins = 4) / length(best_strategy)
  }))

  ceac_data <- data.frame(
    wtp = wtp_thresholds,
    setNames(
      as.data.frame(ceac_probabilities),
      c(
        "prob_status_quo",
        "prob_intervention1",
        "prob_intervention2",
        "prob_intervention3"
      )
    ),
    check.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # Summarize annual model predictions from 2015 to 2035
  # ----------------------------------------------------------------------------

  pred_results <- lapply(simulation_results, `[[`, "pred_result")

  pred_array <- array(
    unlist(pred_results, use.names = FALSE),
    dim = c(
      length(pred_years),
      length(pred_vars),
      length(scen_cols),
      length(pred_results)
    )
  )

  pred_mean   <- apply(pred_array, c(1, 2, 3), mean, na.rm = TRUE)
  pred_median <- apply(pred_array, c(1, 2, 3), median, na.rm = TRUE)

  pred_low <- apply(
    pred_array, c(1, 2, 3),
    stats::quantile,
    probs = 0.025, na.rm = TRUE
  )

  pred_high <- apply(
    pred_array, c(1, 2, 3),
    stats::quantile,
    probs = 0.975, na.rm = TRUE
  )

  pred_data <- do.call(rbind, lapply(seq_along(scen_cols), function(s) {
    do.call(rbind, lapply(seq_along(pred_vars), function(v) {
      data.frame(
        Year     = pred_years,
        Scenario = scen_cols[s],
        Outcome  = pred_vars[v],
        Mean     = pred_mean[, v, s],
        Median   = pred_median[, v, s],
        Low_2.5   = pred_low[, v, s],
        High_97.5 = pred_high[, v, s]
      )
    }))
  }))

  list(
    cea_final = cea_final,
    ceac_data = ceac_data,
    pred_data = pred_data
  )
}

# ------------------------------------------------------------------------------
# 1.7 Complete analysis under one care-engagement assumption
# ------------------------------------------------------------------------------
run_analysis <- function(calibration, base_par, par_ranges,
                         assumptions = "basecase") {
  base_par$p_igra_E <- calibration$p_igra_E
  calib_mean <- setNames(
    calibration$A1_param_sum$Mean, rownames(calibration$A1_param_sum)
  )
  cea_baseline <- run_cea(calib_mean, base_par, assumptions = assumptions)
  chain_posts  <- calibration$chain_posts

  # Include IGRA retention uncertainty in sensitivity analyses.
  usa_ranges <- par_ranges

  usa_data <- run_usa(
    a1_param_sum     = calibration$A1_param_sum,
    base_par        = base_par,
    par_ranges = usa_ranges,
    assumptions = assumptions
  )

  known_par_dist_df <- sample_known_par(
    base_par        = base_par,
    par_ranges = par_ranges,
    n    = nrow(chain_posts),
    seed = settings$seed
  )

  psa_base <- run_psa(
    chain_posts, known_par_dist_df,
    plf2 = base_par$p_lf2, assumptions = assumptions
  )
  psa_sen60 <- run_psa(chain_posts, known_par_dist_df, plf2 = 0.60,
                       assumptions = assumptions)
  psa_sen45 <- run_psa(chain_posts, known_par_dist_df, plf2 = 0.45,
                       assumptions = assumptions)
  psa_sen30 <- run_psa(chain_posts, known_par_dist_df, plf2 = 0.30,
                       assumptions = assumptions)
  psa_sen15 <- run_psa(chain_posts, known_par_dist_df, plf2 = 0.15,
                       assumptions = assumptions)

  list(
    assumptions  = assumptions,
    seed         = settings$seed,
    p_igra_E     = calibration$p_igra_E,
    CEA_baseline = cea_baseline,
    usa_data     = usa_data,
    psa_base     = psa_base,
    psa_sen60    = psa_sen60,
    psa_sen45    = psa_sen45,
    psa_sen30    = psa_sen30,
    psa_sen15    = psa_sen15
  )

}

# Each care assumption has its own posterior and complete analysis.
run_care_scenarios <- function(cases, base_par, par_ranges, calib_settings,
                               worker_budget = parallel_workers) {
  checkpoint_dir <- file.path(save_dir, "process", "care_checkpoints")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  concurrent <- min(length(cases), max(1L, worker_budget %/% 2L))
  workers_each <- max(1L, worker_budget %/% concurrent - 1L)
  cl <- parallel::makeCluster(concurrent, outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c(
    "run_calibration", "initialize_population",
    "tb_model", "resistance_prop", "run_cea", "run_analysis", "run_usa",
    "run_psa", "sample_known_par", "settings"
  ), envir = environment())
  parallel::clusterCall(cl, function(budget) {
    parallel_workers <<- budget
    set.seed(settings$seed)
    NULL
  }, workers_each)
  results <- parallel::parLapply(cl, cases, function(assumption, base, ranges,
                                                   controls, checkpoint) {
    message("Starting care assumption: ", assumption)
    calibration <- run_calibration(0, base, controls, assumption)
    path <- file.path(checkpoint, paste0(assumption, ".rds"))
    saveRDS(list(calibration = calibration, analysis = NULL), path)
    analysis <- run_analysis(calibration, base, ranges, assumption)
    result <- list(calibration = calibration, analysis = analysis)
    saveRDS(result, path)
    result
  }, base = base_par, ranges = par_ranges, controls = calib_settings,
     checkpoint = checkpoint_dir)
  setNames(results, cases)
}

# ==============================================================================
# S2. Model calibration and main analysis
# ==============================================================================
# ------------------------------------------------------------------------------
# 2.1 Model calibration and comparison of IGRA assumptions
# ------------------------------------------------------------------------------
if (run_all) {
  if (!is.null(model_results$calibration) || !is.null(model_results$analysis)) {
    entry <- paste0("previous_run_", length(model_results$history) + 1L)
    model_results$history[[entry]] <- list(
      seed        = model_results$seed,
      calibration = model_results$calibration,
      analysis    = model_results$analysis,
      care        = model_results$care,
      source      = model_results$source,
      saved_at    = model_results$saved_at
    )
  }
  # Record the environment used for computation, rather than later plotting.
  model_results$source <- list(
    engine = "R", R_version = R.version.string, started_at = start_time,
    package_versions = vapply(
      required_packages, function(pkg) as.character(utils::packageVersion(pkg)),
      character(1)
    )
  )
  calibration_results <- calibrate_scenarios(
    igra_scenarios, base_par, settings$calib
  )
}

for (p_igra_e in igra_scenarios) {
  key    <- as.character(p_igra_e)
  result <- calibration_results[[key]]

  if (is.null(result)) {
    stop("Missing calibration result for p_igra_E = ", p_igra_e)
  }
  message("\nCalibration available for p_igra_E = ", p_igra_e)

  if (is.null(result$CEA_baseline)) {
    scenario_par <- base_par
    scenario_par$p_igra_E <- p_igra_e
    scenario_calib <- setNames(
      result$A1_param_sum$Mean, rownames(result$A1_param_sum)
    )
    result$CEA_baseline <- run_cea(scenario_calib, scenario_par)
    calibration_results[[key]] <- result
  }
  cat("\nCalibration diagnostics: p_igra_E =", p_igra_e, "\n")

  par_sum   <- result$A1_param_sum
  rhat_diag <- result$A2_gelman_diag
  ess_diag  <- result$A3_effective_size

  # Diagnostic failures are stored as error lists.
  if (!is.matrix(rhat_diag) || !is.numeric(ess_diag) ||
    !all(c("Point est.", "Upper C.I.") %in% colnames(rhat_diag))) {
    cat("\nConvergence assessment unavailable. Diagnostic outputs:\n")
    print(rhat_diag)
    print(ess_diag)
  } else {
    # Align diagnostics with the parameter-summary order.
    parameters <- rownames(par_sum)
    par_index <- match(parameters, rownames(rhat_diag))
    rhat      <- rhat_diag[par_index, "Point est."]
    rhat_ucl  <- rhat_diag[par_index, "Upper C.I."]
    ess        <- unname(ess_diag[parameters])

    summary_table <- data.frame(
      Parameter = parameters,
      Mean      = round(par_sum$Mean, 4),
      `95% CrI` = sprintf(
        "%.4f\u2013%.4f", par_sum$par_q025, par_sum$par_q975
      ),
      Rhat            = round(rhat, 3),
      Rhat_upper_95CI = round(rhat_ucl, 3),
      ESS             = round(ess),
      row.names       = NULL,
      check.names     = FALSE
    )

    print(summary_table, row.names = FALSE)

    # Assess unrounded diagnostics.
    if (any(!is.finite(rhat)) || any(!is.finite(ess))) {
      cat("\nAssessment: Inconclusive; missing or non-finite diagnostics.\n")
    } else {
      if (all(rhat < 1.01)) {
        cat("\nR-hat: All parameters meet the stricter < 1.01 criterion.\n")
      } else if (all(rhat < 1.10)) {
        cat(
          "\nR-hat: All parameters meet the selected < 1.10 criterion.",
          "\nParameters with R-hat >= 1.01: ",
          paste(parameters[rhat >= 1.01], collapse = ", "),
          "\n",
          sep = ""
        )
      } else {
        cat(
          "\nR-hat: Further convergence review is needed.",
          "\nParameters with R-hat >= 1.10: ",
          paste(parameters[rhat >= 1.10], collapse = ", "),
          "\n",
          sep = ""
        )
      }

      if (any(ess < 400)) {
        cat(
          "\nESS: Parameters below 400: ",
          paste(parameters[ess < 400], collapse = ", "),
          "\n",
          sep = ""
        )
      } else {
        cat("\nESS: All parameters have ESS >= 400.\n")
      }
    }

    cat(
      "\nReview trace and autocorrelation plots alongside these diagnostics.",
      "\nAssess model fit separately using the calibration plot.\n"
    )
  }
}

# Save newly computed calibration before analysis; keep earlier result files.
if (run_all) {
  model_results$seed <- settings$seed
  model_results$calibration <- list(
    igra_scenarios = igra_scenarios,
    seed           = settings$seed,
    saved_at       = Sys.time(),
    results        = calibration_results
  )
  model_results$analysis <- tb_analysis_results
  model_results$saved_at <- Sys.time()
  results_file <- file.path(save_dir, paste0(
    "TB_results_seed", settings$seed, "_",
    format(Sys.Date(), "%Y%m%d"), ".rds"
  ))
  saveRDS(model_results, results_file)
  message("Saved: ", results_file)
}
# Use the primary IGRA assumption for the full analysis.
selected_p_igra_e <- base_par$p_igra_E

cat(
  "\nPrimary analysis: p_igra_E =",
  selected_p_igra_e, "\n"
)

# Select the primary calibration for diagnostics and subsequent analyses.
tb_calibration_results <- calibration_results[[
  as.character(selected_p_igra_e)
]]

a1_param_sum <- tb_calibration_results$A1_param_sum
mcmc_chains  <- tb_calibration_results$mcmc_chains
chain_posts  <- tb_calibration_results$chain_posts
cea_baseline <- tb_calibration_results$CEA_baseline

base_par$p_igra_E <- selected_p_igra_e

# ------------------------------------------------------------------------------
# 2.2 Run model analysis
# ------------------------------------------------------------------------------

if (run_all) {
  message("Starting model analysis. This may take a long time.")
  tb_analysis_results <- run_analysis(
    tb_calibration_results, base_par, par_ranges
  )
  model_results$analysis <- tb_analysis_results
  model_results$saved_at <- Sys.time()
  saveRDS(model_results, results_file)
} else {
  message("Using saved model-analysis results.")
}

# ------------------------------------------------------------------------------
# 2.3 Calibrate and analyse best and worst care-engagement assumptions
# ------------------------------------------------------------------------------
if (run_all) model_results$care <- NULL
care_cases <- c("best", "worst")
missing_care <- care_cases[vapply(care_cases, function(case) {
  is.null(model_results$care[[case]]$calibration) ||
    is.null(model_results$care[[case]]$analysis)
}, logical(1))]
if (length(missing_care)) {
  care_results <- run_care_scenarios(
    missing_care, base_par, par_ranges, settings$calib
  )
  model_results$care[missing_care] <- care_results
  model_results$saved_at <- Sys.time()
  saveRDS(model_results, results_file)
}

cat(
  "\nElapsed time through S2:",
  round(difftime(Sys.time(), start_time, units = "mins"), 2),
  "minutes\n"
)

# ==============================================================================
# S3. Create tables and figures
# ==============================================================================
# ------------------------------------------------------------------------------
# 3.1 Create and export result tables
# ------------------------------------------------------------------------------

# Compare cost-effectiveness across IGRA assumptions.
calibration_comparison <- do.call(
  rbind,
  lapply(igra_scenarios, function(p_igra_e) {
    result <- calibration_results[[as.character(p_igra_e)]]

    data.frame(
      p_igra_E = p_igra_e,
      ICER_50  = result$CEA_baseline[
        "ICER_USD_per_DALY_averted", "incremental_50"
      ],
      ICER_80 = result$CEA_baseline[
        "ICER_USD_per_DALY_averted", "incremental_80"
      ],
      ICER_100 = result$CEA_baseline[
        "ICER_USD_per_DALY_averted", "incremental_100"
      ],
      row.names = NULL
    )
  })
)

cat("\nICER comparison across IGRA assumptions:\n")
print(calibration_comparison, row.names = FALSE)

# Prepare numeric comparison data for the shared Excel workbook.
scenario_cols <- c(
  "status_quo", "scenario_50", "scenario_80", "scenario_100"
)
incremental_cols <- c(
  "incremental_50", "incremental_80", "incremental_100"
)
comparison_table <- do.call(
  rbind,
  lapply(igra_scenarios, function(p_igra_e) {
    cea <- calibration_results[[as.character(p_igra_e)]]$CEA_baseline
    data.frame(
      p_igra_E                  = p_igra_e,
      Strategy                  = scenario_cols,
      Total_DALYs_million = unname(cea["Total_DALYs", scenario_cols]) / 1e6,
      Total_cost_USD_billion =
        unname(cea["Total_discounted_cost", scenario_cols]) / 1e9,
      ICER_USD_per_DALY_averted = c(
        NA_real_, unname(cea["ICER_USD_per_DALY_averted", incremental_cols])
      ),
      row.names = NULL
    )
  })
)

# Summarize IGRA assumptions for Supplementary Table 11.
table11 <- do.call(rbind, lapply(igra_scenarios, function(p) {
  values <- comparison_table[comparison_table$p_igra_E == p, ]
  section <- data.frame(
    Strategy = paste0("Assumed p_igra_E = ", p),
    DALYs = "", Cost = "", ICER = ""
  )
  rows <- data.frame(
    Strategy = c("Status quo", "50% coverage", "80% coverage", "100% coverage"),
    DALYs = sprintf("%.2f", values$Total_DALYs_million),
    Cost = sprintf("%.2f", values$Total_cost_USD_billion),
    ICER = ifelse(is.na(values$ICER_USD_per_DALY_averted), "NA",
                  formatC(values$ICER_USD_per_DALY_averted,
                          format = "f", digits = 2, big.mark = ","))
  )
  rbind(section, rows)
}))
names(table11) <- c(
  "Strategy", "Total DALYs (million)", "Total cost (US$ billion)",
  "ICER (US$ per DALY averted)"
)

# Combine parameter diagnostics for all IGRA assumptions.
calibration_diagnostic_tables <- lapply(
  calibration_results,
  function(diagnostic_result) {
    par_sum   <- diagnostic_result$A1_param_sum
    rhat_diag <- diagnostic_result$A2_gelman_diag
    ess_diag  <- diagnostic_result$A3_effective_size

    parameters <- rownames(par_sum)
    # Preserve parameter summaries even when convergence diagnostics fail.
    rhat              <- rhat_ucl <- ess <- rep(NA_real_, length(parameters))
    diagnostic_errors <- character()
    if (is.matrix(rhat_diag) &&
      all(c("Point est.", "Upper C.I.") %in% colnames(rhat_diag))) {
      diagnostic_index <- match(parameters, rownames(rhat_diag))
      rhat             <- rhat_diag[diagnostic_index, "Point est."]
      rhat_ucl         <- rhat_diag[diagnostic_index, "Upper C.I."]
    } else {
      diagnostic_errors <- c(
        diagnostic_errors,
        paste("R-hat diagnostic output:",
              paste(capture.output(str(rhat_diag)), collapse = " "))
      )
    }
    if (is.numeric(ess_diag)) {
      ess <- unname(ess_diag[parameters])
    } else {
      diagnostic_errors <- c(
        diagnostic_errors,
        paste("ESS diagnostic output:",
              paste(capture.output(str(ess_diag)), collapse = " "))
      )
    }

    # Format the exported parameter table.
    selected_diagnostics <- data.frame(
      Parameter = parameters,
      Mean      = sprintf("%.4f", par_sum$Mean),
      `95%CrI` = sprintf(
        "%.4f\u2013%.4f",
        par_sum$par_q025,
        par_sum$par_q975
      ),
      Rhat            = sprintf("%.3f", rhat),
      Rhat_upper_95CI = sprintf("%.3f", rhat_ucl),
      ESS             = round(ess),
      row.names       = NULL,
      check.names     = FALSE
    )

    diagnostic_notes <- c(
      paste(
        "Criterion: R-hat < 1.1 suggests adequate convergence;",
        "< 1.01 indicates stronger evidence of convergence."
      ),
      "Criterion: ESS >= 400 indicates adequate sampling precision.",
      if (any(!is.finite(rhat))) {
        paste(
          "R-hat: Assessment unavailable;",
          "some values are missing or non-finite."
        )
      } else if (all(rhat < 1.01)) {
        "R-hat: All parameters meet the stricter < 1.01 criterion."
      } else if (all(rhat < 1.10)) {
        "R-hat: All parameters meet the selected < 1.10 criterion."
      } else {
        "R-hat: Further convergence review is needed."
      },
      if (any(!is.finite(ess))) {
        paste(
          "ESS: Assessment unavailable;",
          "some values are missing or non-finite."
        )
      } else if (all(ess >= 400)) {
        "ESS: All parameters have ESS >= 400."
      } else {
        "ESS: Some parameters have ESS < 400."
      }
    )

    bounds <- diagnostic_result$parameter_bounds
    if (!is.null(bounds)) {
      bounds <- bounds[match(parameters, rownames(bounds)), , drop = FALSE]
      relative_mean <- (par_sum$Mean - bounds$lower) /
        (bounds$upper - bounds$lower)
      near_boundary <- parameters[
        which(relative_mean <= 0.05 | relative_mean >= 0.95)
      ]
      if (length(near_boundary)) {
        diagnostic_notes <- c(diagnostic_notes, paste(
          "Review posterior means within 5% of prior bounds:",
          paste(near_boundary, collapse = ", ")
        ))
      }
    } else {
      diagnostic_notes <- c(
        diagnostic_notes,
        "Boundary assessment unavailable: no saved parameter bounds."
      )
    }
    diagnostic_notes <- c(diagnostic_notes, diagnostic_errors)
    # Append a blank row followed by diagnostic notes.
    note_rows <- selected_diagnostics[
      rep(NA_integer_, length(diagnostic_notes) + 1L), ,
      drop = FALSE
    ]
    note_rows$Parameter <- c("", diagnostic_notes)

    rbind(selected_diagnostics, note_rows)
  }
)
# Combine the three calibration scenarios in one diagnostic table.
calibration_diagnostics <- do.call(rbind, lapply(
  seq_along(calibration_diagnostic_tables), function(i) {
    data.frame(
      p_igra_E = as.numeric(names(calibration_results)[i]),
      calibration_diagnostic_tables[[i]],
      row.names = NULL,
      check.names = FALSE
    )
  }
))

# Use one workbook for calibration and analysis tables.
results_excel <- file.path(save_dir, paste0(
  format(Sys.Date(), "%Y%m%d"), "_TB_results_seed", settings$seed, ".xlsx"
))
calibration_sheets <- list(Calibration_diagnostics = calibration_diagnostics)
# Extract results for tables and figures.
cea_baseline <- tb_analysis_results$CEA_baseline
usa_data     <- tb_analysis_results$usa_data
psa_base     <- tb_analysis_results$psa_base
psa_sen60    <- tb_analysis_results$psa_sen60
psa_sen45    <- tb_analysis_results$psa_sen45
psa_sen30    <- tb_analysis_results$psa_sen30
psa_sen15    <- tb_analysis_results$psa_sen15

# Export baseline CEA and PSA uncertainty results to one workbook.
cea_baseline_table <- data.frame(
  Outcome = sub("infection", "case", rownames(cea_baseline), fixed = TRUE),
  cea_baseline,
  row.names   = NULL,
  check.names = FALSE
)

# Build Table 1 and the corresponding sensitivity-analysis tables.
# Costs, DALYs and CEA point estimates use the deterministic baseline.
# Other point estimates and all displayed 95% UIs use the PSA summary.
table1_rows <- c(
  "Epidemiological outcomes" = "",
  "TB incidence in 2035 (per 100,000)" = "Incidence_2035",
  "Reduction in TB incidence compared with 2015 (%)" =
    "Incidence_reduction_2015_2035",
  "TB deaths in 2035 (millions)" = "Deaths_2035",
  "Reduction in TB deaths compared with 2015 (%)" =
    "Death_reduction_2015_2035",
  "Cumulative incident TB cases, 2026–2035 (millions)" =
    "Cumulative_infections_averted",
  "Cumulative TB deaths, 2026–2035 (millions)" =
    "Cumulative_deaths_averted",
  "Care engagement outcomes" = "",
  "Cumulative first-visit diagnosis and treatment initiation" =
    "Direct_treatment_difference",
  "Cumulative referrals for TB cases, 2026–2035 (millions)" =
    "Reduced_referrals",
  "Cumulative losses to follow-up for TB cases, 2026–2035 (millions)" =
    "Reduced_loss_to_follow_up",
  "Mean health system delay in 2035 (days)" = "Mean_system_delay_2035",
  "Reduction in mean health system delay compared with 2015 (%)" =
    "D_sys_reduction_2015_2035",
  "Cost, 2026–2035 (US$, billions)" = "",
  "Total cost" = "Total_discounted_cost",
  "Cost for TB diagnosis and treatment" = "Discounted_cost_DT",
  "Cost for implementation" = "Discounted_expansion_cost",
  "DALYs, 2026–2035 (millions)" = "",
  "Total" = "Total_DALYs",
  "YLDs" = "DALYs_from_YLDs",
  "YLLs" = "DALYs_from_YLLs",
  "Cost-effectiveness outcomes" = "",
  "Incremental cost per TB case averted (US$)" =
    "Cost_per_infection_averted",
  "Incremental cost per TB death averted (US$)" = "Cost_per_death_averted",
  "Incremental net monetary benefit (US$, billions)" =
    "Net_monetary_benefit",
  "ICER (US$ per DALY averted)" = "ICER_USD_per_DALY_averted",
  "Frontier ICER (US$ per DALY averted)" =
    "Sequential_ICER_USD_per_DALY_averted"
)

names(table1_rows)[table1_rows == "Direct_treatment_difference"] <- paste0(
  "Cumulative first-visit diagnosis and treatment initiation for TB cases, ",
  "2026–2035 (millions)"
)

# Conversion factors and display precision follow the manuscript table.
cea_keys <- c(
  "Cost_per_infection_averted", "Cost_per_death_averted",
  "Net_monetary_benefit", "ICER_USD_per_DALY_averted",
  "Sequential_ICER_USD_per_DALY_averted"
)
cost_keys <- c(
  "Total_discounted_cost", "Discounted_cost_DT", "Discounted_expansion_cost"
)
daly_keys <- c("Total_DALYs", "DALYs_from_YLDs", "DALYs_from_YLLs")
cum_keys <- c("Cumulative_infections_averted", "Cumulative_deaths_averted")
care_keys <- c(
  "Reduced_referrals", "Reduced_loss_to_follow_up",
  "Direct_treatment_difference"
)
pct_keys <- c(
  "Incidence_reduction_2015_2035", "Death_reduction_2015_2035",
  "D_sys_reduction_2015_2035"
)
# Apply the same source and display rules to all three analysis tables.
make_results_table <- function(analysis, title) {
  cea_baseline <- analysis$CEA_baseline
  psa_base     <- analysis$psa_base
  table1_cols <- colnames(cea_baseline)

  # Keep changes paired with outcomes from the same PSA simulations.
  table1_values <- lapply(unname(table1_rows), function(key) {
    if (key == "") return(rep("", length(table1_cols)))
    scale <- if (key %in% c(cost_keys, "Net_monetary_benefit")) {
      1e9
    } else if (key %in% c(daly_keys, cum_keys, care_keys)) {
      1e6
    } else if (key == "Deaths_2035") {
      1e6
    } else {
      1
    }
    digits <- if (key %in% c(cum_keys, "Deaths_2035")) 2L else 1L
    psa_row <- psa_base$cea_final[
      match(key, psa_base$cea_final$Outcome), , drop = FALSE
    ]
    fmt <- function(x) {
      formatC(x / scale, format = "f", digits = digits, big.mark = ",")
    }
    vapply(table1_cols, function(col) {
      is_change <- startsWith(col, "incremental_")
      if (key == "Discounted_expansion_cost" && col == "status_quo") {
        return("0")
      }
      if (key %in% cea_keys && !is_change) return("—")
      if (key == "Sequential_ICER_USD_per_DALY_averted" &&
          col == "incremental_50") return("—")
      point <- if (key %in% c(cost_keys, daly_keys, cea_keys)) {
        cea_baseline[key, col]
      } else {
        psa_row[[paste0("mean_", col)]]
      }
      value <- fmt(point)
      if (key %in% pct_keys && !is_change) value <- paste0(value, "%")
      if (key %in% cea_keys) return(value)
      paste0(
        value, " (", fmt(psa_row[[paste0("low_2.5_", col)]]), "–",
        fmt(psa_row[[paste0("high_97.5_", col)]]), ")"
      )
    }, character(1))
  })

  # Preserve the manuscript's outcome order and seven result columns.
  table1 <- data.frame(
    Outcome = names(table1_rows),
    do.call(rbind, table1_values),
    row.names = NULL, check.names = FALSE
  )
  names(table1)[-1] <- c(
    "Status quo (13%): Outcome",
    "50%: Outcome", "50%: Change vs status quo",
    "80%: Outcome", "80%: Change vs status quo",
    "100%: Outcome", "100%: Change vs status quo"
  )
  table1_title <- title

  table1_notes <- c(
    table1_title,
    paste0("Notes: Parentheses indicate 95% uncertainty intervals. Coverage ",
           "continues along its existing trend under the status quo; ",
           "expansion reaches 50%, 80% or 100% by 2030, then stays constant."),
    paste0("Changes are relative to the status quo; differences in percentage ",
           "reductions are expressed in percentage points."),
    paste0("First-visit diagnosis and treatment initiation means TB is ",
           "diagnosed and treatment started at the first care-seeking visit."),
    paste0("Net monetary benefit assumes US$40,335 per DALY averted. Frontier ",
           "ICERs compare 80% versus 50% and 100% versus 80% coverage."),
    paste0("Point estimates for costs, DALYs and CEA use CEA_baseline; other ",
           "point estimates and all UIs use PSA_cea_final."),
    paste0("Abbreviations: TB, tuberculosis; ",
           "DALY, disability-adjusted life-year; ",
           "YLD, year lived with disability; YLL, year of life lost; ",
           "ICER, incremental cost-effectiveness ratio.")
  )
  if (identical(analysis$assumptions, "best")) {
    table1$Outcome[table1$Outcome == "Frontier ICER (US$ per DALY averted)"] <-
      "Sequential ICER (US$ per DALY averted)"
    table1_notes <- sub("Frontier ICERs", "Sequential ICERs", table1_notes,
                        fixed = TRUE)
    table1_notes <- c(table1_notes, paste0(
      "The 50% and 80% coverage scenarios are subject to extended dominance; ",
      "the cost-effectiveness frontier connects status quo and 100% coverage."
    ))
  }
  note_rows <- as.data.frame(matrix(
    "", nrow = length(table1_notes) + 1L, ncol = ncol(table1)
  ))
  names(note_rows) <- names(table1)
  note_rows$Outcome <- c("", table1_notes)
  rbind(table1, note_rows)
}

table1 <- make_results_table(tb_analysis_results, paste0(
  "Table 1. Projected health outcomes and cost-effectiveness of expanding ",
  "TB service coverage of general hospitals in China, 2026–2035."
))
care_tables <- list(
  Table12_Best = make_results_table(model_results$care$best$analysis, paste0(
    "Supplementary Table 12. Projected health outcomes and cost-effectiveness ",
    "under the best care-engagement assumption."
  )),
  Table13_Worst = make_results_table(model_results$care$worst$analysis, paste0(
    "Supplementary Table 13. Projected health outcomes and cost-effectiveness ",
    "under the worst care-engagement assumption."
  ))
)

psa_export <- psa_base$cea_final
psa_export$Outcome <- sub("infection", "case", psa_export$Outcome, fixed = TRUE)
writexl::write_xlsx(
  c(
    list(Table1 = table1, Table11 = table11),
    care_tables,
    calibration_sheets,
    list(
      CEA_baseline           = cea_baseline_table,
      PSA_cea_final          = psa_export
    )
  ),
  path = results_excel
)

message("Saved: ", results_excel)

# ------------------------------------------------------------------------------
# 3.2 Generate calibration figures: Supplementary Figures 3-6
# ------------------------------------------------------------------------------

autocorrelation_data <- tb_calibration_results$autocorrelation_plot_data
if (is.null(autocorrelation_data)) {
  stop("The calibration results are missing cached autocorrelation data.")
}

parnames     <- colnames(mcmc_chains[[1]]$pars)
nchains      <- length(mcmc_chains)
chain_labels <- paste("Chain", seq_len(nchains))
chain_colors <- setNames(
  grDevices::colorRampPalette(
    c("#2878B5", "#D9534F", "#39A77D", "#9467BD")
  )(nchains),
  chain_labels
)

# Shared plotting styles.
diagnostic_theme <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    legend.position  = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

figure_text_theme <- ggplot2::theme(
  text         = ggplot2::element_text(size = 16),
  axis.title   = ggplot2::element_text(size = 14),
  axis.text    = ggplot2::element_text(size = 11),
  strip.text   = ggplot2::element_text(size = 14),
  plot.title   = ggplot2::element_text(size = 16, hjust = 0),
  plot.caption = ggplot2::element_text(size = 12),
  legend.title = ggplot2::element_text(size = 14),
  legend.text  = ggplot2::element_text(size = 14),
  plot.margin  = ggplot2::margin(10, 18, 14, 12, unit = "pt")
)

# ESS combines information from all four chains.
ess_labels <- data.frame(
  Parameter = parnames,
  Label     = paste0(
    "ESS = ",
    formatC(
      unname(tb_calibration_results$A3_effective_size[parnames]),
      format = "f", digits = 0, big.mark = ","
    )
  )
)

# Autocorrelation plots.
message("Preparing autocorrelation plots from saved data...")

autocorrelation_data$Parameter <- factor(
  autocorrelation_data$Parameter,
  levels = parnames
)
autocorrelation_data$Chain <- factor(
  autocorrelation_data$Chain,
  levels = chain_labels
)

max_lag <- max(autocorrelation_data$Lag, na.rm = TRUE)

p_autocorrelation <- ggplot2::ggplot(
  autocorrelation_data,
  ggplot2::aes(Lag, ACF, colour = Chain)
) +
  ggplot2::geom_hline(
    yintercept = 0, colour = "grey60", linewidth = 0.3
  ) +
  ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
  ggplot2::scale_colour_manual(
    values = chain_colors, breaks = chain_labels
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(0, max_lag / 2, max_lag),
    expand = ggplot2::expansion(mult = c(0.02, 0.08))
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(0, 0.5, 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    x = "Lag", y = "Autocorrelation", colour = NULL
  ) +
  diagnostic_theme

# Trace plots and posterior distributions.
message("Preparing parameter traces and distributions...")

trace_data <- do.call(rbind, lapply(seq_len(nchains), function(j) {
  draws <- as.matrix(mcmc_chains[[j]]$pars)

  data.frame(
    Iteration = rep(seq_len(nrow(draws)), times = ncol(draws)),
    Value     = as.vector(draws),
    Parameter = rep(colnames(draws), each = nrow(draws)),
    Chain     = chain_labels[j]
  )
}))

trace_data$Parameter <- factor(
  trace_data$Parameter,
  levels = parnames
)
trace_data$Chain <- factor(
  trace_data$Chain,
  levels = chain_labels
)

p_trace <- ggplot2::ggplot(
  trace_data, ggplot2::aes(Iteration, Value, colour = Chain)
) +
  ggplot2::geom_line(linewidth = 0.25, alpha = 0.65) +
  ggplot2::scale_colour_manual(
    values = chain_colors, breaks = chain_labels
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_pretty(n = 3),
    labels = function(x) {
      format(
        x,
        scientific = FALSE, trim = TRUE, big.mark = ","
      )
    },
    expand = ggplot2::expansion(mult = c(0.02, 0.08))
  ) +
  ggplot2::labs(
    x = "Iteration", y = "Parameter value", colour = NULL
  ) +
  diagnostic_theme

# Use the same 10,000 posterior draws as subsequent simulations.
posterior_data <- data.frame(
  Value     = as.vector(as.matrix(chain_posts)),
  Parameter = factor(
    rep(colnames(chain_posts), each = nrow(chain_posts)),
    levels = parnames
  )
)
p_distribution <- ggplot2::ggplot(
  posterior_data, ggplot2::aes(Value)
) +
  ggplot2::geom_histogram(
    bins = 40, fill = "#A9CCE3", colour = "white", linewidth = 0.2
  ) +
  ggplot2::facet_wrap(~Parameter, scales = "free", ncol = 4) +
  ggplot2::labs(
    x     = "Parameter value", y = "Frequency",
    title = "Marginal posterior parameter distributions"
  ) +
  diagnostic_theme

# Arrange parameter panels and place the legend in the spare cells.
parameter_legend_layout <- function(
  plot, plot_data, plot_title, show_ess = FALSE
) {
  n_parameters <- length(parnames)
  n_columns    <- 4L

  panels <- lapply(seq_along(parnames), function(i) {
    parameter <- parnames[i]

    panel_data <- plot_data[
      as.character(plot_data$Parameter) == parameter, ,
      drop = FALSE
    ]

    # Show axis titles only on the left and bottom exposed panels.
    left_column  <- (i - 1L) %% n_columns == 0L
    bottom_panel <- i + n_columns > n_parameters

    panel <- (plot + panel_data) +
      ggplot2::labs(
        title = parameter,
        x     = if (bottom_panel) {
          if (show_ess) "Lag" else "Iteration"
        } else {
          NULL
        },
        y = if (left_column) {
          if (show_ess) "Autocorrelation" else "Parameter value"
        } else {
          NULL
        }
      )

    if (show_ess) {
      panel <- panel +
        ggplot2::geom_text(
          data = ess_labels[
            ess_labels$Parameter == parameter, ,
            drop = FALSE
          ],
          ggplot2::aes(label = Label),
          x     = Inf, y = Inf,
          hjust = 1.08, vjust = 1.3,
          size  = 6,
          inherit.aes = FALSE
        )
    }

    panel
  })

  n_cells     <- ceiling((n_parameters + 1L) / n_columns) * n_columns
  panel_codes <- c(
    LETTERS[seq_len(n_parameters)],
    rep(LETTERS[n_parameters + 1L], n_cells - n_parameters)
  )

  design <- paste(
    vapply(seq.int(1L, n_cells, by = n_columns), function(i) {
      paste0(
        panel_codes[i:(i + n_columns - 1L)],
        collapse = ""
      )
    }, character(1)),
    collapse = "\n"
  )

  combined <- patchwork::wrap_plots(
    c(panels, list(patchwork::guide_area())),
    design = design,
    guides = "collect"
  ) +
    patchwork::plot_annotation(title = plot_title)

  combined &
    figure_text_theme &
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "horizontal"
    )
}

p_autocorrelation <- parameter_legend_layout(
  p_autocorrelation,
  autocorrelation_data,
  "Autocorrelation of posterior parameter samples",
  show_ess = TRUE
)

p_trace <- parameter_legend_layout(
  p_trace,
  trace_data,
  "Trace plots for the MCMC chains"
)

p_distribution <- p_distribution + figure_text_theme

# Save the three diagnostic figures.
calibration_figures <- list(
  Supplementary_Fig3_MCMC_trace                            = p_trace,
  Supplementary_Fig4_MCMC_autocorrelation = p_autocorrelation,
  Supplementary_Fig5_calibrated_parameter_posteriors       = p_distribution
)

for (figure_name in names(calibration_figures)) {
  figure <- calibration_figures[[figure_name]]
  message("Rendering and saving: ", figure_name)

  print(figure)

  filename <- paste0(figure_name, ".png")

  ggplot2::ggsave(
    filename = file.path(save_dir, filename),
    plot     = figure,
    width    = 14,
    height   = 16,
    units    = "in",
    dpi      = 300,
    bg       = "white"
  )

  message("Saved: ", filename)
}

model_data <- tb_calibration_results$calibration_plot_data
model_data$Indicator <- sub(
  "Mean system delay", "Mean health system delay", model_data$Indicator
)
if (is.null(model_data)) {
  stop("The calibration results are missing cached model predictions.")
}

target_names <- c(
  "L1: TB notifications",
  "L2: Reported TB mortality",
  "L3: IGRA positivity, age 5+",
  "L4: Active TB prevalence, age 15+",
  "L5: Subclinical TB proportion",
  "L6: Mean health system delay",
  "L7: Smear-positive TB proportion",
  "L8: Active case-finding proportion"
)
target_units <- c(
  "Rate per 100,000",
  "Rate per 100,000",
  "Proportion (%)",
  "Prevalence per 100,000",
  "Proportion (%)",
  "Days",
  "Proportion (%)",
  "Proportion (%)"
)

if (!all(target_names %in% model_data$Indicator)) {
  stop("Saved model indicator names do not match the plotting labels.")
}

# L1: assumed normal error intervals; L2: reported confidence intervals.
# L3/L4/L5/L7: exact binomial confidence intervals.
# L6: working normal interval for the inferred mean, not the reported IQR.

notice_obs <- c(
  106.64, 105.12, 99.69, 89.14, 81.94, 75.97, 73.76, 69.43,
  68.08, 65.10, 61.67, 60.74, 59.38, 55.57, 47.58
)
mortality_obs <- c(7.03, 7.13, 5.43, 5.04, 4.04, 4.33, 4.38)
mortality_lo  <- c(6.84, 6.93, 5.26, 4.88, 3.90, 4.19, 4.23)
mortality_hi  <- c(7.23, 7.32, 5.60, 5.20, 4.18, 4.48, 4.52)

# L8 is a working normal calibration target based on the national
# estimate that approximately 3% of notifications were actively detected.

acf_target <- 0.03
acf_sd     <- 0.0015

binomial_target <- function(indicator, years, counts, sizes, scale) {
  intervals <- t(vapply(seq_along(counts), function(i) {
    as.numeric(stats::binom.test(counts[i], sizes[i])$conf.int)
  }, numeric(2)))

  data.frame(
    Year      = years,
    Mean      = counts / sizes * scale,
    Low       = intervals[, 1] * scale,
    High      = intervals[, 2] * scale,
    Indicator = indicator
  )
}

target_data <- rbind(
  data.frame(
    Year      = 2006:2020, Mean = notice_obs,
    Low       = notice_obs * (1 - 1.96 * 0.05),
    High      = notice_obs * (1 + 1.96 * 0.05),
    Indicator = target_names[1]
  ),
  data.frame(
    Year      = 2004:2010, Mean = mortality_obs,
    Low       = mortality_lo, High = mortality_hi,
    Indicator = target_names[2]
  ),
  binomial_target(target_names[3], 2013, 3955, 21022, 100),
  binomial_target(target_names[4], 2010, 1310, 252940, 1e5),
  binomial_target(target_names[5], 2010, 561, 1301, 100),
  data.frame(
    Year      = 2019, Mean = 99.88,
    Low       = 99.88 - 1.96 * 21.66,
    High      = 99.88 + 1.96 * 21.66,
    Indicator = target_names[6]
  ),
  binomial_target(
    target_names[7], c(2000, 2010), c(447, 181), c(1340, 1310), 100
  ),
  data.frame(
    Year      = 2023,
    Mean      = acf_target * 100,
    Low       = pmax(0, (acf_target - 1.96 * acf_sd) * 100),
    High      = pmin(100, (acf_target + 1.96 * acf_sd) * 100),
    Indicator = target_names[8]
  )
)

message("Preparing model-fit plots from saved predictions...")
calibration_panels <- lapply(seq_along(target_names), function(i) {
  model    <- subset(model_data, Indicator == target_names[i])
  observed <- subset(target_data, Indicator == target_names[i])

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = model,
      ggplot2::aes(Year, ymin = Low, ymax = High),
      fill = "#2878B5", alpha = 0.18
    ) +
    ggplot2::geom_line(
      data = model,
      ggplot2::aes(
        Year, Mean,
        colour = "Model mean (95% uncertainty interval)"
      ),
      linewidth = 0.7
    ) +
    ggplot2::geom_errorbar(
      data = observed,
      ggplot2::aes(
        Year,
        ymin   = Low, ymax = High,
        colour = "Calibration targets"
      ),
      width = 0.5, linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = observed,
      ggplot2::aes(
        Year, Mean,
        colour = "Calibration targets"
      ),
      size = 2
    ) +
    ggplot2::scale_colour_manual(values = c(
      "Model mean (95% uncertainty interval)" = "#2878B5",
      "Calibration targets" = "#D9534F"
    )) +
    ggplot2::scale_x_continuous(
      breaks = seq(2000, 2025, 5), limits = c(2000, 2025)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      x     = "Year", y = target_units[i],
      title = target_names[i], colour = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )
})

# Display the shared legend inside the first calibration panel.
calibration_panels[[1]] <- calibration_panels[[1]] +
  ggplot2::theme(
    legend.position      = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction     = "vertical",
    legend.background    = ggplot2::element_blank(),
    legend.key.height    = grid::unit(0.35, "cm"),
    legend.text          = ggplot2::element_text(size = 9)
  )

calibration_panels[-1] <- lapply(
  calibration_panels[-1],
  function(panel) panel + ggplot2::theme(legend.position = "none")
)

p_calibration <- patchwork::wrap_plots(
  calibration_panels,
  ncol = 2
) +
  patchwork::plot_annotation(
    caption = paste(
      "L6 is an inferred mean under a lognormal assumption;",
      "its interval uses the working normal calibration approximation."
    )
  )

p_calibration <- p_calibration & figure_text_theme

print(p_calibration)

filename <- "Supplementary_Fig6_model_fit.png"

ggplot2::ggsave(
  filename = file.path(save_dir, filename),
  plot     = p_calibration,
  width    = 14,
  height   = 16,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

message("Saved: ", filename)

# ------------------------------------------------------------------------------
# 3.3 Plot model predictions: Figure 1
# ------------------------------------------------------------------------------

plot_fig1_pred <- function(model_pred) {
  end_year <- max(settings$pred_years)

  scenario_labels <- c(
    status_quo   = "Status quo (13% coverage)",
    scenario_50  = "Expansion with 50% coverage",
    scenario_80  = "Expansion with 80% coverage",
    scenario_100 = "Expansion with 100% coverage"
  )

  color_map <- c(
    "Status quo (13% coverage)"    = "black",
    "Expansion with 50% coverage"  = "purple",
    "Expansion with 80% coverage"  = "green",
    "Expansion with 100% coverage" = "orange"
  )

  plot_specs <- data.frame(
    Outcome = c(
      "inci_rate", "inci_num", "mort_rate", "death_num",
      "direct_tx", "ref_num", "ltfu_num", "delay_days"
    ),
    Scale = c(1, 1 / 1000, 1, 1 / 1000, 1e-6, 1e-6, 1e-6, 1),
    Ylab = c(
      "Incidence/100,000", "Number of cases (thousand)",
      "Mortality/100,000", "Number of deaths (thousand)",
      "Number (million)", "Number (million)", "Number (million)", "Days"
    ),
    Title = c(
      "a. Predicted active TB incidence",
      "b. Predicted new active TB cases",
      "c. Estimated TB mortality",
      "d. Estimated TB deaths",
      "e. First-visit diagnosis and treatment",
      "f. Referrals",
      "g. Losses to follow-up",
      "h. Mean health system delay"
    ),
    Ymax = c(100, NA, NA, 300, 1, 3, NA, NA),
    stringsAsFactors = FALSE
  )

  plot_one <- function(spec, ratios = NULL) {
    df <- model_pred[
      model_pred$Outcome == spec$Outcome &
        model_pred$Year >= 2015 & model_pred$Year <= end_year,
    ]

    df$Mean <- df$Mean * spec$Scale
    df$Low <- df$Low_2.5 * spec$Scale
    df$High <- df$High_97.5 * spec$Scale
    df$group <- unname(scenario_labels[df$Scenario])
    df$group <- factor(df$group, levels = names(color_map))

    # Use specified limits where supplied; otherwise allow data-based headroom.
    ymax <- spec$Ymax
    if (is.na(ymax)) {
      ymax <- max(pretty(c(0, max(df$Mean, na.rm = TRUE) * 1.08), n = 4))
    }
    y_breaks <- pretty(c(0, ymax), n = 4)
    y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= ymax]

    p <- ggplot2::ggplot(
      df, ggplot2::aes(x = Year, y = Mean, colour = group, fill = group)
    ) +
      ggplot2::geom_line(
        data      = df[df$group != "Status quo (13% coverage)", , drop = FALSE],
        linewidth = 1
      ) +
      ggplot2::geom_line(
        data      = df[df$group == "Status quo (13% coverage)", , drop = FALSE],
        linewidth = 1
      ) +
      ggplot2::scale_colour_manual(
        values = color_map,
        breaks = names(color_map)
      ) +
      ggplot2::scale_fill_manual(
        values = color_map,
        breaks = names(color_map)
      ) +
      ggplot2::scale_x_continuous(
        limits = c(2015, end_year),
        breaks = seq(2015, end_year, 5),
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        breaks = y_breaks,
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(ylim = c(0, ymax)) +
      ggplot2::labs(
        x      = "Year", y = spec$Ylab, title = spec$Title,
        colour = NULL, fill = NULL
      ) +
      ggplot2::theme_bw(
        base_size   = 22,
        base_family = "Times New Roman"
      ) +
      ggplot2::theme(
        legend.position = "none",
        plot.title      = ggplot2::element_text(size = 24, face = "bold"),
        axis.title      = ggplot2::element_text(size = 22),
        axis.text       = ggplot2::element_text(size = 21),
        plot.margin     = ggplot2::margin(10, 12, 10, 12, unit = "pt"),
        panel.grid      = ggplot2::element_blank()
      )

    if (!is.null(ratios)) {
      base_value <- df$Mean[
        df$Scenario == "status_quo" & df$Year == 2015
      ][1]

      target_points <- data.frame(
        Year  = c(2025, 2030, 2035),
        Value = base_value * (1 - ratios),
        Label = paste0(ratios * 100, "%\nreduction")
      )

      p <- p +
        ggplot2::geom_segment(
          data = target_points,
          ggplot2::aes(x = 2015, xend = Year, y = Value, yend = Value),
          inherit.aes = FALSE, colour = "red", linetype = "dashed",
          alpha       = 0.7, linewidth = 0.5
        ) +
        ggplot2::geom_segment(
          data = target_points,
          ggplot2::aes(x = Year, xend = Year, y = 0, yend = Value),
          inherit.aes = FALSE, colour = "red", linetype = "dashed",
          alpha       = 0.7, linewidth = 0.5
        ) +
        ggplot2::geom_point(
          data = target_points,
          ggplot2::aes(x = Year, y = Value),
          inherit.aes = FALSE, shape = 23, fill = "red",
          colour      = "red", size = 2.8, stroke = 1.2
        ) +
        ggplot2::geom_text(
          data = target_points,
          ggplot2::aes(
            x = Year - 3.7, y = Value + ymax / 14, label = Label
          ),
          inherit.aes = FALSE, colour = "red", size = 5.0,
          hjust = 0, lineheight = 0.9
        )
    }


    p
  }

  # Place the shared legend in the unused lower-right panel.
  legend_data <- data.frame(
    Label = names(color_map), Y = c(0.82, 0.67, 0.52, 0.37)
  )
  legend_panel <- ggplot2::ggplot(legend_data) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0.02, xend = 0.14, y = Y, yend = Y, colour = Label),
      linewidth = 1
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0.18, y = Y, label = Label),
      hjust = 0, size = 7.5, family = "Times New Roman"
    ) +
    ggplot2::scale_colour_manual(values = color_map, guide = "none") +
    ggplot2::annotate(
      "point", x = 0.08, y = 0.16, shape = 23,
      size = 3, fill = "red", colour = "red"
    ) +
    ggplot2::annotate(
      "text", x = 0.18, y = 0.16,
      label = "WHO End TB targets and\nmilestones compared with 2015",
      hjust = 0, size = 7.0, family = "Times New Roman"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    ggplot2::theme_void()

  panels <- list(
    plot_one(plot_specs[1, ], c(0.50, 0.80, 0.90)),
    plot_one(plot_specs[2, ]),
    plot_one(plot_specs[3, ]),
    plot_one(plot_specs[4, ], c(0.75, 0.90, 0.95)),
    plot_one(plot_specs[5, ]),
    plot_one(plot_specs[6, ]),
    plot_one(plot_specs[7, ]),
    plot_one(plot_specs[8, ]),
    legend_panel
  )

  patchwork::wrap_plots(panels, ncol = 3) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        plot.margin = ggplot2::margin(14, 18, 14, 18, unit = "pt")
      )
    )
}

p_fig1_pred <- plot_fig1_pred(psa_base$pred_data)

ggplot2::ggsave(
  file.path(save_dir, "Fig1_Projected_health_outcomes.png"),
  p_fig1_pred,
  width = 18, height = 16, units = "in",
  dpi   = 300, bg = "white"
)

# ------------------------------------------------------------------------------
# 3.4 Plot the cost-effectiveness plane: Figure 2
# ------------------------------------------------------------------------------

plot_fig2_cea_plane <- function(cea_baseline) {
  incr_cols <- c("incremental_50", "incremental_80", "incremental_100")

  points <- data.frame(
    Delta_E = c(0, cea_baseline["Total_DALYs", incr_cols] / 1e6),
    Delta_C = c(0, cea_baseline["Total_discounted_cost", incr_cols] / 1e9),
    ICER    = c(0, cea_baseline["ICER_USD_per_DALY_averted", incr_cols]),
    ICER2   = c(
      0, 0,
      cea_baseline["Sequential_ICER_USD_per_DALY_averted", "incremental_80"],
      cea_baseline["Sequential_ICER_USD_per_DALY_averted", "incremental_100"]
    ),
    Name = c(
      "Status quo (13% coverage)",
      "Expansion with 50% coverage",
      "Expansion with 80% coverage",
      "Expansion with 100% coverage"
    )
  )

  icer_fmt  <- format(round(points$ICER), big.mark = ",")
  icer2_fmt <- format(round(points$ICER2), big.mark = ",")

  points$Label <- with(points, c(
    Name[1],
    paste0(Name[2], ":\n$", icer_fmt[2], "/DALY averted (vs. status quo)"),
    paste0(
      Name[3], ":\n$", icer2_fmt[3],
      "/DALY averted (vs. expansion with 50% coverage);\n$",
      icer_fmt[3], "/DALY averted (vs. status quo)"
    ),
    paste0(
      Name[4], ":\n$", icer2_fmt[4],
      "/DALY averted (vs. expansion with 80% coverage);\n$",
      icer_fmt[4], "/DALY averted (vs. status quo)"
    )
  ))

  ggplot2::ggplot(points, ggplot2::aes(Delta_E, Delta_C)) +
    ggplot2::geom_abline(
      slope    = 13445 / 1000, intercept = 0,
      linetype = "dashed", colour = "#F2B5C0", linewidth = 0.7
    ) +
    ggplot2::geom_abline(
      slope    = 13445 * 3 / 1000, intercept = 0,
      linetype = "dashed", colour = "#F2B5C0", linewidth = 0.7
    ) +
    ggplot2::geom_line(colour = "#F2B5C0", linewidth = 0.8) +
    ggplot2::geom_point(colour = "#F2B5C0", size = 5) +
    ggplot2::geom_point(colour = "#792B37", size = 2) +
    ggplot2::geom_text(ggplot2::aes(label = Label),
      hjust = 0, vjust = 1.3,
      size  = 4.7, family = "Times New Roman"
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", colour = "grey50"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 20), breaks = seq(0, 20, 5),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-2, 25), breaks = seq(0, 25, 5),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = "DALYs averted (million)",
      y = "Incremental cost (US$, billion)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Times New Roman", size = 16),
      plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
      panel.grid  = ggplot2::element_blank()
    ) +
    ggplot2::annotate(
      "text",
      x      = 0.85, y = 15, label = "Three-time per capita GDP",
      colour = "black", angle = 87, size = 4.7, hjust = 0.1, vjust = -0.2,
      family = "Times New Roman"
    ) +
    ggplot2::annotate(
      "text",
      x      = 1.57, y = 15, label = "One-time per capita GDP",
      colour = "black", angle = 81, size = 4.7, hjust = 0.1, vjust = -0.2,
      family = "Times New Roman"
    )
}

p_fig2_cea <- plot_fig2_cea_plane(cea_baseline)

ggplot2::ggsave(
  file.path(
    save_dir,
    "Fig2_CEA_panel.png"
  ),
  p_fig2_cea,
  width = 10.7, height = 7, units = "in",
  dpi   = 300, bg = "white"
)

# ------------------------------------------------------------------------------
# 3.5 Plot the cost-effectiveness acceptability curve: Figure 3
# ------------------------------------------------------------------------------

plot_fig3_ceac <- function(ceac_data, fig_title = NULL, show_legend = TRUE,
                           legend_position = c(0.52, 0.80)) {
  colors <- c(
    "Status quo (13% coverage)"    = "black",
    "Expansion with 50% coverage"  = "purple",
    "Expansion with 80% coverage"  = "green",
    "Expansion with 100% coverage" = "orange"
  )

  ceac_data <- ceac_data[ceac_data$wtp <= 15000, ]
  ceac_long <- data.frame(
    wtp         = rep(ceac_data$wtp, 4),
    probability = unlist(ceac_data[c(
      "prob_status_quo", "prob_intervention1",
      "prob_intervention2", "prob_intervention3"
    )], use.names = FALSE),
    strategy = factor(
      rep(names(colors), each = nrow(ceac_data)),
      levels = names(colors)
    )
  )

  ggplot2::ggplot(
    ceac_long,
    ggplot2::aes(wtp, probability, colour = strategy)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(
      xintercept = 13445, linetype = "dashed", linewidth = 0.7
    ) +
    ggplot2::annotate(
      "text",
      x      = 13445 * 1.02, y = 0.18,
      label  = "One-time per capita GDP",
      angle  = 90, hjust = 0, vjust = 1.2, size = 5.5,
      family = "Times New Roman"
    ) +
    ggplot2::scale_colour_manual(values = colors, breaks = names(colors)) +
    ggplot2::scale_x_continuous(
      limits = c(0, 15000),
      breaks = seq(0, 15000, by = 3000),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(
      title  = fig_title,
      x      = "Willingness-to-pay (US$ per DALY averted)",
      y      = "Probability of being cost-effective",
      colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "Times New Roman", size = 16
      ),
      legend.position   = if (show_legend) "inside" else "none",
      legend.position.inside = legend_position,
      legend.direction  = "vertical",
      legend.key.height = grid::unit(0.65, "cm"),
      legend.text = ggplot2::element_text(
        family = "Times New Roman", size = 16
      ),
      legend.background = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(0.7, 1.0, 0.7, 0.7, unit = "cm"),
      panel.grid        = ggplot2::element_blank()
    )
}

p_fig3_ceac <- plot_fig3_ceac(psa_base$ceac_data)

ggplot2::ggsave(
  file.path(save_dir, "Fig3_CEAC.png"),
  p_fig3_ceac,
  width = 9, height = 6, units = "in",
  dpi   = 300, bg = "white"
)

# ------------------------------------------------------------------------------
# 3.6 Plot supplementary figures: Figures 7-9
# ------------------------------------------------------------------------------
# Map hospital coverage to visit coverage under the three assumptions.
plot_fig7_coverage <- function() {
  hosp_high <- c(0, cumsum(c(9.8, 23.7, 46.6, 19.8)))
  visit_high <- c(0, cumsum(c(64.7, 23.4, 8.3, 3.5)))
  hosp_low <- c(0, cumsum(c(19.8, 46.6, 23.7, 9.8)))
  visit_low <- c(0, cumsum(c(3.5, 8.3, 23.4, 64.7)))
  labels <- c(
    "Baseline: proportional visit coverage",
    "Worst case: visit coverage above hospital coverage",
    "Best case: visit coverage below hospital coverage"
  )
  curves <- rbind(
    data.frame(Hospital = c(0, 100), Visits = c(0, 100), Case = labels[1]),
    data.frame(Hospital = hosp_high, Visits = visit_high, Case = labels[2]),
    data.frame(Hospital = hosp_low, Visits = visit_low, Case = labels[3])
  )
  curves$Case <- factor(curves$Case, levels = labels)
  hosp_2020 <- 11.11
  markers <- data.frame(
    Hospital = hosp_2020,
    Visits = c(
      hosp_2020,
      approx(hosp_high, visit_high, hosp_2020)$y,
      approx(hosp_low, visit_low, hosp_2020)$y
    ),
    Case = factor(labels, levels = labels)
  )
  markers$Label <- sprintf(
    "Achieved %.1f%% TB service coverage in 2020", markers$Visits
  )
  ggplot2::ggplot(curves, ggplot2::aes(Hospital, Visits, colour = Case)) +
    ggplot2::geom_line(ggplot2::aes(linetype = Case), linewidth = 0.9) +
    ggplot2::geom_point(
      data = markers, shape = 21, fill = "white", colour = "black", size = 3
    ) +
    ggplot2::geom_text(
      data = markers, ggplot2::aes(label = Label),
      hjust = 0, nudge_x = 1, nudge_y = c(0, 0, -2),
      size = 3, family = "Times New Roman", show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c("blue", "red", "darkgreen")) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, 10)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 10)) +
    ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    ggplot2::labs(
      x = "General hospitals providing TB services (%)",
      y = "General-hospital visits covered by TB services (%)",
      colour = NULL, linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = 14, base_family = "Times New Roman") +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.99),
      legend.justification = c(0, 1), legend.direction = "vertical",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 9, face = "bold"),
      legend.key.height = grid::unit(0.4, "cm"),
      legend.key.width = grid::unit(1.2, "cm"),
      panel.grid.minor = ggplot2::element_blank()
    )
}
p_fig7_coverage <- plot_fig7_coverage()
ggplot2::ggsave(
  file.path(save_dir, "Supplementary_Fig7_hospital_visit_coverage.png"),
  p_fig7_coverage, width = 10, height = 6, dpi = 300, bg = "white"
)

# Compare the probability that each strategy is cost-effective across
# alternative assumptions about loss to follow-up during referral from GA.
p_fig_s8_ceac <- (
  plot_fig3_ceac(
    psa_sen60$ceac_data, "p_lf2 = 0.60",
    show_legend = TRUE, legend_position = c(0.56, 0.56)
  ) +
    plot_fig3_ceac(
      psa_sen45$ceac_data, "p_lf2 = 0.45",
      show_legend = FALSE
    ) +
    plot_fig3_ceac(
      psa_sen30$ceac_data, "p_lf2 = 0.30",
      show_legend = FALSE
    ) +
    plot_fig3_ceac(
      psa_sen15$ceac_data, "p_lf2 = 0.15",
      show_legend = FALSE
    )
) +
  patchwork::plot_layout(ncol = 2)

ggplot2::ggsave(
  file.path(save_dir, "Supplementary_Fig8_CEAC_loss_to_followup.png"),
  p_fig_s8_ceac,
  width = 15, height = 12, units = "in",
  dpi   = 300, bg = "white"
)

# Plot one-way sensitivity analysis of ICERs for the three expansion scenarios.
# Each tornado plot shows the 10 parameters with the largest ICER ranges.
plot_fig_s9_usa <- function(usa_data) {
  format_value <- function(x) {
    ifelse(x == 0 | x > 10,
           formatC(x, format = "f", big.mark = ",", digits = 0),
      ifelse(x < 0.015, sprintf("%.4f", x), sprintf("%.2f", x))
    )
  }

  process_icer <- function(df, icer_col, parameter) {
    icer  <- df[[icer_col]]
    valid <- which(is.finite(icer))
    if (!length(valid)) {
      return(NULL)
    }

    min_i <- valid[which.min(icer[valid])]
    max_i <- valid[which.max(icer[valid])]
    min_p <- df$par_value[min_i]
    max_p <- df$par_value[max_i]

    data.frame(
      parameter  = parameter, baseline = icer[df$value_n == 4][1],
      min_icer   = icer[min_i], max_icer = icer[max_i],
      range_icer = icer[max_i] - icer[min_i],
      min_par    = min_p, max_par = max_p,
      direction  = ifelse(max_p > min_p, "positive", "negative")
    )
  }

  get_top <- function(icer_col) {
    results <- do.call(rbind, lapply(names(usa_data), function(parameter) {
      process_icer(usa_data[[parameter]], icer_col, parameter)
    }))
    head(results[order(results$range_icer, decreasing = TRUE), ], 10)
  }

  tornado_plot <- function(top, title) {
    top$parameter <- factor(
      top$parameter,
      levels = top$parameter[order(top$range_icer)]
    )

    ggplot2::ggplot(top, ggplot2::aes(x = parameter)) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = as.numeric(parameter) - 0.3,
          xmax = as.numeric(parameter) + 0.3,
          ymin = pmin(baseline, max_icer), ymax = pmax(baseline, max_icer),
          fill = ifelse(direction == "positive",
            "High parameter value", "Low parameter value"
          )
        ),
        alpha = 0.65
      ) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = as.numeric(parameter) - 0.3,
          xmax = as.numeric(parameter) + 0.3,
          ymin = pmin(baseline, min_icer), ymax = pmax(baseline, min_icer),
          fill = ifelse(direction == "positive",
            "Low parameter value", "High parameter value"
          )
        ),
        alpha = 0.65
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = baseline), colour = "red", size = 1
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          y      = min_icer, label = format_value(min_par),
          colour = ifelse(direction == "positive", "lightblue", "lightcoral")
        ),
        hjust = 1.15, size = 3.2, family = "Times New Roman",
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          y      = max_icer, label = format_value(max_par),
          colour = ifelse(direction == "positive", "lightcoral", "lightblue")
        ),
        hjust = -0.15, size = 3.2, family = "Times New Roman",
        show.legend = FALSE
      ) +
      ggplot2::geom_hline(
        yintercept = top$baseline[1], linetype = "dashed", colour = "grey50"
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Low parameter value" = "lightblue",
          "High parameter value" = "lightcoral"
        ), name = NULL
      ) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_y_continuous(
        breaks = scales::breaks_pretty(n = 5),
        expand = c(0, 0)
      ) +
      ggplot2::coord_flip(ylim = c(0, max(4500, top$max_icer * 1.55))) +
      ggplot2::labs(
        title = title, x = "Parameter", y = "ICER (US$ per DALY averted)"
      ) +
      ggplot2::theme_bw(base_size = 12, base_family = "Times New Roman") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 11
        ),
        axis.title         = ggplot2::element_text(size = 10),
        axis.text          = ggplot2::element_text(size = 9),
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.08),
        legend.justification = c(1, 0),
        legend.direction = "vertical",
        legend.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 9),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey90")
      )
  }

  p1 <- tornado_plot(
    get_top("ICER1"),
    "Top 10 Sensitive Parameters (Expansion with 50% coverage vs. Status quo)"
  )

  p2 <- tornado_plot(
    get_top("ICER2"),
    "Top 10 Sensitive Parameters (Expansion with 80% coverage vs. Status quo)"
  )

  p3 <- tornado_plot(
    get_top("ICER3"),
    "Top 10 Sensitive Parameters (Expansion with 100% coverage vs. Status quo)"
  )

  p1 / p2 / p3
}

p_fig_s9_usa <- plot_fig_s9_usa(usa_data)

ggplot2::ggsave(
  file.path(save_dir, "Supplementary_Fig9_oneway_sensitivity.png"),
  p_fig_s9_usa,
  width = 7.5, height = 9.2, units = "in",
  dpi   = 300, bg = "white"
)

# ==============================================================================

cat("Total runtime:",
    round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes\n")
