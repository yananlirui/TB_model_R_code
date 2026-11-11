# TB transmission and cost-effectiveness model

R code accompanying Li et al., *Evaluating the impact of expanding tuberculosis
diagnostic and treatment availability of general hospitals in China*.

## Overview

The model represents pulmonary tuberculosis (TB) across all ages in China,
including infection, subclinical and clinical disease, care-seeking, referral,
treatment, recovery, relapse and reinfection. Extrapulmonary TB is not modelled.
TB disease incidence and case notifications are separate model outcomes:
incidence counts episodes entering active disease, whereas notifications count
flows into treatment through passive care-seeking or active case-finding.

The evaluation covers 2026–2035. Expansion scenarios reach 50%, 80% or 100% of
general hospitals providing TB services by 2030 and remain constant thereafter.
Under the status quo, hospital coverage follows its existing trend until it
reaches approximately 14.5%.

## Files

```text
TB_model_R_code/
├── TB_model_main.R    # Model, calibration, analyses and output generation
├── README.md
└── output/           # Saved results, Excel workbook and figures
```

Model inputs and calibration observations are embedded in the script.

## Requirements

The code has been checked with **R 4.6.1**. The package versions used are also
listed in section 0.1 of the script:

| Package | Version | Package | Version |
|---|---|---|---|
| coda | 0.19.4.1 | cowplot | 1.2.0 |
| deSolve | 1.42 | dplyr | 1.2.1 |
| FME | 1.3.6.4 | ggplot2 | 4.0.3 |
| grid | 4.6.1 | gridExtra | 2.3.1 |
| gtools | 3.9.5 | lubridate | 1.9.5 |
| parallel | 4.6.1 | patchwork | 1.3.2 |
| pbapply | 1.7.5 | pracma | 2.4.6 |
| rstudioapi | 0.19.0 | triangle | 1.1.0 |
| writexl | 2.0.1 | | |

`grid` and `parallel` are included with R. The script attempts to install missing
packages, so an initial installation requires internet access and a configured
CRAN repository.
Times New Roman should be installed for figures that request this font.
Linux servers also need a working headless graphics device, such as Cairo.

## Running the script

From the directory containing the script:

```bash
Rscript TB_model_main.R
```

Alternatively, open the saved script in RStudio and source it, or run:

```r
source("TB_model_main.R")
```

Use a fresh R session: the script clears objects in the global environment.
Output paths are resolved relative to the script, rather than the working
directory. Close the output Excel workbook before regenerating it.

## Settings and workflow

Section **0.3** defines the seed, study periods, MCMC controls and calibration
plotting option. Defaults are seed **2026**, four chains of **250,000** iterations,
50% burn-in and retention of every 50th post-burn-in draw for simulations. This
yields 10,000 joint posterior parameter sets per calibration. Parameter summaries
and convergence diagnostics use all retained post-burn-in draws before thinning.

The full workflow contains:

1. Three calibrations assuming that 0%, 50% or 100% of individuals with cleared
   infection remain IGRA-positive. The 0% assumption is used in the primary analysis.
2. Two additional calibrations under alternative relationships between hospital
   coverage and the proportion of general-hospital visits covered by TB services.
3. Three complete analyses: primary, best and worst care-engagement assumptions.
   Each includes baseline cost-effectiveness analysis, one-way sensitivity
   analysis and probabilistic sensitivity analyses. The 50% and 100% IGRA
   calibrations also provide baseline cost-effectiveness comparisons.

Probabilistic sensitivity analyses use a general-hospital referral
loss-to-follow-up probability of 0.75 in the primary setting and 0.60, 0.45,
0.30 or 0.15 in alternative scenarios.

Equipment costs are annualised over the assumed equipment life. Equipment and
programme costs continue while expanded services operate; training costs are
charged once when each additional hospital starts providing services.

Parallel worker allocation is defined in section **1.4**. By default, the script
uses the detected number of logical cores minus two, with a minimum of one worker.
Reduce this budget manually on a shared server if necessary. A complete fresh run
can take many hours; elapsed time depends on the hardware and available resources.

## Saved results and fresh runs

Saved RDS results are available from the
[GitHub Releases page](https://github.com/yananlirui/TB_model_R_code/releases).
Place the downloaded RDS file in `output/` to regenerate tables and figures.

The script selects the latest dated file matching the current seed:

```text
output/TB_results_seed2026_YYYYMMDD.rds
```

If the three IGRA calibrations and primary analysis are present for that seed,
they are reused. Missing best/worst calibration and analysis pairs are computed
separately.
Generating tables and figures from complete saved results does not rewrite the RDS.

To recompute the entire workflow, use an empty `output/` directory after moving
existing results somewhere safe, or choose a seed without a matching saved RDS.

## Outputs

- **RDS:** five calibration results and three complete analyses, with diagnostics,
  posterior samples, model predictions and cost-effectiveness outputs.
- **Excel:** a dated workbook with seven sheets: `Table1`, `Table11`,
  `Table12_Best`, `Table13_Worst`, `Calibration_diagnostics`, `CEA_baseline`
  and `PSA_cea_final`.
- **Main figures:** `Fig1_Projected_health_outcomes.png`, `Fig2_CEA_panel.png`
  and `Fig3_CEAC.png`.
- **Supplementary figures 3–9:** MCMC traces, autocorrelation, posterior
  distributions, calibration fit, hospital/visit coverage relationships,
  loss-to-follow-up sensitivity analyses and one-way sensitivity analyses.

Set `plot_calib = FALSE` to skip the calibration diagnostic figures. Set
`verify_resistance = TRUE` to rerun the `nlminb` resistance-trend fit and generate
`Supplementary_Fig2_rifampicin_resistance.png`. Set `verify_environment = TRUE`
to print the R and package versions. These switches are defined in the script.

Table 1 and its best/worst counterparts use deterministic baseline point estimates
for costs, DALYs and cost-effectiveness outcomes. Other point estimates and all
displayed uncertainty intervals use the corresponding probabilistic summaries.
Intervals are the 2.5th and 97.5th percentiles. Costs and DALYs are discounted at
3% annually in the base case; costs are expressed in 2024 US dollars.

Figures use stable filenames and are replaced when regenerated. The Excel filename
includes the export date and seed. Intermediate checkpoints are stored under
`output/process/`.

## Code organisation

| Section | Contents |
|---|---|
| S0 | Environment, settings, parameters and resistance-trend verification |
| S1 | Population, model, calibration, CEA and sensitivity-analysis functions |
| S2 | Calibration and analysis execution, diagnostics and checkpoints |
| S3 | Excel tables and figures |
