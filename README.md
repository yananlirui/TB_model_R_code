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
├── LICENSE           # MIT License for the model code and documentation
└── output/           # Saved results, Excel workbook and figures
```

Model inputs and calibration observations are embedded in the script.

## Requirements

The code has been checked with **R 4.6.1** on Windows 11 and Ubuntu 20.04.
The package versions used are also
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
No specialised hardware is required. The workflow adapts parallel execution to
the available logical cores; fewer cores increase computation time.

With the listed packages already installed, no additional installation is needed.
First-time package installation depends on internet speed and whether packages
must be built from source; this installation time has not been benchmarked.

## Running the script

**Before running, [download the saved RDS results](https://github.com/yananlirui/TB_model_R_code/releases/download/v1.1.0/TB_results_seed2026_20261112.rds)
and place the file in `output/`.** At approximately 300 MB, the RDS file exceeds
GitHub's regular repository file-size limit and is provided as a Release
attachment instead. If complete saved results are not found, the script reruns
the calibration and analyses, which can take several days. Downloading the RDS
allows the script to reuse the saved results and regenerate tables and figures.

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

For a short demonstration, download the RDS and run the script as above. This
regenerates the Excel workbook and Figures 1-3 and Supplementary Figures 3-9.
On the authors' Windows computer with R 4.6.1 and packages already installed,
this saved-results workflow took approximately 50 seconds. This excludes the RDS
download and initial package installation; timing varies by computer.
Supplementary Figure 2 is generated separately using `verify_resistance = TRUE`.

## Settings and workflow

Section **0.3** defines the seed, study periods and MCMC controls.
Defaults are seed **2026**, four chains of **250,000** iterations,
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

Section **2.1** runs all five calibrations within a shared CPU budget, concurrently
when sufficient cores are available and in smaller batches or serially otherwise.
Section **2.2** uses these fitted results for the three complete analyses;
the best and worst analyses do not repeat calibration.

Probabilistic sensitivity analyses use a general-hospital referral
loss-to-follow-up probability of 0.75 in the primary setting and 0.60, 0.45,
0.30 or 0.15 in alternative scenarios.

Equipment costs are annualised over the assumed equipment life. Equipment and
programme costs continue while expanded services operate; training costs are
charged once when each additional hospital starts providing services.

Parallel worker allocation is defined in section **1.4**. By default, the script
uses the detected number of logical cores minus two, with a minimum of one worker.
Reduce this budget manually on a shared server if necessary. A complete fresh run
can take several days; elapsed time depends on the hardware and available resources.

## Saved results and fresh runs

**[Download the saved RDS results](https://github.com/yananlirui/TB_model_R_code/releases/download/v1.1.0/TB_results_seed2026_20261112.rds)**
(`TB_results_seed2026_20261112.rds`, approximately 300 MB).
The file is distributed as a Release attachment and does not appear in the
repository's `output/` file list. Place the downloaded RDS file in `output/`
to regenerate tables and figures.

The script selects the latest dated file matching the current seed:

```text
output/TB_results_seed2026_YYYYMMDD.rds
```

Complete saved results are reused for the selected seed. Missing calibrations are scheduled together, subject to the available worker budget; analyses then use their corresponding calibrated parameters.
Generating tables and figures from complete saved results does not rewrite the RDS.

To recompute the entire workflow, use an empty `output/` directory after moving
existing results somewhere safe, or choose a seed without a matching saved RDS.

## Outputs

- **RDS:** five calibration results and three complete analyses, with diagnostics,
  posterior samples, model predictions and cost-effectiveness outputs.
  [Download the RDS file](https://github.com/yananlirui/TB_model_R_code/releases/download/v1.1.0/TB_results_seed2026_20261112.rds)
  and save it in `output/`.
- **Excel:** a dated workbook with seven sheets: `Table1`, `Table11`,
  `Table12_Best`, `Table13_Worst`, `Calibration_diagnostics`, `CEA_baseline`
  and `PSA_cea_final`.
- **Figures:** main Figures 1–3 show projected health outcomes, the
  cost-effectiveness plane and cost-effectiveness acceptability curves.
  Supplementary Figures 2–9 show the rifampicin-resistance trend, MCMC traces,
  autocorrelation, posterior distributions, calibration fit, hospital/visit
  coverage relationships, loss-to-follow-up sensitivity analyses and one-way
  sensitivity analyses.

Section **0.4.4** contains the resistance-trend fit and Supplementary Figure 2.
Set `verify_resistance = TRUE` to rerun the `nlminb` fit and generate
`Supplementary_Fig2_rifampicin_resistance.png`.
Section **3.2** generates Supplementary Figures 3–6 directly.
Set `verify_environment = TRUE` to print the R and package versions.

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

## License

The model source code and accompanying documentation are released under the
[MIT License](LICENSE). R and third-party packages retain their own licenses.
