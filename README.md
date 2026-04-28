# Quantifying Renewable Energy Resource Uncertainty: A Nonparametric and Copula-Based Approach

**Nattirat Mayer**  
Bureau of Theoretical and Applied Economics (BETA UMR 7522)  
University of Strasbourg, CNRS  
Contact: npromwang@unistra.fr

---

## Overview

This repository contains the R code used in the paper *Quantifying Renewable Energy Resource Uncertainty: A Nonparametric and Copula-Based Approach*. The paper proposes a methodology for uncertainty quantification that combines kernel density estimation (KDE) for marginal distributions with Joe-Clayton copulas for multivariate joint dependence, applied to 23 years of hourly meteorological data at Mont Blanc, France (latitude 45°N, longitude 6.75°E).

---

## Repository Structure

```
├── 01_x_s_summary.R          # Solar irradiance: summary statistics and time series (Figures 1a, 2; Table 2)
├── 02_x_w_summary.R          # Wind speed: summary statistics and time series (Figures 1b, 3; Table 2)
├── 03_x_r_summary.R          # Total precipitation: summary statistics and time series (Figures 1c, 4; Table 2)
├── 04_x_s_marginal.R         # Solar irradiance: KDE vs lognormal, Monte Carlo, CDF (Figures 5)
├── 05_x_w_marginal.R         # Wind speed: KDE vs Rayleigh, Monte Carlo, CDF (Figures 6)
├── 06_x_r_marginal.R         # Total precipitation: KDE vs Gumbel, Monte Carlo, mixture CDF (Figures 7)
├── 07_copula_daytime.R       # Daytime copula fitting, tail dependence, bootstrap uncertainty (Figures 8 (a-c), 9(a), 10; Table 3)
├── 08_copula_nighttime.R     # Nighttime copula fitting, tail dependence, bootstrap uncertainty (Figure 9(b), 11; Table 4)
└── README.md
```

Scripts should be run in order. Each script reads data from the paths specified in the configuration section at the top and writes outputs to the working directory.

---

## Data

Meteorological data were downloaded manually from the Copernicus Climate Data Store (CDS):

- **Dataset:** ERA5 hourly data on single levels  
- **URL:** https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels  
- **Variables:**
  - Surface solar radiation downwards (`ssrd`) → solar irradiance X_s
  - 100m wind speed (`var_100_metre_wind_speed`) → wind speed X_w
  - Total precipitation (`tp`) → precipitation X_r
- **Location:** latitude 45°N, longitude 6.75°E (Mont Blanc, France)
- **Period:** 2000–2022 (hourly)
- **Format:** NetCDF, converted to CSV prior to analysis

The raw NetCDF files are not hosted in this repository due to file size. The processed CSV files used in the analysis are available at: [insert Zenodo or Figshare DOI here].

A Data Availability Statement is included in the paper.

---

## Requirements

All scripts are written in R. The following packages are required:

```r
install.packages(c(
  "dplyr", "readr", "lubridate", "ggplot2", "gridExtra",
  "scales", "zoo", "kde1d", "MASS", "VineCopula",
  "purrr", "stringr", "tidyr", "future", "future.apply",
  "plot3D"
))
```

Package versions used in the analysis are listed below. Results may differ with substantially different versions.

| Package       | Version |
|---------------|---------|
| R             | 4.x.x   |
| VineCopula    |         |
| kde1d         | 0.3.4   |
| future        |         |
| future.apply  |         |
| dplyr         |         |
| ggplot2       |         |

*Full session info available on request via `sessionInfo()` output.*

---

## Reproducing the Analysis

### Step 1 — Set data paths

At the top of each script, update the `folder_path` (scripts 01–06) or `path_r`, `path_s`, `path_w` (scripts 07–08) variables to point to the location of the processed CSV files on your machine.

### Step 2 — Run marginal distribution scripts (04–06)

These scripts fit KDE and parametric distributions (lognormal, Rayleigh, Gumbel) to each meteorological variable, run Monte Carlo simulations to compare goodness-of-fit, and export CDF (probability integral transform) CSV files. The exported CDF files are the inputs for scripts 07–08.

Scripts 04 and 05 export CDF files to:
- `x_s_CDF_daytime/`, `x_s_CDF_nighttime/`
- `x_w_CDF/`

Script 06 exports to:
- `x_r_CDF/`, `x_r_CDF_exports/`

### Step 3 — Run copula scripts (07–08)

Update the `path_r`, `path_s`, `path_w` variables in scripts 07 and 08 to point to the CDF export folders from Step 2. These scripts fit Joe-Clayton (BB7) copulas for each season-year combination, compute tail dependence coefficients, run bootstrap resampling of Kendall's tau, and produce the final figures and result CSVs.

Output CSVs:
- `copula_summary_bb7_daytime.csv`
- `copula_uncertainty_bb7_daytime.csv`
- `copula_summary_bb7_nighttime.csv`
- `copula_uncertainty_bb7_nighttime.csv`

---

## AI Use Disclosure

The author used Claude (claude.ai, Anthropic) to assist with R code writing and English language editing during the preparation of this manuscript. The author takes full responsibility for the content, methodology, and results presented.

---

## Citation

If you use this code, please cite the paper:

> Mayer, N. (forthcoming). *Quantifying Renewable Energy Resource Uncertainty: A Nonparametric and Copula-Based Approach*. [Journal name].

---

## Acknowledgements

The author thanks Bertrand Koebel for his valuable comments and feedback. This gratitude is extended to Nathalie Picard, Ralf Brüggemann, Michael Kupper, Marcus Rockel, and the participants at the Action versus Inaction Facing Climate Change conference 2025 in Lausanne, Switzerland, for their helpful suggestions and discussions. All remaining errors are the author's own.

---

## License

This code is shared for reproducibility purposes in accordance with the data and code availability policy of the Journal of the Royal Statistical Society Series C.
