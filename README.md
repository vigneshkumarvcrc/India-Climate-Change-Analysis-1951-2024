# India Climate Change Analysis, 1951–2024

This repository contains the original R scripts and computational
workflows used for the climate analysis presented in the manuscript.

## Contents

- `climate_analysis.R`
  Main analysis workflow for annual climate trends, anomalies,
  spatial metrics, risk indices, decadal summaries and maps.

- `audit_climate_inputs.R`
  Quality-control audit of input raster and vector datasets.

- `decadal_tables.R`
  Generation of decadal climate summary tables.

- `create_combined_risk_map.R`
  Generation of the combined climate-hazard index and map.

- `summarize_results.R`
  Generation of spatial summary statistics and state/district summaries.

## Software

The analysis was performed using R.

Required packages include:

- terra
- ggplot2
- dplyr
- scales
- viridisLite

## Input data

The analysis uses annual gridded temperature and precipitation
datasets covering 1951–2024, together with India, state and district
boundary datasets.

The original input datasets are not redistributed in this repository.
Their sources and instructions for obtaining the data are provided
in the manuscript and/or supplementary documentation.

## Running the analysis

1. Download or obtain the required input datasets.
2. Set the input and output directory paths in `climate_analysis.R`.
3. Run `audit_climate_inputs.R` to check the input datasets.
4. Run `climate_analysis.R`.
5. Run `decadal_tables.R` if required to regenerate the decadal tables.
6. Run `create_combined_risk_map.R` to generate the combined risk map.
7. Run `summarize_results.R` to generate the spatial summary tables.

## Reproducibility

The scripts provided here correspond to the computational workflow
used in the manuscript. The repository is intended to allow reviewers
to inspect the original analysis code and computational procedures.

## Code availability

The repository is publicly accessible and does not require personal
information, login credentials, or passwords to access the code.
