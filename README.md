# aerosol-heron_island

## Aerosol size distribution analysis workflow from Heron Island (Great Barrier Reef)

This repository contains example R scripts to process, analyse, and visualise airborne aerosol data collected with a DJI Matrice 600 drone and a Brechtel Model 9405 Miniaturised Optical Particle Counter (mOPC). 

https://www.brechtel.com/product/miniaturized-optical-particle-counter/

The workflow demonstrates how to synchronise OPC data with drone flight logs, merge the datasets to derive altitude-resolved particle concentrations, and fit bimodal lognormal curves to parameterise aerosol size distributions.

---

## Overview

The project showcases:

- **Data ingestion**: Parsing example flight logs from the DJI Matrice 600 and example particle measurement files from the Brechtel mOPC.
- **Data merging**: Synchronising and merging datasets to associate aerosol observations with corresponding altitude measurements.
- **Aerosol size distribution modelling**: Fitting bimodal lognormal functions to the merged data to characterise aerosol populations.

---

## Scripts

### 1. `01_clean_dji_flightlogs.R`

- Reads and processes example flight log data from the DJI Matrice 600 exported from Airdata.com.
- Cleans altitude information for later merging with aerosol particle data.

### 2. `02_clean_mopc_merge_dji.R`

- Parses example data from the Brechtel Model 9405 mOPC.
- Calculates number size distributions based on particle counts across multiple size bins.
- Merges mOPC data with clean DJI flightlog data

### 3. `03_2mode_lognormal_fit.R`

- Merges the drone flight data with OPC measurements using timestamp alignment.
- Computes altitude for each aerosol observation.
- Fits a bimodal lognormal model to the resulting aerosol size distributions.
- Outputs summary plots optimised for storage (e.g., `.svg` or compressed `.png`) with clear labelling of observed and modelled distributions.

---

## Data

The repository includes:

- Sample DJI Matrice 600 flight log (CSV format).
- Sample Brechtel mOPC data file (CSV format).
- Merged dataset and derived plots.

> Note: The sample data has been anonymised and downsampled for demonstration purposes.

---

## Output

The final output consists of:

- **Merged datasets** with altitude-resolved particle size distribution data.
- **Parameterised aerosol models** for two flight events.
- **Publication-ready visualisations** saved in .svg lightweight formats suitable for sharing or versioning.

---

## Dependencies

### 📦 Required R Packages


install.packages(c(
  "data.table",
  "dplyr",
  "here",
  "lubridate",
  "minpack.lm",
  "pracma",
  "svglite",
  "tidyverse",
  "tools"
))

