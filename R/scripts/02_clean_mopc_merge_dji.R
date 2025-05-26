################################################################
# Author: Christian Eckert
# Contact: c.eckert.10@student.scu.edu.au
# Project: Heron Island Atmospheric Aerosol Analysis
# Description: Clean and process mOPC data and align with DJI flight logs
# Last Updated: 07-12-2024
################################################################



# Clear the environment -----------------------------------------------------
rm(list = ls())

# Load required libraries ---------------------------------------------------
library(data.table)
library(tidyverse)
library(lubridate)
library(here)
library(tools)

# Set working directory -----------------------------------------------------
project_dir <- here::here()
setwd(project_dir)

# Define input/output paths -------------------------------------------------
input_dji_path <- here("data", "processed", "dji_flightlogs_clean")
input_opc_path <- here("data", "raw", "mopc")
output_fine_path <- here("data", "processed", "mopc_dji_clean", "fine")
output_coarse_path <- here("data", "processed", "mopc_dji_clean", "coarse")
output_total_path <- here("data", "processed", "mopc_dji_clean", "total")

# Create output directories if they don't exist -----------------------------
dir.create(output_fine_path, recursive = TRUE, showWarnings = FALSE)
dir.create(output_coarse_path, recursive = TRUE, showWarnings = FALSE)
dir.create(output_total_path, recursive = TRUE, showWarnings = FALSE)

# Load DJI logs and extract tide condition from filenames -------------------
dji_files <- list.files(path = input_dji_path, pattern = "*.csv", full.names = TRUE)

dji_list <- lapply(dji_files, function(file) {
  dji <- fread(file)
  dji <- dji %>%
    rename(epoch = dji_epoch) %>%
    mutate(epoch = as.numeric(epoch))
    return(dji)
})

# Combine all DJI flight logs into one table --------------------------------
dji_all <- bind_rows(dji_list)

# Load mOPC files -----------------------------------------------------------
opc_files <- list.files(path = input_opc_path, pattern = "*.txt", full.names = TRUE)

# Define bin limits and compute reference table -----------------------------
bin_limits <- c(165.0, 167.8, 170.6, 173.4, 176.3, 179.3, 182.3, 185.3, 188.4, 191.6,
                194.8, 198.0, 201.4, 204.7, 208.2, 211.6, 215.2, 218.8, 222.4, 226.2,
                229.9, 233.8, 237.7, 241.7, 245.7, 249.8, 254.0, 258.3, 262.6, 267.0,
                271.5, 276.1, 280.8, 285.7, 290.7, 295.8, 301.2, 306.7, 312.6, 318.7,
                325.2, 332.1, 339.5, 347.5, 356.2, 365.6, 375.9, 387.4, 400.0, 414.2,
                430.0, 447.9, 468.1, 491.1, 517.2, 547.2, 581.5, 620.9, 666.4, 718.9,
                779.6, 849.9, 931.5, 1026.2, 1136.3, 1264.3, 1413.5, 1587.2, 1789.7,
                2025.9, 2301.4, 2623.0, 3000.0)

bin_mid <- sqrt(bin_limits[-length(bin_limits)] * bin_limits[-1])
vol_factor <- 1e-9 * (pi * (bin_mid^3)) / 6
dN_dlogD_factor <- 1 / log10(bin_limits[-1] / bin_limits[-length(bin_limits)])

bin_reference <- data.frame(
  BinID = seq_along(dN_dlogD_factor),
  dia = bin_mid,
  vol_factor = vol_factor,
  dN_dlogD_factor = dN_dlogD_factor
)

# Process each mOPC file ----------------------------------------------------
for (file in opc_files) {
  
  # Read mOPC file (skipping header) ---------------------------------------
  opc_raw <- fread(file, skip = 66, fill = TRUE) %>%
    rename(epoch = `#YY/MM/DD`, Time = `HR:MN:SC`) %>%
    mutate(
      Time = ymd_hms(paste(epoch, Time), tz = "UTC"),
      epoch = as.numeric(Time)
    )
  
  # Convert bin columns to numeric concentrations --------------------------
  bin_cols <- names(opc_raw)[18:89]
  opc_raw[, (bin_cols) := lapply(.SD, as.numeric), .SDcols = bin_cols]
  
  # Convert to concentration
  opc_raw[, (bin_cols) := lapply(.SD, function(x) x / (1000 / 60 * as.numeric(sample_flw))), .SDcols = bin_cols]
  
  # Add metadata and reshape
  opc_data <- opc_raw %>%
    mutate(
      total_conc2 = rowSums(across(all_of(bin_cols)), na.rm = TRUE),
      diff_tconc = total_conc - total_conc2
    ) %>%
    select(epoch, Time, total_conc, total_conc2, diff_tconc, sample_flw, sheath_flw, samp_rh, all_of(bin_cols))
  
  opc_long <- opc_data %>%
    pivot_longer(
      cols = starts_with("bin"),
      names_to = "BinID",
      names_prefix = "bin",
      values_to = "dN"
    ) %>%
    mutate(BinID = as.integer(BinID), dN = as.numeric(dN)) %>%
    left_join(bin_reference, by = "BinID")
  
  # Merge with DJI and tide info using 'epoch' -----------------------------
  opc_dji <- opc_long %>%
    left_join(dji_all, by = "epoch") %>%
    filter(!is.na(dji_lat)) %>%
    mutate(
      dN_dlogDp = dN * dN_dlogD_factor,
      volume = dN_dlogDp * vol_factor
    )
  
  # Save full data (includes tide info now) -------------------------------
  total_file <- file.path(output_total_path, paste0("total_", file_path_sans_ext(basename(file)), ".csv"))
  write_csv(opc_dji, total_file)
  
  # Split by particle size ------------------------------------------------
  opc_fine <- opc_dji %>% filter(BinID >= 3 & BinID <= 49)
  opc_coarse <- opc_dji %>% filter(BinID >= 49 & BinID <= 72)
  
  # Save fine and coarse data separately ----------------------------------
  fine_file <- file.path(output_fine_path, paste0("fine_", file_path_sans_ext(basename(file)), ".csv"))
  coarse_file <- file.path(output_coarse_path, paste0("coarse_", file_path_sans_ext(basename(file)), ".csv"))
  
  write_csv(opc_fine, fine_file)
  write_csv(opc_coarse, coarse_file)
}

# End of script -------------------------------------------------------------
