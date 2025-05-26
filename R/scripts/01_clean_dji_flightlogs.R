################################################################
# Author: Christian Eckert
# Contact: c.eckert.10@student.scu.edu.au
# Project: Heron Island Atmospheric Aerosol Analysis
# Description: Clean DJI drone flight logs exported from AirData.com.
# Includes tide information extracted from filenames.
# Last Updated: 05-12-2024
################################################################

# Clear environment and set options
rm(list = ls())
options(digits = 10)

# Load libraries
library(data.table)
library(dplyr)
library(lubridate)
library(here)

# Set working directory -----------------------------------------------------
project_dir <- here::here()  # project root
setwd(project_dir)

# Define input/output directories
input_dir <- here("data", "raw", "dji_flightlogs")
output_dir <- here("data", "processed", "dji_flightlogs_clean")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Get list of flight log CSVs
file_list <- list.files(path = input_dir, pattern = "\\.csv$", full.names = TRUE)

# Define column selection and renaming
col_select <- c("datetime(utc)", "latitude", "longitude", 
                "height_above_takeoff(meters)", 
                "compass_heading(degrees)", "pitch(degrees)", "roll(degrees)")
col_names <- c("dji_datetime", "dji_lat", "dji_long", "dji_alt", 
               "dji_compass", "dji_pitch", "dji_roll")

# Define sampling locations
location_data <- data.frame(
  Location = c("island", "lagoon", "crest", "background"),
  Longitude = c(151.917885, 151.918148, 151.918625, 151.918813),
  Latitude  = c(-23.444112, -23.446360, -23.450425, -23.452012)
)

# Tolerance for location matching
tolerance <- 0.00003

# Altitude change threshold
altitude_threshold <- 0.42

# Process each flight log
for (file in file_list) {
  
  # Read and rename data
  dji <- fread(file, select = col_select, col.names = col_names,
               colClasses = list(character = 2))
  
  # Extract metadata from filename
  filename <- basename(file)
  tide_match <- regmatches(filename, regexpr("High|Low", filename, ignore.case = TRUE))
  tide <- ifelse(length(tide_match) > 0, tide_match, NA)
  
  # Convert datetime to POSIXct and epoch
  dji$dji_datetime <- parse_date_time(dji$dji_datetime, orders = c("dmy HMS", "ymd HMS"), tz = "UTC")
  dji$dji_epoch <- as.numeric(dji$dji_datetime)
  
  # Round compass and location data
  dji <- dji %>%
    mutate(
      dji_compass_rounded = ifelse(round(dji_compass) %in% c(360, 359, 0, 1), 0, round(dji_compass)),
      dji_deviation = dji_compass - round(dji_compass),
      dji_lat = round(dji_lat, 6),
      dji_long = round(dji_long, 6)
    )
  
  # Aggregate by epoch second
  dji_agg <- dji %>%
    group_by(dji_epoch) %>%
    summarise(
      dji_datetime = min(dji_datetime, na.rm = TRUE),
      dji_lat = mean(dji_lat, na.rm = TRUE),
      dji_long = mean(dji_long, na.rm = TRUE),
      dji_alt = mean(dji_alt, na.rm = TRUE),
      dji_pitch = mean(dji_pitch, na.rm = TRUE),
      dji_roll = mean(dji_roll, na.rm = TRUE),
      dji_compass_rounded = round(mean(dji_compass_rounded, na.rm = TRUE)),
      dji_deviation = mean(dji_deviation, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Match sampling locations
  dji_agg <- dji_agg %>%
    rowwise() %>%
    mutate(
      location = location_data$Location[
        which.min((dji_long - location_data$Longitude)^2 +
                    (dji_lat - location_data$Latitude)^2)
      ],
      match_error = sqrt((dji_long - location_data$Longitude[which.min((dji_long - location_data$Longitude)^2 +
                                                                         (dji_lat - location_data$Latitude)^2)])^2 +
                           (dji_lat - location_data$Latitude[which.min((dji_long - location_data$Longitude)^2 +
                                                                         (dji_lat - location_data$Latitude)^2)])^2),
      location = ifelse(match_error < tolerance, location, NA)
    ) %>%
    ungroup()
  
  # Filter: Valid locations, compass near 0, and sufficient duration (11+ seconds)
  dji_agg <- dji_agg %>%
    mutate(valid = location %in% location_data$Location,
           group = cumsum(c(TRUE, diff(valid) != 0))) %>%
    group_by(group) %>%
    mutate(consec_count = sum(valid)) %>%
    ungroup() %>%
    filter(consec_count >= 11 & valid & dji_compass_rounded < 1) %>%
    select(-valid, -group, -consec_count)
  
  # Detect descent and hover transitions
  dji_clean <- dji_agg %>%
    mutate(
      alt_diff = abs(dji_alt - lag(dji_alt, default = first(dji_alt))),
      is_descending = dji_alt < lag(dji_alt, default = first(dji_alt)) &
        alt_diff > altitude_threshold,
      cut_off = cumsum(is_descending)
    ) %>%
    filter(cut_off == 0 | row_number() == 1) %>%
    select(-cut_off, -is_descending) %>%
    mutate(
      alt_diff = abs(dji_alt - lag(dji_alt, default = first(dji_alt))),
      is_moving = alt_diff > altitude_threshold,
      transition = lag(is_moving, default = FALSE) == FALSE & is_moving == TRUE,
      hover_label = 0
    )
  
  # Label hover blocks (15 rows before transition)
  for (i in which(dji_clean$transition)) {
    start_row <- max(1, i - 15)
    avg_alt <- round(mean(dji_clean$dji_alt[start_row:(i - 1)], na.rm = TRUE) / 5) * 5
    dji_clean$hover_label[start_row:(i - 1)] <- avg_alt
  }
  
  # Label last 15 rows
  last_row <- nrow(dji_clean)
  start_row <- max(1, last_row - 14)
  avg_alt_last <- round(mean(dji_clean$dji_alt[start_row:last_row], na.rm = TRUE) / 5) * 5
  dji_clean$hover_label[start_row:last_row] <- avg_alt_last
  
  # Final clean output, include tide info
  final_clean <- dji_clean %>%
    filter(hover_label != 0) %>%
    select(-alt_diff, -is_moving, -transition) %>%
    mutate(tide = tide)
  
  # Export cleaned CSV
  out_name <- here(output_dir, paste0(tools::file_path_sans_ext(filename), "_clean.csv"))
  fwrite(final_clean, out_name)
}

# End script---------------------------------------------------------------------