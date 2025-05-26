################################################################
# Author: Christian Eckert
# Contact: c.eckert.10@student.scu.edu.au
# Project: Heron Island Atmospheric Aerosol Analysis
# Description: Fit and plot bimodal lognormal distributions of OPC data
# Last Updated: 12-01-2025
################################################################

# --- Environment Setup ---
rm(list = ls()); gc()
options(timeout = 10000)

# Load required libraries ---------------------------------------------------
library(data.table)
library(tidyverse)
library(minpack.lm)
library(pracma)
library(here)

# Define base project directory ---------------------------------------------
project_dir <- here::here()
setwd(project_dir)

# Define consistent input/output paths --------------------------------------
input_path <- file.path(project_dir, "data", "processed", "mopc_dji_clean", "total")
output_plot_path <- file.path(project_dir, "R", "plots", "logfit")
output_results_path <- file.path(project_dir, "data", "processed", "logfit")

# Create output folders if they don't exist ---------------------------------
dir.create(output_plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(output_results_path, recursive = TRUE, showWarnings = FALSE)

# --- User-defined Parameters ---
event_number <- 40  # Change this number for other events
event <- event_number %/% 10

# Locate matching OPC DJI file ----------------------------------------------
file_pattern <- paste0("total_", event_number, "_mOPC_")
opc_file_list <- list.files(path = input_path, pattern = file_pattern, full.names = TRUE)

if (length(opc_file_list) == 0) {
  stop("No OPC DJI file found for event ", event_number)
} else {
  opc_dji <- fread(opc_file_list[1])
}

opc_dji <- opc_dji %>% mutate(event = event)

# Prepare data --------------------------------------------------------------
filter_opc_dji_1 <- opc_dji %>% filter(BinID < 62)

locations <- c("island", "lagoon", "crest", "background")
hover_labels <- seq(5, 105, by = 5)

all_results <- list()

# Loop through each location and hover --------------------------------------
for (loc in locations) {
  data_loc <- filter_opc_dji_1 %>% filter(location == loc)
  
  for (hover in hover_labels) {
    data_hover <- data_loc %>%
      filter(hover_label == hover) %>%
      select(event, hover_label, dia, dN_dlogDp, location, tide)
    
    opc_avg <- data_hover %>%
      group_by(dia) %>%
      summarise(mean_dN_dlogDp = mean(dN_dlogDp, na.rm = TRUE), .groups = "drop")
    
    start_params <- list(
      N1 = max(opc_avg$mean_dN_dlogDp, na.rm = TRUE) * 0.8,
      Dg1 = log(190.1),
      sigma1 = 0.11,
      N2 = max(opc_avg$mean_dN_dlogDp, na.rm = TRUE) * 0.2,
      Dg2 = log(249.9),
      sigma2 = 0.126
    )
    
    fit <- tryCatch({
      nls(mean_dN_dlogDp ~ 
            N1 * dlnorm(dia, meanlog = Dg1, sdlog = sigma1) +
            N2 * dlnorm(dia, meanlog = Dg2, sdlog = sigma2),
          data = opc_avg,
          start = start_params,
          algorithm = "port",
          lower = c(0, log(160), 0.09, 0.001 * max(opc_avg$mean_dN_dlogDp), log(235), 0.09),
          upper = c(Inf, log(240), 0.128, Inf, log(530), 0.143))
    }, error = function(e) {
      message(paste("Skipping location:", loc, "hover:", hover, "– fit failed:", e$message))
      return(NULL)
    })
    
    if (is.null(fit)) next
    
    params <- coef(fit)
    modes <- data.frame(
      meanlog = c(params["Dg1"], params["Dg2"]),
      N = c(params["N1"], params["N2"]),
      sigma = c(params["sigma1"], params["sigma2"])
    ) %>% arrange(meanlog)
    
    mode_dia <- exp(modes$meanlog)
    mode_peak <- modes$N * dlnorm(mode_dia, meanlog = modes$meanlog, sdlog = modes$sigma)
    
    residuals_fit <- residuals(fit)
    rmse <- sqrt(mean(residuals_fit^2))
    ss_total <- sum((opc_avg$mean_dN_dlogDp - mean(opc_avg$mean_dN_dlogDp))^2)
    r_squared <- 1 - sum(residuals_fit^2) / ss_total
    
    result <- data_hover %>%
      select(event, hover_label, location) %>%
      mutate(
        m1_dia = mode_dia[1], m2_dia = mode_dia[2],
        m1_peak = mode_peak[1], m2_peak = mode_peak[2],
        m1_sigma = modes$sigma[1], m2_sigma = modes$sigma[2],
        m1_N = modes$N[1], m2_N = modes$N[2],
        m1_meanlog = modes$meanlog[1], m2_meanlog = modes$meanlog[2],
        R2 = r_squared, rmse = rmse
      )
    
    all_results[[length(all_results) + 1]] <- result
    
    dia_seq <- seq(min(opc_avg$dia), max(opc_avg$dia), length.out = 200)
    pred_total <- predict(fit, newdata = data.frame(dia = dia_seq))
    
    curve1 <- modes$N[1] * dlnorm(dia_seq, meanlog = modes$meanlog[1], sdlog = modes$sigma[1])
    curve2 <- modes$N[2] * dlnorm(dia_seq, meanlog = modes$meanlog[2], sdlog = modes$sigma[2])
    
    p <- ggplot() +
      geom_line(data = opc_avg, aes(x = dia, y = mean_dN_dlogDp), colour = "#440154", linewidth = 1.6) +
      geom_line(aes(x = dia_seq, y = curve1), colour = "#30678D", linewidth = 1, alpha = 0.8) +
      geom_line(aes(x = dia_seq, y = curve2), colour = "orange", linewidth = 1, alpha = 0.8) +
      geom_line(aes(x = dia_seq, y = pred_total), colour = "black", linetype = "dashed") +
      labs(title = paste("Event", event_number, "-", loc, "hover", hover),
           x = expression(Diameter~(nm)), y = expression(dN/dlogD[p])) +
      theme_minimal(base_size = 14)
    
    plot_file <- file.path(output_plot_path, paste0("event_", event_number, "_", loc, "_hover_", hover, ".svg"))
    ggsave(plot_file, plot = p, width = 8, height = 5, dpi = 96)
  }
}

# Combine and save results --------------------------------------------------
combined_results <- bind_rows(all_results)
output_file <- file.path(output_results_path, paste0("logfit_results_event_", event_number, ".csv"))
fwrite(combined_results, output_file)

message("Finished fitting and plotting for event ", event_number)

# End script---------------------------------------------------------------------
