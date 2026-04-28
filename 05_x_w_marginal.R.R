# =============================================================
# Wind Speed (X_w): KDE vs Rayleigh, Monte Carlo, CDF
# Paper: Quantifying Renewable Energy Resource Uncertainty
# =============================================================
rm(list=ls())
# --- Libraries -----------------------------------------------
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(scales)
library(kde1d)
library(zoo)

# --- 1. Data loading -----------------------------------------

folder_path <- "D:/PhD/thesis/chapter_1/datasets/x_w/x_w_all"

years      <- 2000:2022
file_names <- paste0("x_w_fra_", years, ".csv")
file_paths <- file.path(folder_path, file_names)

filtered_data <- list()
for (file in file_paths) {
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    subset <- data %>%
      filter(latitude == 45, longitude == 6.75) %>%   # Mont Blanc
      dplyr::select(time, x_w = var_100_metre_wind_speed) %>%
      mutate(time = ymd_hms(time))
    filtered_data[[file]] <- subset
  } else {
    warning(paste("File does not exist:", file))
  }
}

# Check number of valid files loaded
length(filtered_data)
filtered_lengths <- sapply(filtered_data, nrow)
summary(filtered_lengths)

# --- 2. Data preprocessing -----------------------------------

x_w_df <- bind_rows(filtered_data) %>%
  arrange(time) %>%
  mutate(
    x_w    = ifelse(x_w < 0, 0, x_w),
    roll30 = rollmean(x_w, k = 24 * 30, fill = NA),   # 30-day rolling mean
    year   = year(time),
    month  = month(time),
    season = case_when(
      month %in% c(12, 1, 2)  ~ "winter",
      month %in% c(3, 4, 5)   ~ "spring",
      month %in% c(6, 7, 8)   ~ "summer",
      month %in% c(9, 10, 11) ~ "autumn"
    )
  )

# --- 3. Summary statistics (Table 2) -------------------------

summary_stats <- x_w_df %>%
  summarise(
    mean = mean(x_w, na.rm = TRUE),
    sd   = sd(x_w, na.rm = TRUE),
    p25  = quantile(x_w, 0.25, na.rm = TRUE),
    p50  = quantile(x_w, 0.50, na.rm = TRUE),
    p75  = quantile(x_w, 0.75, na.rm = TRUE)
  )

print(summary_stats)

# --- 4. Time series plot with rolling mean (Figure 1b) -------

x_w_df_sampled <- x_w_df %>%
  slice(seq(1, n(), by = 24 * 30))

p_raw <- ggplot(x_w_df_sampled, aes(x = time)) +
  geom_line(
    aes(y = x_w, color = "Raw data", linetype = "Raw data"),
    size = 0.7
  ) +
  geom_line(
    aes(y = roll30, color = "30-day rolling mean", linetype = "30-day rolling mean"),
    size = 1
  ) +
  labs(
    x        = "Year",
    y        = expression(X[w]),
    color    = "Legend",
    linetype = "Legend"
  ) +
  scale_color_manual(values = c(
    "Raw data"             = "black",
    "30-day rolling mean"  = "blue"
  )) +
  scale_linetype_manual(values = c(
    "Raw data"             = "solid",
    "30-day rolling mean"  = "solid"
  )) +
  scale_x_datetime(labels = date_format("%Y"), date_breaks = "2 years") +
  theme_minimal()

p_raw

# --- 5. Seasonal PDF plots (Figure 3) ------------------------
x_w_2000 <- x_w_df %>% filter(year(time) == 2000)

winter_wind <- x_w_2000 %>%
  filter(
    time >= ymd_hms("2000-12-21 00:00:00"),
    time <= ymd_hms("2001-03-20 23:59:59")
  )

summer_wind <- x_w_2000 %>%
  filter(
    time >= ymd_hms("2000-06-21 00:00:00"),
    time <= ymd_hms("2000-09-20 23:59:59")
  )

p_winter <- ggplot(winter_wind, aes(x = x_w)) +
  geom_density(fill = "grey", alpha = 0.3) +
  annotate("segment",
           x = 0, xend = max(winter_wind$x_w, na.rm = TRUE),
           y = 0, yend = 0,
           color = "black", size = 0.6) +
  coord_cartesian(xlim = c(0, max(winter_wind$x_w, na.rm = TRUE))) +
  labs(x = expression(X[w]), y = "Density", title = "Winter") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

p_summer <- ggplot(summer_wind, aes(x = x_w)) +
  geom_density(fill = "grey", alpha = 0.3) +
  annotate("segment",
           x = 0, xend = max(summer_wind$x_w, na.rm = TRUE),
           y = 0, yend = 0,
           color = "black", size = 0.6) +
  coord_cartesian(xlim = c(0, max(summer_wind$x_w, na.rm = TRUE))) +
  labs(x = expression(X[w]), y = "Density", title = "Summer") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

grid.arrange(p_winter, p_summer, ncol = 2)

# --- 6. Helper functions -------------------------------------

# Rayleigh PDF
d_rayleigh <- function(x, sigma) {
  (x / sigma^2) * exp(-x^2 / (2 * sigma^2))
}

# Rayleigh CDF
p_rayleigh <- function(x, sigma) {
  1 - exp(-x^2 / (2 * sigma^2))
}

# Rayleigh sampler
r_rayleigh <- function(n, sigma) {
  sigma * sqrt(-2 * log(runif(n)))
}

# Empirical CDF evaluated on a grid
ecdf_on_grid <- function(x, grid) {
  ecdf(x)(grid)
}

# --- 7. Side-by-side PDF plot: KDE vs Rayleigh (Figure 6) ----

plot_wind_side_by_side_mc <- function(data_vector, kde_fit, sigma_hat, winner) {
  data_vector <- data_vector[!is.na(data_vector) & data_vector > 0]
  if (length(data_vector) < 5) return(NULL)
  
  emp_density <- density(data_vector, from = 0)
  emp_df      <- data.frame(x = emp_density$x, y = emp_density$y)
  
  kde_df <- data.frame(x = kde_fit$grid_points, y = kde_fit$values)
  
  xgrid       <- seq(0, max(kde_df$x), length.out = 300)
  rayleigh_df <- data.frame(x = xgrid, y = d_rayleigh(xgrid, sigma = sigma_hat))
  
  # Plot 1: KDE vs Empirical
  p1 <- ggplot() +
    geom_line(data = kde_df,
              aes(x = x, y = y, color = "KDE"), size = 0.8) +
    geom_line(data = emp_df,
              aes(x = x, y = y, color = "Empirical"),
              linetype = "dashed", size = 0.5) +
    labs(x = expression(X[w]), y = "Density", color = "Legend") +
    scale_color_manual(values = c("KDE" = "darkgreen", "Empirical" = "black")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_line(color = "grey90", size = 0.25)
    ) +
    xlim(0, max(kde_df$x))
  
  # Plot 2: Rayleigh vs Empirical
  p2 <- ggplot() +
    geom_line(data = rayleigh_df,
              aes(x = x, y = y, color = "Rayleigh"), size = 0.8) +
    geom_line(data = emp_df,
              aes(x = x, y = y, color = "Empirical"),
              linetype = "dashed", size = 0.5) +
    labs(x = expression(X[w]), y = "Density", color = "Legend") +
    scale_color_manual(values = c("Rayleigh" = "cornflowerblue", "Empirical" = "black")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_line(color = "grey90", size = 0.25)
    ) +
    xlim(0, max(kde_df$x))
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

# --- 8. CDF plot function-------------------------
# Note: unlike precipitation, wind speed has no zero-inflation issue.
# The CDF is plotted over the full range of observed values including zero.

plot_wind_cdf_mc <- function(data_vector, kde_fit, sigma_hat, winner) {
  data_vector <- data_vector[!is.na(data_vector)]   # keep zeros
  if (length(data_vector) < 5) return(NULL)
  
  xplot <- sort(data_vector)
  Fvals <- if (winner == "KDE") pkde1d(xplot, kde_fit) else p_rayleigh(xplot, sigma_hat)
  
  ggplot(data.frame(x = xplot, F = Fvals), aes(x = x, y = F)) +
    geom_line(size = 0.9) +
    labs(
      x = expression(X[w]),
      y = expression(F(X[w]))
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey90")
    )
}

# --- 9. Monte Carlo: KDE vs Rayleigh MSE comparison ----------

set.seed(123)

MC     <- 500
grid_n <- 100

seasons      <- c("winter", "spring", "summer", "autumn")
mse_results  <- tibble()
results_list <- list()

output_folder <- "x_w_CDF"
dir.create(output_folder, showWarnings = FALSE)

for (yy in years) {
  for (ss in seasons) {
    
    cat("MC Processing:", yy, ss, "\n")
    
    sub   <- x_w_df %>% filter(year == yy, season == ss, x_w > 0)
    x_obs <- sub$x_w
    n     <- length(x_obs)
    
    if (n < 20) {
      cat("  --> skipped (not enough data)\n")
      next
    }
    
    # Fit models
    kde_fit   <- kde1d(x_obs)
    sigma_hat <- sqrt(sum(x_obs^2) / (2 * n))   # MLE for Rayleigh sigma
    
    # Evaluation grid
    x_grid <- seq(min(x_obs), max(x_obs), length.out = grid_n)
    F_obs  <- ecdf_on_grid(x_obs, x_grid)
    
    mse_kde_mc <- numeric(MC)
    mse_ray_mc <- numeric(MC)
    
    for (b in 1:MC) {
      x_kde_sim <- rkde1d(n, kde_fit)
      F_kde_sim <- ecdf_on_grid(x_kde_sim, x_grid)
      
      x_ray_sim <- r_rayleigh(n, sigma_hat)
      F_ray_sim <- ecdf_on_grid(x_ray_sim, x_grid)
      
      mse_kde_mc[b] <- mean((F_obs - F_kde_sim)^2)
      mse_ray_mc[b] <- mean((F_obs - F_ray_sim)^2)
    }
    
    mse_kde <- mean(mse_kde_mc)
    mse_ray <- mean(mse_ray_mc)
    winner  <- ifelse(mse_kde < mse_ray, "KDE", "Rayleigh")
    
    mse_results <- bind_rows(
      mse_results,
      tibble(
        year         = yy,
        season       = ss,
        mse_kde      = mse_kde,
        mse_rayleigh = mse_ray,
        winner       = winner
      )
    )
    
    # CDF function (MC-selected model)
    cdf_fun <- if (winner == "KDE") {
      function(xx) pkde1d(xx, kde_fit)
    } else {
      function(xx) p_rayleigh(xx, sigma_hat)
    }
    
    # PIT transform and CSV export
    sub <- sub %>% mutate(x_w_u = cdf_fun(x_w))
    write_csv(
      sub %>% dplyr::select(time, x_w, year, season, x_w_u),
      file.path(output_folder, paste0("x_w_cdf_", yy, "_", ss, ".csv"))
    )
    
    results_list[[paste0(yy, "_", ss)]] <- list(
      winner    = winner,
      kde       = kde_fit,
      sigma     = sigma_hat
    )
  }
}

# Summary of Monte Carlo results
mse_results %>% count(winner)

cat("\n=== DONE === CDF files saved to:", output_folder, "\n")

# --- 10. Example plots: Winter 2000 (Figures 6) --------

key       <- "2000_winter"
sub_ex    <- x_w_df %>% filter(year == 2000, season == "winter")   # full range for CDF
kde_fit   <- results_list[[key]]$kde
sigma_hat <- results_list[[key]]$sigma
winner    <- results_list[[key]]$winner

plot_wind_side_by_side_mc(sub_ex$x_w, kde_fit, sigma_hat, winner)
plot_wind_cdf_mc(sub_ex$x_w, kde_fit, sigma_hat, winner)

