# =============================================================
# Sun Irradiance (X_s): KDE vs Lognormal, Monte Carlo, CDF
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
library(MASS)   # for fitdistr

# --- 1. Data loading -----------------------------------------

folder_path <- "D:/PhD/thesis/chapter_1/datasets/x_s/x_s_all"

years <- 2000:2022
file_names <- paste0("x_s_fra_", years, ".csv")
file_paths <- file.path(folder_path, file_names)

filtered_data <- list()
for (file in file_paths) {
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    subset <- data %>%
      filter(latitude == 45, longitude == 6.75) %>%   # Mont Blanc
      dplyr::select(time, x_s = ssrd) %>%
      mutate(time = ymd_hms(time))
    filtered_data[[file]] <- subset
  } else {
    warning(paste("File does not exist:", file))
  }
}

# Check number of valid files loaded
length(filtered_data)

# --- 2. Data preprocessing -----------------------------------

x_s_df <- bind_rows(filtered_data) %>%
  arrange(time) %>%
  mutate(
    x_s    = pmax(x_s / 100, 0),   # scale and remove negatives
    roll30 = rollmean(x_s, k = 24 * 30, fill = NA),  # 30-day rolling mean
    year   = as.integer(lubridate::year(time)),
    month  = lubridate::month(time),
    season = case_when(
      month %in% c(12, 1, 2)  ~ "winter",
      month %in% c(3, 4, 5)   ~ "spring",
      month %in% c(6, 7, 8)   ~ "summer",
      month %in% c(9, 10, 11) ~ "autumn"
    )
  )

# --- 3. Summary statistics (Table 2) -------------------------

summary_stats <- x_s_df %>%
  summarise(
    mean = mean(x_s, na.rm = TRUE),
    sd   = sd(x_s, na.rm = TRUE),
    p25  = quantile(x_s, 0.25, na.rm = TRUE),
    p50  = quantile(x_s, 0.50, na.rm = TRUE),
    p75  = quantile(x_s, 0.75, na.rm = TRUE)
  )

print(summary_stats)

# --- 4. Time series plot with rolling mean (Figure 1a) -------

x_s_df_sampled <- x_s_df %>%
  slice(seq(1, n(), by = 24 * 30))

p_raw <- ggplot(x_s_df_sampled, aes(x = time)) +
  geom_line(
    aes(y = x_s, color = "Raw data", linetype = "Raw data"),
    size = 0.7
  ) +
  geom_line(
    aes(y = roll30, color = "30-day rolling mean", linetype = "30-day rolling mean"),
    size = 1
  ) +
  labs(
    x = "Year",
    y = expression(X[s]),
    color = "Legend",
    linetype = "Legend"
  ) +
  scale_color_manual(values = c(
    "Raw data"          = "black",
    "30-day rolling mean" = "blue"
  )) +
  scale_linetype_manual(values = c(
    "Raw data"          = "solid",
    "30-day rolling mean" = "solid"
  )) +
  scale_x_datetime(labels = date_format("%Y"), date_breaks = "2 years") +
  theme_minimal()

p_raw

# --- 5. Helper functions -------------------------------------

# Daytime hours by season
daytime_hours <- list(
  winter = 8:17,
  spring = 6:18,
  summer = 6:20,
  autumn   = 7:17
)

# Extract season-year daytime subset using meteorological season boundaries
get_season_data <- function(df, year, season = c("winter", "spring", "summer", "autumn")) {
  season <- match.arg(season)

  if (season == "winter") {
    start <- ymd_hms(paste0(year - 1, "-12-21 00:00:00"))
    end   <- ymd_hms(paste0(year,     "-03-20 23:59:59"))
    hours <- 8:17
  } else if (season == "spring") {
    start <- ymd_hms(paste0(year, "-03-21 00:00:00"))
    end   <- ymd_hms(paste0(year, "-06-20 23:59:59"))
    hours <- 6:18
  } else if (season == "summer") {
    start <- ymd_hms(paste0(year, "-06-21 00:00:00"))
    end   <- ymd_hms(paste0(year, "-09-22 23:59:59"))
    hours <- 6:20
  } else if (season == "autumn") {
    start <- ymd_hms(paste0(year, "-09-23 00:00:00"))
    end   <- ymd_hms(paste0(year, "-12-20 23:59:59"))
    hours <- 7:17
  }

  df %>%
    filter(time >= start & time <= end) %>%
    mutate(hour = hour(time)) %>%
    filter(hour %in% hours)
}

# Empirical CDF evaluated on a grid
ecdf_on_grid <- function(x, grid) {
  ecdf(x)(grid)
}

# Lognormal sampler for Monte Carlo
r_lognormal <- function(n, meanlog, sdlog) {
  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}

# --- 6. Side-by-side PDF plot: KDE vs Lognormal (Figure 5) --

plot_season_side_by_side <- function(obj, title = NULL) {
  if (is.null(obj$kde) | is.null(obj$lognormal)) return(NULL)

  # KDE grid
  kde_df <- data.frame(
    x = obj$kde$grid_points,
    y = obj$kde$values
  )

  # Empirical density
  emp    <- density(obj$kde$x, from = 0)
  emp_df <- data.frame(x = emp$x, y = emp$y)

  # Lognormal density
  xgrid  <- seq(0, max(kde_df$x), length.out = 300)
  log_y  <- dlnorm(xgrid,
                   meanlog = obj$lognormal$estimate["meanlog"],
                   sdlog   = obj$lognormal$estimate["sdlog"])
  log_df <- data.frame(x = xgrid, y = log_y)

  # Plot 1: KDE vs Empirical
  p1 <- ggplot() +
    geom_line(data = kde_df, aes(x = x, y = y, color = "KDE"),
              size = 0.8) +
    geom_line(data = emp_df, aes(x = x, y = y, color = "Empirical"),
              linetype = "dashed", size = 0.5) +
    labs(x = expression(X[s]), y = "Density", color = "Legend") +
    scale_color_manual(values = c("KDE" = "darkgreen", "Empirical" = "black")) +
    theme_minimal() +
    theme(
      plot.title  = element_text(face = "plain", size = 14),
      axis.title  = element_text(face = "plain", size = 12),
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_line(color = "grey90", size = 0.25)
    ) +
    xlim(0, max(kde_df$x))

  # Plot 2: Lognormal vs Empirical
  p2 <- ggplot() +
    geom_line(data = log_df, aes(x = x, y = y, color = "Lognormal"),
              size = 0.8) +
    geom_line(data = emp_df, aes(x = x, y = y, color = "Empirical"),
              linetype = "dashed", size = 0.5) +
    labs(x = expression(X[s]), y = "Density", color = "Legend") +
    scale_color_manual(values = c("Lognormal" = "cornflowerblue", "Empirical" = "black")) +
    theme_minimal() +
    theme(
      plot.title  = element_text(face = "plain", size = 14),
      axis.title  = element_text(face = "plain", size = 12),
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_line(color = "grey90", size = 0.25)
    ) +
    xlim(0, max(kde_df$x))

  grid.arrange(p1, p2, ncol = 2)
}

# Example: Figure 5 (Winter 2000)
winter2000_fit <- list(
  kde      = kde1d(get_season_data(x_s_df, 2000, "winter")$x_s %>% .[. > 0]),
  lognormal = fitdistr(
    get_season_data(x_s_df, 2000, "winter")$x_s %>% .[. > 0],
    "lognormal"
  )
)
plot_season_side_by_side(winter2000_fit, title = "Winter 2000")

# --- 7. Monte Carlo: KDE vs Lognormal MSE comparison ---------

set.seed(123)

MC     <- 500
grid_n <- 100

mc_results  <- tibble()
results_list <- list()

for (yr in 2000:2022) {
  for (season in c("winter", "spring", "summer", "autumn")) {

    cat("MC Processing:", yr, season, "\n")

    season_df <- get_season_data(x_s_df, yr, season)
    x_obs     <- season_df$x_s
    x_obs     <- x_obs[x_obs > 0 & !is.na(x_obs)]

    n <- length(x_obs)
    if (n < 20) {
      cat("  --> skipped (not enough data)\n")
      next
    }

    # Fit models
    kde_fit <- kde1d(x_obs)
    log_fit <- fitdistr(x_obs, "lognormal")
    meanlog <- log_fit$estimate["meanlog"]
    sdlog   <- log_fit$estimate["sdlog"]

    # Evaluation grid
    x_grid <- seq(min(x_obs), max(x_obs), length.out = grid_n)
    F_obs  <- ecdf_on_grid(x_obs, x_grid)

    mse_kde_mc  <- numeric(MC)
    mse_logn_mc <- numeric(MC)

    for (b in 1:MC) {
      x_kde_sim    <- rkde1d(n, kde_fit)
      F_kde_sim    <- ecdf_on_grid(x_kde_sim, x_grid)

      x_log_sim    <- r_lognormal(n, meanlog, sdlog)
      F_log_sim    <- ecdf_on_grid(x_log_sim, x_grid)

      mse_kde_mc[b]  <- mean((F_obs - F_kde_sim)^2)
      mse_logn_mc[b] <- mean((F_obs - F_log_sim)^2)
    }

    mse_kde  <- mean(mse_kde_mc)
    mse_logn <- mean(mse_logn_mc)
    winner   <- ifelse(mse_kde < mse_logn, "KDE", "Lognormal")

    mc_results <- bind_rows(
      mc_results,
      tibble(
        year          = yr,
        season        = season,
        mse_kde       = mse_kde,
        mse_lognormal = mse_logn,
        winner        = winner
      )
    )

    results_list[[paste0(yr, "_", season)]] <- list(
      winner    = winner,
      kde       = kde_fit,
      lognormal = log_fit
    )
  }
}

# Summary of Monte Carlo results
mc_results %>% count(winner)

# --- 8. CDF via probability integral transform ---------------

# Trapezoid integration of KDE to produce CDF function
kde_to_cdf <- function(kde_obj) {
  x   <- kde_obj$grid_points
  pdf <- kde_obj$values

  dx  <- diff(x)
  cdf <- c(0, cumsum(0.5 * dx * (pdf[-1] + pdf[-length(pdf)])))
  cdf <- cdf / max(cdf)   # normalise to [0, 1]

  function(xx) approx(x, cdf, xout = xx, rule = 2)$y
}

# --- 9. Daytime CDF export (season-year) ---------------------

output_folder <- "x_s_CDF_daytime"
dir.create(output_folder, showWarnings = FALSE)

cdf_results <- list()

for (yy in years) {
  for (ss in c("winter", "spring", "summer", "autumn")) {

    cat("Processing daytime CDF:", yy, ss, "...\n")

    sub <- x_s_df %>%
      dplyr::filter(
        year   == yy,
        season == ss,
        lubridate::hour(time) %in% daytime_hours[[ss]]
      )

    if (nrow(sub) < 5 | sum(sub$x_s > 0) < 5) {
      cat("  --> skipped (not enough data)\n")
      next
    }

    kde_fit <- kde1d(sub$x_s)
    cdf_fun <- kde_to_cdf(kde_fit)

    sub <- sub %>%
      dplyr::mutate(x_s_u = cdf_fun(x_s)) %>%
      dplyr::select(time, x_s, season, year, x_s_u)

    outfile <- file.path(output_folder, paste0("x_s_cdf_day_", yy, "_", ss, ".csv"))
    readr::write_csv(sub, outfile)

    cdf_results[[paste0(yy, "_", ss)]] <- list(kde = kde_fit, cdf_fun = cdf_fun)
  }
}

cat("\n=== DONE === Daytime CDFs saved to:", output_folder, "\n")

# --- 10. Nighttime CDF export (season-year) ------------------

night_output_folder <- "x_s_CDF_nighttime"
dir.create(night_output_folder, showWarnings = FALSE)

cdf_results_night <- list()

for (yy in years) {
  for (ss in c("winter", "spring", "summer", "autumn")) {

    cat("Processing nighttime CDF:", yy, ss, "...\n")

    sub <- x_s_df %>%
      dplyr::filter(
        year   == yy,
        season == ss,
        !lubridate::hour(time) %in% daytime_hours[[ss]]
      )

    if (nrow(sub) < 5) {
      cat("  --> skipped (not enough data)\n")
      next
    }

    kde_fit <- kde1d(sub$x_s)
    cdf_fun <- kde_to_cdf(kde_fit)

    sub <- sub %>%
      dplyr::mutate(x_s_u = cdf_fun(x_s)) %>%
      dplyr::select(time, x_s, season, year, x_s_u)

    outfile <- file.path(night_output_folder, paste0("x_s_cdf_night_", yy, "_", ss, ".csv"))
    readr::write_csv(sub, outfile)

    cdf_results_night[[paste0(yy, "_", ss)]] <- list(kde = kde_fit, cdf_fun = cdf_fun)
  }
}

cat("\n=== DONE === Nighttime CDFs saved to:", night_output_folder, "\n")

# --- 11. CDF plot---------------------------------

plot_sun_cdf_simple <- function(data_vector, title = NULL) {
  data_vector <- data_vector[data_vector > 0 & !is.na(data_vector)]
  if (length(data_vector) < 2) return(NULL)

  kde_fit  <- kde1d(data_vector)
  dx       <- diff(kde_fit$grid_points)
  cdf_vals <- c(0, cumsum(0.5 * dx * (kde_fit$values[-1] + kde_fit$values[-length(kde_fit$values)])))
  cdf_vals <- cdf_vals / max(cdf_vals)
  xplot    <- kde_fit$grid_points

  qplot(xplot, cdf_vals, geom = "line") +
    labs(
      x     = expression(X[s]),
      y     = expression(F(X[s])),
      title = title
    ) +
    theme_bw() +
    theme(
      panel.grid.major  = element_line(color = "grey90"),
      panel.grid.minor  = element_line(color = "grey90"),
      panel.background  = element_rect(fill = "white"),
      plot.title        = element_text(hjust = 0.5)
    )
}

# Example: Winter 2000 daytime CDF
winter2000_day <- x_s_df %>%
  dplyr::filter(
    year   == 2000,
    season == "winter",
    lubridate::hour(time) %in% daytime_hours[["winter"]],
    x_s > 0
  ) %>%
  dplyr::pull(x_s)

plot_sun_cdf_simple(winter2000_day, title = "")
