# =============================================================
# Total Precipitation (X_r): KDE vs Gumbel, Monte Carlo, CDF
# Paper: Quantifying Renewable Energy Resource Uncertainty
# =============================================================
rm(list=ls())
# --- Libraries -----------------------------------------------
library(dplyr)
library(readr)
library(lubridate)
library(kde1d)
library(ggplot2)
library(gridExtra)

# --- 1. Data loading -----------------------------------------

folder_path <- "D:/PhD/thesis/chapter_1/datasets/x_r/x_r_all"

years      <- 2000:2022
file_names <- paste0("x_r_fra_", years, ".csv")
file_paths <- file.path(folder_path, file_names)

filtered_data <- list()
for (file in file_paths) {
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    subset <- data %>%
      filter(latitude == 45, longitude == 6.75) %>%   # Mont Blanc grid cell
      dplyr::select(time, x_r = tp) %>%
      mutate(
        time = ymd_hms(time),
        x_r  = pmax(x_r, 0),                           # remove negatives
        x_r  = x_r * 10000                             # scale to paper units
      )
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

x_r_df <- bind_rows(filtered_data) %>%
  arrange(time) %>%
  mutate(
    year   = year(time),
    month  = month(time),
    season = case_when(
      month %in% c(12, 1, 2)  ~ "winter",
      month %in% c(3, 4, 5)   ~ "spring",
      month %in% c(6, 7, 8)   ~ "summer",
      month %in% c(9, 10, 11) ~ "autumn"
    )
  ) %>%
  filter(year %in% years)

# Confirm no negative values remain
min(x_r_df$x_r)

# --- 3. Helper functions -------------------------------------

# Gumbel MOM estimator (method of moments)
fit_gumbel <- function(x) {
  b <- sd(x) * sqrt(6) / pi
  a <- mean(x) - b * 0.5772156649   # Euler-Mascheroni constant
  list(a = a, b = b)
}

# Gumbel PDF
d_gumbel <- function(x, a, b) {
  (1 / b) * exp(-(x - a) / b) * exp(-exp(-(x - a) / b))
}

# Gumbel CDF
p_gumbel <- function(x, a, b) {
  exp(-exp(-(x - a) / b))
}

# Gumbel sampler
r_gumbel <- function(n, a, b) {
  a - b * log(-log(runif(n)))
}

# Empirical CDF evaluated on a grid
ecdf_on_grid <- function(x, grid) {
  ecdf(x)(grid)
}

# Trapezoidal integration of KDE grid to produce CDF values
kde_to_cdf_grid <- function(xgrid, ygrid) {
  dx       <- diff(xgrid)
  areas    <- (ygrid[-1] + ygrid[-length(ygrid)]) / 2 * dx
  cdf_vals <- c(0, cumsum(areas))
  cdf_vals <- cdf_vals / max(cdf_vals, na.rm = TRUE)   # normalise to [0, 1]
  data.frame(x = xgrid, cdf = cdf_vals)
}

# Convert kde1d or density() object to an interpolating CDF function
make_kde_cdf <- function(kde_obj) {
  if (!is.null(kde_obj$grid_points)) {
    xgrid <- kde_obj$grid_points
    ygrid <- kde_obj$values
  } else {
    xgrid <- kde_obj$x
    ygrid <- kde_obj$y
  }

  o     <- order(xgrid)
  xgrid <- xgrid[o]
  ygrid <- ygrid[o]

  cdf_df <- kde_to_cdf_grid(xgrid, ygrid)

  function(q) {
    q   <- as.numeric(q)
    res <- approx(x = cdf_df$x, y = cdf_df$cdf, xout = q, rule = 2)$y
    res[q < min(cdf_df$x)] <- 0
    res[q > max(cdf_df$x)] <- 1
    res
  }
}

# Zero-inflated mixture CDF: point mass at zero + continuous KDE part
make_mixture_cdf_zero_inflated <- function(kde_obj, p_zero) {
  cdf_cont <- make_kde_cdf(kde_obj)
  function(q) {
    q    <- as.numeric(q)
    atom <- ifelse(q >= 0, p_zero, 0)
    cont <- (1 - p_zero) * cdf_cont(q)
    pmin(pmax(atom + cont, 0), 1)
  }
}

# --- 4. Monte Carlo: KDE vs Gumbel MSE comparison ------------

set.seed(123)

seasons      <- c("winter", "spring", "summer", "autumn")
MC           <- 500
grid_n       <- 100

mse_results  <- tibble()
results_list <- list()

output_folder <- "x_r_CDF"
dir.create(output_folder, showWarnings = FALSE)

for (yy in years) {
  for (ss in seasons) {

    cat("MC Processing:", yy, ss, "\n")

    sub_all <- x_r_df %>% filter(year == yy, season == ss)
    x_obs   <- sub_all$x_r[sub_all$x_r > 0]
    n       <- length(x_obs)

    if (n < 20) {
      cat("  --> skipped (not enough data)\n")
      next
    }

    # Fit models
    kde_fit    <- kde1d(x_obs)
    gumbel_fit <- fit_gumbel(x_obs)
    a <- gumbel_fit$a
    b <- gumbel_fit$b

    # Evaluation grid
    x_grid <- seq(min(x_obs), max(x_obs), length.out = grid_n)
    F_obs  <- ecdf_on_grid(x_obs, x_grid)

    mse_kde_mc    <- numeric(MC)
    mse_gumbel_mc <- numeric(MC)

    for (b_mc in 1:MC) {
      x_kde_sim <- rkde1d(n, kde_fit)
      F_kde_sim <- ecdf_on_grid(x_kde_sim, x_grid)

      x_gum_sim <- r_gumbel(n, a, b)
      F_gum_sim <- ecdf_on_grid(x_gum_sim, x_grid)

      mse_kde_mc[b_mc]    <- mean((F_obs - F_kde_sim)^2)
      mse_gumbel_mc[b_mc] <- mean((F_obs - F_gum_sim)^2)
    }

    mse_kde    <- mean(mse_kde_mc)
    mse_gumbel <- mean(mse_gumbel_mc)
    winner     <- ifelse(mse_kde < mse_gumbel, "KDE", "Gumbel")

    mse_results <- bind_rows(
      mse_results,
      tibble(
        year       = yy,
        season     = ss,
        mse_kde    = mse_kde,
        mse_gumbel = mse_gumbel,
        winner     = winner
      )
    )

    results_list[[paste0(yy, "_", ss)]] <- list(
      winner = winner,
      kde    = kde_fit,
      gumbel = gumbel_fit
    )
  }
}

# Summary of Monte Carlo results
mse_results %>% count(winner)

# --- 5. Side-by-side PDF plot: KDE vs Gumbel (Figure 7) ------

plot_rain_side_by_side_mc <- function(data_vector, kde_fit, gumbel_fit, winner) {
  data_vector <- data_vector[!is.na(data_vector) & data_vector > 0]
  if (length(data_vector) < 5) return(NULL)

  emp_density <- density(data_vector, from = 0)
  emp_df      <- data.frame(x = emp_density$x, y = emp_density$y)

  kde_df <- data.frame(x = kde_fit$grid_points, y = kde_fit$values)

  xgrid     <- seq(0, max(kde_df$x), length.out = 300)
  gumbel_df <- data.frame(
    x = xgrid,
    y = d_gumbel(xgrid, gumbel_fit$a, gumbel_fit$b)
  )

  # Plot 1: KDE vs Empirical
  p1 <- ggplot() +
    geom_line(data = kde_df,
              aes(x = x, y = y, color = "KDE"), size = 0.8) +
    geom_line(data = emp_df,
              aes(x = x, y = y, color = "Empirical"),
              linetype = "dashed", size = 0.5) +
    labs(x = expression(X[r]), y = "Density", color = "Legend") +
    scale_color_manual(values = c("KDE" = "darkgreen", "Empirical" = "black")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_line(color = "grey90", size = 0.25)
    ) +
    xlim(0, max(kde_df$x))

  # Plot 2: Gumbel vs Empirical
  p2 <- ggplot() +
    geom_line(data = gumbel_df,
              aes(x = x, y = y, color = "Gumbel"), size = 0.8) +
    geom_line(data = emp_df,
              aes(x = x, y = y, color = "Empirical"),
              linetype = "dashed", size = 0.5) +
    labs(x = expression(X[r]), y = "Density", color = "Legend") +
    scale_color_manual(values = c("Gumbel" = "cornflowerblue", "Empirical" = "black")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_line(color = "grey90", size = 0.25)
    ) +
    xlim(0, max(kde_df$x))

  grid.arrange(p1, p2, ncol = 2)
}

# Example: Figure 7 (Winter 2000)
key        <- "2000_winter"
sub_ex     <- x_r_df %>% filter(year == 2000, season == "winter", x_r > 0)
kde_fit_ex <- results_list[[key]]$kde
gum_fit_ex <- results_list[[key]]$gumbel
winner_ex  <- results_list[[key]]$winner

plot_rain_side_by_side_mc(sub_ex$x_r, kde_fit_ex, gum_fit_ex, winner_ex)

# --- 6. Mixture CDF export (season-year) ---------------------

cdf_results <- list()

cdf_output_folder <- "x_r_CDF"
dir.create(cdf_output_folder, showWarnings = FALSE)

for (yy in years) {
  for (ss in seasons) {

    cat("Processing mixture CDF:", yy, ss, "...\n")

    sub <- x_r_df %>% filter(year == yy, season == ss)

    if (nrow(sub) < 5) {
      cat("  --> skipped (not enough data)\n")
      next
    }

    p_zero        <- mean(sub$x_r == 0)
    positive_data <- sub$x_r[sub$x_r > 0]

    if (length(positive_data) < 5) {
      cat("  --> skipped (not enough positive values)\n")
      next
    }

    # KDE on positive values only
    kde_pos <- kde1d(positive_data)

    # Zero-inflated mixture CDF
    cdf_fun <- make_mixture_cdf_zero_inflated(kde_obj = kde_pos, p_zero = p_zero)

    # PIT transform: randomise over [0, p_zero] for zero observations
    sub <- sub %>%
      mutate(
        x_r_u = ifelse(
          x_r == 0,
          runif(n(), 0, p_zero),
          p_zero + (1 - p_zero) * make_kde_cdf(kde_pos)(x_r)
        )
      )

    outfile <- file.path(cdf_output_folder, paste0("x_r_cdf_", yy, "_", ss, ".csv"))
    write_csv(sub %>% dplyr::select(time, x_r, year, season, x_r_u), outfile)

    cdf_results[[paste0(yy, "_", ss)]] <- list(
      cdf_fun = cdf_fun,
      p_zero  = p_zero,
      kde     = kde_pos
    )
  }
}

cat("\n=== DONE === Mixture CDFs saved to:", cdf_output_folder, "\n")

# --- 7. CDF plot example -------------------------

season_data   <- x_r_df %>% filter(year == 2000, season == "winter") %>% pull(x_r)
p_zero_ex     <- mean(season_data == 0)
positive_ex   <- season_data[season_data > 0]
kde_pos_ex    <- kde1d(positive_ex)
cdf_fun_ex    <- make_mixture_cdf_zero_inflated(kde_obj = kde_pos_ex, p_zero = p_zero_ex)

xplot    <- seq(0, max(positive_ex, na.rm = TRUE), length.out = 200)
cdf_vals <- cdf_fun_ex(xplot)

qplot(xplot, cdf_vals, geom = "line") +
  labs(
    x = expression(X[r]),
    y = expression(F(X[r]))
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.background = element_rect(fill = "white")
  )
