# =============================================================
# Total Precipitation (X_r): Summary Statistics and Figures
# Paper: Quantifying Renewable Energy Resource Uncertainty
# Note: raw data files use the column name 'tp'; renamed to x_r throughout
# =============================================================
rm(list=ls())
# --- Libraries -----------------------------------------------
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)

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

x_r_df <- bind_rows(filtered_data) %>%
  arrange(time) %>%
  mutate(
    x_r    = pmax(x_r, 0),             # remove negatives
    x_r    = x_r * 10000,              # scale to paper units
    roll30 = rollmean(x_r, k = 24 * 30, fill = NA)   # 30-day rolling mean
  )

# --- 3. Summary statistics (Table 2) -------------------------

summary_stats <- x_r_df %>%
  summarise(
    mean = mean(x_r, na.rm = TRUE),
    sd   = sd(x_r,   na.rm = TRUE),
    p25  = quantile(x_r, 0.25, na.rm = TRUE),
    p50  = quantile(x_r, 0.50, na.rm = TRUE),
    p75  = quantile(x_r, 0.75, na.rm = TRUE)
  )

print(summary_stats)

# --- 4. Time series plot with rolling mean (Figure 1c) -------

x_r_df_sampled <- x_r_df %>%
  slice(seq(1, n(), by = 24 * 30))

p_raw <- ggplot(x_r_df_sampled, aes(x = time)) +
  geom_line(
    aes(y = x_r, color = "Raw data", linetype = "Raw data"),
    size = 0.7
  ) +
  geom_line(
    aes(y = roll30, color = "30-day rolling mean", linetype = "30-day rolling mean"),
    size = 1
  ) +
  labs(
    x        = "Year",
    y        = expression(X[r]),
    color    = "Legend",
    linetype = "Legend"
  ) +
  scale_color_manual(values = c(
    "Raw data"            = "black",
    "30-day rolling mean" = "blue"
  )) +
  scale_linetype_manual(values = c(
    "Raw data"            = "solid",
    "30-day rolling mean" = "solid"
  )) +
  scale_x_datetime(labels = date_format("%Y"), date_breaks = "2 years") +
  theme_minimal()

p_raw

# --- 5. Discrete-continuous mixture plot (Figure 4) ----------
# Summer 2000: zero-bar (left axis = Occurrence) +
#              continuous density (right axis = Density)

summer_precip <- x_r_df %>%
  filter(
    time >= ymd_hms("2000-06-21 00:00:00"),
    time <= ymd_hms("2000-09-20 23:59:59")
  ) %>%
  mutate(x_r = as.numeric(x_r))

zeros_count <- sum(summer_precip$x_r == 0, na.rm = TRUE)
pos         <- summer_precip %>% filter(x_r > 0)
n_pos       <- nrow(pos)

# Density for positive values
if (n_pos > 1) {
  d       <- density(pos$x_r, from = 0)
  dens_df <- data.frame(x = d$x, y = d$y)
} else {
  dens_df <- data.frame(x = numeric(0), y = numeric(0))
}

# Scale zero bar to same y-axis as density
Dmax     <- max(dens_df$y, na.rm = TRUE)
zero_bar <- data.frame(
  x        = 0,
  y_scaled = zeros_count / max(zeros_count, 1) * Dmax
)

ggplot() +
  # Zero-rainfall bar (left axis: Occurrence)
  geom_col(
    data  = zero_bar,
    aes(x = x, y = y_scaled),
    width = 0.5,
    fill  = "gray80", color = "black"
  ) +
  # Continuous density for positive values (right axis: Density)
  {if (nrow(dens_df) > 0)
    geom_area(
      data  = dens_df,
      aes(x = x, y = y),
      fill  = "gray70", alpha = 0.5
    )
  } +
  scale_y_continuous(
    name   = "Occurrence",
    breaks = seq(0, Dmax, length.out = 5),
    labels = function(v) round(v / Dmax * zeros_count),
    sec.axis = sec_axis(~ ., name = "Density")
  ) +
  labs(x = expression(X[r])) +
  theme_minimal() +
  theme(
    plot.title  = element_text(hjust = 0.5),
    panel.grid  = element_blank()
  )
