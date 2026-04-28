# =============================================================
# Solar Irradiance (X_s): Summary Statistics and Figures
# Paper: Quantifying Renewable Energy Resource Uncertainty
# =============================================================
rm(list=ls())
# --- Libraries -----------------------------------------------
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(scales)
library(zoo)

# --- 1. Data loading -----------------------------------------

folder_path <- "D:/PhD/thesis/chapter_1/datasets/x_s/x_s_all"
years      <- 2000:2022
file_names <- paste0("x_s_fra_", years, ".csv")
file_paths <- file.path(folder_path, file_names)

filtered_data <- list()
for (file in file_paths) {
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    subset <- data %>%
      filter(latitude == 45, longitude == 6.75) %>%   # Mont Blanc grid cell
      dplyr::select(time, x_s = ssrd) %>%
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

x_s_df <- bind_rows(filtered_data) %>%
  arrange(time) %>%
  mutate(
    x_s    = pmax(x_s / 100, 0),                       # scale and remove negatives
    roll30 = rollmean(x_s, k = 24 * 30, fill = NA)     # 30-day rolling mean
  )

# --- 3. Summary statistics (Table 2) -------------------------

summary_stats <- x_s_df %>%
  summarise(
    mean = mean(x_s, na.rm = TRUE),
    sd   = sd(x_s,   na.rm = TRUE),
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
    x        = "Year",
    y        = expression(X[s]),
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

# --- 5. Zero-inflation PDF plot (Figure 2) -------------------
# Full summer 2000 data (daytime + nighttime) to show mass at zero

summer_all <- x_s_df %>%
  filter(
    time >= ymd_hms("2000-06-21 00:00:00"),
    time <= ymd_hms("2000-09-22 23:59:59")
  )

ggplot(summer_all, aes(x = x_s)) +
  geom_density(fill = "grey", alpha = 0.3) +
  labs(
    x = expression(X[s]),
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(hjust = 0.5),
    panel.grid  = element_blank()
  )
