# =============================================================
# Copula Fitting: Joe-Clayton (BB7), Tail Dependence,
# Bootstrap Uncertainty — Nighttime
# Paper: Quantifying Renewable Energy Resource Uncertainty
# Note: X_s is predominantly zero at night; the primary pair
#       of interest is therefore (X_r, X_w).
# =============================================================
rm(list = ls())

# --- Libraries -----------------------------------------------
library(VineCopula)
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)
library(future)
library(future.apply)
library(MASS)    # for kde2d
library(plot3D)  # for persp3D

# --- 1. Configuration ----------------------------------------

path_r <- "C:/Users/natti/OneDrive/Desktop/PhD Dissertation/chapter1_uncertainty/submission_jrss/code/x_r_CDF"
path_s <- "C:/Users/natti/OneDrive/Desktop/PhD Dissertation/chapter1_uncertainty/submission_jrss/code/x_s_CDF_nighttime"
path_w <- "C:/Users/natti/OneDrive/Desktop/PhD Dissertation/chapter1_uncertainty/submission_jrss/code/x_w_CDF"

parallel_run        <- TRUE
n_workers           <- max(1, parallel::detectCores(logical = FALSE) - 1)
min_sample_size     <- 20
n_boot              <- 500
min_sample_for_boot <- 10

# --- 2. Helper functions -------------------------------------

# Load all CSVs from a folder and tag with variable name
load_cdf_folder <- function(path, varname) {
  files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No CSV files found in: ", path)
  map_df(files, function(f) {
    read_csv(f, show_col_types = FALSE) %>% mutate(variable = varname)
  })
}

# Row-order merge of two season-year subsets
merge_two_marginals_rowwise <- function(df1, df2, u1_name, u2_name) {
  n <- min(nrow(df1), nrow(df2))
  if (n == 0) return(tibble(x1_u = numeric(0), x2_u = numeric(0)))
  tibble(x1_u = df1[[u1_name]][1:n], x2_u = df2[[u2_name]][1:n])
}

# Copula fit: family selected based on sign of Kendall's tau
# family 9  = BB7 (Joe-Clayton) for positive dependence
# family 24 = Rotated Gumbel 90 degrees for negative dependence
fit_joe_clayton <- function(u, v, refine_mle = FALSE) {
  ok <- which(!is.na(u) & !is.na(v))
  u  <- as.numeric(u[ok])
  v  <- as.numeric(v[ok])
  if (length(u) < min_sample_size) return(NULL)

  tau_hat <- cor(u, v, method = "kendall")
  fam     <- if (tau_hat >= 0) 9 else 24

  tryCatch(
    BiCopEst(u, v, family = fam, method = "mle"),
    error = function(e) { message("MLE failed: ", e$message); NULL }
  )
}

# Extract fit summary as a tibble row
extract_summary <- function(fit) {
  if (is.null(fit)) return(tibble(family = NA_integer_, tau = NA_real_, par = NA_real_, par2 = NA_real_, loglik = NA_real_))
  tibble(
    family = fit$family,
    tau    = BiCopPar2Tau(fit$family, fit$par, fit$par2),
    par    = fit$par,
    par2   = fit$par2,
    loglik = fit$loglik
  )
}

# Empirical tail dependence at quantile level p
empirical_tail_dep <- function(u1, u2, p = 0.1) {
  ok <- which(!is.na(u1) & !is.na(u2))
  u1 <- u1[ok]; u2 <- u2[ok]
  n  <- length(u1)
  if (n == 0) return(list(lambdaL = NA_real_, lambdaU = NA_real_))
  denomL  <- sum(u1 <= p)
  denomU  <- sum(u1 >= 1 - p)
  lambdaL <- if (denomL == 0) NA_real_ else sum(u1 <= p & u2 <= p) / denomL
  lambdaU <- if (denomU == 0) NA_real_ else sum(u1 >= 1 - p & u2 >= 1 - p) / denomU
  list(lambdaL = lambdaL, lambdaU = lambdaU)
}

# Compute tau + theoretical and empirical tail dependence for a fitted pair
compute_tails <- function(fit, u1, u2) {
  if (is.null(fit) || length(u1) < 2) return(rep(NA_real_, 5))

  tail_theo <- tryCatch(
    BiCopPar2TailDep(fit$family, fit$par, fit$par2),
    error = function(e) list(lower = NA_real_, upper = NA_real_)
  )

  tail_emp <- empirical_tail_dep(u1, u2, p = 0.1)
  c(
    BiCopPar2Tau(fit$family, fit$par, fit$par2),
    tail_theo$lower,
    tail_theo$upper,
    tail_emp$lambdaL,
    tail_emp$lambdaU
  )
}

# Get cleaned (u, v) pair for a given year-season and pair label
get_clean_pair <- function(year_val, season_val, pair) {
  r_sub <- cdf_r %>% filter(year == year_val, season == season_val)
  s_sub <- cdf_s %>% filter(year == year_val, season == season_val)
  w_sub <- cdf_w %>% filter(year == year_val, season == season_val)

  tmp <- switch(pair,
    rs = merge_two_marginals_rowwise(r_sub, s_sub, "x_r_u", "x_s_u"),
    rw = merge_two_marginals_rowwise(r_sub, w_sub, "x_r_u", "x_w_u"),
    sw = merge_two_marginals_rowwise(s_sub, w_sub, "x_s_u", "x_w_u"),
    stop("Unknown pair: ", pair)
  )

  ok <- which(
    !is.na(tmp$x1_u) & !is.na(tmp$x2_u) &
    is.finite(tmp$x1_u) & is.finite(tmp$x2_u) &
    tmp$x1_u > 0 & tmp$x1_u < 1 &
    tmp$x2_u > 0 & tmp$x2_u < 1
  )
  if (length(ok) == 0) return(list(u = numeric(0), v = numeric(0)))
  list(u = as.numeric(tmp$x1_u[ok]), v = as.numeric(tmp$x2_u[ok]))
}

# Bootstrap Kendall's tau with 95% CI
bootstrap_pair_stats <- function(u, v, n_boot = 500) {
  n <- length(u)
  if (n < min_sample_for_boot) return(list(mean = NA_real_, sd = NA_real_, lo = NA_real_, hi = NA_real_))

  tau_boot <- future.apply::future_lapply(seq_len(n_boot), function(b) {
    idx <- sample.int(n, size = n, replace = TRUE)
    uu  <- u[idx]; vv <- v[idx]
    if (length(unique(uu)) <= 1 || length(unique(vv)) <= 1) return(NA_real_)
    tryCatch(as.numeric(cor(uu, vv, method = "kendall", use = "complete.obs")), error = function(e) NA_real_)
  }, future.seed = TRUE)

  tau_boot <- unlist(tau_boot)
  tau_boot <- tau_boot[!is.na(tau_boot)]
  if (length(tau_boot) == 0) return(list(mean = NA_real_, sd = NA_real_, lo = NA_real_, hi = NA_real_))

  list(
    mean = mean(tau_boot),
    sd   = sd(tau_boot),
    lo   = as.numeric(quantile(tau_boot, 0.025, na.rm = TRUE)),
    hi   = as.numeric(quantile(tau_boot, 0.975, na.rm = TRUE))
  )
}

# Save 3D copula density plot to PNG

save_persp <- function(z, u, v, title, pair, yy, ss) {
  dir.create(file.path("copula_plots", pair), recursive = TRUE, showWarnings = FALSE)
  file <- file.path("copula_plots", pair, paste0(pair, "_", yy, "_", ss, ".png"))
  png(file, width = 1400, height = 1000)
  plot3D::persp3D(x = u, y = v, z = z,
                  theta = 40, phi = 25,
                  main  = title,
                  xlab  = "u1", ylab = "u2", zlab = "c(u1, u2)")
  dev.off()
}

# Nighttime copula density plot - (X_r, X_w) pair only
plot_copula_density_night <- function(yy, ss) {
  r_sub <- cdf_r %>% filter(year == yy, season == ss)
  w_sub <- cdf_w %>% filter(year == yy, season == ss)
  dat   <- merge_two_marginals_rowwise(r_sub, w_sub, "x_r_u", "x_w_u")
  if (nrow(dat) < 2) return(invisible(NULL))
  z <- kde2d(dat$x1_u, dat$x2_u, n = 40)
  save_persp(z$z, z$x, z$y,
             title = paste0("RW Night (", yy, " ", ss, ")"),
             pair  = "RW_night", yy = yy, ss = ss)
}

# --- 3. Data loading and alignment ---------------------------

cdf_r <- load_cdf_folder(path_r, "x_r")
cdf_s <- load_cdf_folder(path_s, "x_s")
cdf_w <- load_cdf_folder(path_w, "x_w")

# Add row IDs within year-season groups
cdf_r <- cdf_r %>% group_by(year, season) %>% mutate(row_id = row_number())
cdf_s <- cdf_s %>% group_by(year, season) %>% mutate(row_id = row_number())
cdf_w <- cdf_w %>% group_by(year, season) %>% mutate(row_id = row_number())

# Standardise season names to title case
cdf_r <- cdf_r %>% mutate(season = str_to_title(season))
cdf_s <- cdf_s %>% mutate(season = str_to_title(season))
cdf_w <- cdf_w %>% mutate(season = str_to_title(season))

# Align r and w to nighttime timestamps from s
nighttime_keys <- cdf_s %>% dplyr::select(time, year, season)
cdf_r <- cdf_r %>% semi_join(nighttime_keys, by = c("time", "year", "season"))
cdf_w <- cdf_w %>% semi_join(nighttime_keys, by = c("time", "year", "season"))

message("Rows after alignment - r: ", nrow(cdf_r), "  s: ", nrow(cdf_s), "  w: ", nrow(cdf_w))

# Build year-season combos present in all three datasets
years   <- sort(Reduce(intersect, list(unique(cdf_r$year), unique(cdf_s$year), unique(cdf_w$year))))
seasons <- sort(Reduce(intersect, list(unique(cdf_r$season), unique(cdf_s$season), unique(cdf_w$season))))
combos  <- expand.grid(year = years, season = seasons, stringsAsFactors = FALSE) %>% arrange(year, season)

message("Processing ", nrow(combos), " year-season combinations")

# Example: Spring 2010 nighttime (X_r, X_w) density
plot_copula_density_night(2010, "Spring")

# --- 4. Copula fitting loop ----------------------------------

if (parallel_run && n_workers > 1) plan(multisession, workers = n_workers) else plan(sequential)

res_list <- future.apply::future_lapply(seq_len(nrow(combos)), function(i) {
  yy <- combos$year[i]
  ss <- combos$season[i]

  r_sub <- cdf_r %>% filter(year == yy, season == ss)
  s_sub <- cdf_s %>% filter(year == yy, season == ss)
  w_sub <- cdf_w %>% filter(year == yy, season == ss)

  message("Fitting: ", yy, " ", ss, " | n: r=", nrow(r_sub), " s=", nrow(s_sub), " w=", nrow(w_sub))

  if (nrow(r_sub) == 0 || nrow(s_sub) == 0 || nrow(w_sub) == 0)
    return(list(key = paste0(yy, "_", ss), fits = NULL))

  rs <- merge_two_marginals_rowwise(r_sub, s_sub, "x_r_u", "x_s_u")
  rw <- merge_two_marginals_rowwise(r_sub, w_sub, "x_r_u", "x_w_u")
  sw <- merge_two_marginals_rowwise(s_sub, w_sub, "x_s_u", "x_w_u")

  fit_rs <- tryCatch(fit_joe_clayton(rs$x1_u, rs$x2_u), error = function(e) NULL)
  fit_rw <- tryCatch(fit_joe_clayton(rw$x1_u, rw$x2_u), error = function(e) NULL)
  fit_sw <- tryCatch(fit_joe_clayton(sw$x1_u, sw$x2_u), error = function(e) NULL)

  list(key = paste0(yy, "_", ss), fits = list(rs = fit_rs, rw = fit_rw, sw = fit_sw))
}, future.seed = TRUE)

plan(sequential)

# Convert to named list
results <- list()
for (el in res_list) results[[el$key]] <- el$fits

# --- 5. Build copula summary dataframe -----------------------

summary_df <- map_df(names(results), function(name) {
  tmp        <- results[[name]]
  year_val   <- as.numeric(strsplit(name, "_")[[1]][1])
  season_val <- strsplit(name, "_")[[1]][2]

  rs_s <- extract_summary(tmp$rs)
  rw_s <- extract_summary(tmp$rw)
  sw_s <- extract_summary(tmp$sw)

  tibble(
    year_season = name,
    year        = year_val,
    season      = season_val,
    rs_tau      = rs_s$tau,
    rs_par      = rs_s$par,
    rs_par2     = rs_s$par2,
    rs_loglik   = rs_s$loglik,
    rw_tau      = rw_s$tau,
    sw_tau      = sw_s$tau
  )
})

summary(summary_df)

# --- 6. Tail dependence for all pairs ------------------------

unc_list <- list()

for (name in names(results)) {
  tmp        <- results[[name]]
  year_val   <- as.numeric(strsplit(name, "_")[[1]][1])
  season_val <- strsplit(name, "_")[[1]][2]

  r_sub <- cdf_r %>% filter(year == year_val, season == season_val)
  s_sub <- cdf_s %>% filter(year == year_val, season == season_val)
  w_sub <- cdf_w %>% filter(year == year_val, season == season_val)

  row <- tibble(
    year_season = name, year = year_val, season = season_val,
    rs_tau = NA_real_, rs_lambdaL_theo = NA_real_, rs_lambdaU_theo = NA_real_,
    rs_lambdaL_emp = NA_real_, rs_lambdaU_emp = NA_real_,
    rw_tau = NA_real_, rw_lambdaL_theo = NA_real_, rw_lambdaU_theo = NA_real_,
    rw_lambdaL_emp = NA_real_, rw_lambdaU_emp = NA_real_,
    sw_tau = NA_real_, sw_lambdaL_theo = NA_real_, sw_lambdaU_theo = NA_real_,
    sw_lambdaL_emp = NA_real_, sw_lambdaU_emp = NA_real_
  )

  if (!is.null(tmp$rs)) {
    n_rs <- min(nrow(r_sub), nrow(s_sub))
    if (n_rs >= 2) {
      res <- compute_tails(tmp$rs, r_sub$x_r_u[1:n_rs], s_sub$x_s_u[1:n_rs])
      row$rs_tau <- res[1]; row$rs_lambdaL_theo <- res[2]; row$rs_lambdaU_theo <- res[3]
      row$rs_lambdaL_emp <- res[4]; row$rs_lambdaU_emp <- res[5]
    }
  }
  if (!is.null(tmp$rw)) {
    n_rw <- min(nrow(r_sub), nrow(w_sub))
    if (n_rw >= 2) {
      res <- compute_tails(tmp$rw, r_sub$x_r_u[1:n_rw], w_sub$x_w_u[1:n_rw])
      row$rw_tau <- res[1]; row$rw_lambdaL_theo <- res[2]; row$rw_lambdaU_theo <- res[3]
      row$rw_lambdaL_emp <- res[4]; row$rw_lambdaU_emp <- res[5]
    }
  }
  if (!is.null(tmp$sw)) {
    n_sw <- min(nrow(s_sub), nrow(w_sub))
    if (n_sw >= 2) {
      res <- compute_tails(tmp$sw, s_sub$x_s_u[1:n_sw], w_sub$x_w_u[1:n_sw])
      row$sw_tau <- res[1]; row$sw_lambdaL_theo <- res[2]; row$sw_lambdaU_theo <- res[3]
      row$sw_lambdaL_emp <- res[4]; row$sw_lambdaU_emp <- res[5]
    }
  }

  unc_list[[name]] <- row
}

uncertainty_df <- bind_rows(unc_list) %>%
  mutate(across(ends_with("_emp"), ~ ifelse(is.nan(.), NA, .)))

# Seasonal averages — empirical upper tails, theoretical lower tails enforced as 0
seasonal_tail <- uncertainty_df %>%
  group_by(season) %>%
  summarise(
    lambda_U_rs = mean(rs_lambdaU_emp, na.rm = TRUE),
    lambda_U_rw = mean(rw_lambdaU_emp, na.rm = TRUE),
    lambda_U_sw = mean(sw_lambdaU_emp, na.rm = TRUE),
    lambda_L_rs = 0,
    lambda_L_rw = 0,
    lambda_L_sw = mean(sw_lambdaL_emp, na.rm = TRUE)
  ) %>%
  ungroup()

print(seasonal_tail, width = Inf)

# --- 7. Save copula results to CSV ---------------------------

write.csv(summary_df,     "copula_summary_bb7_nighttime.csv",     row.names = FALSE)
write.csv(uncertainty_df, "copula_uncertainty_bb7_nighttime.csv", row.names = FALSE)

message("Saved summary (", nrow(summary_df), " rows) and uncertainty (", nrow(uncertainty_df), " rows)")

# --- 8. Bootstrap Kendall's tau for all 3 pairs --------------

if (parallel_run && n_workers > 1) plan(multisession, workers = n_workers) else plan(sequential)

bootstrap_results <- vector("list", length = length(names(results)))
names(bootstrap_results) <- names(results)

for (name in names(results)) {
  year_val   <- as.numeric(strsplit(name, "_")[[1]][1])
  season_val <- strsplit(name, "_")[[1]][2]

  rs_uv <- get_clean_pair(year_val, season_val, "rs")
  rw_uv <- get_clean_pair(year_val, season_val, "rw")
  sw_uv <- get_clean_pair(year_val, season_val, "sw")

  res_rs <- bootstrap_pair_stats(rs_uv$u, rs_uv$v, n_boot = n_boot)
  res_rw <- bootstrap_pair_stats(rw_uv$u, rw_uv$v, n_boot = n_boot)
  res_sw <- bootstrap_pair_stats(sw_uv$u, sw_uv$v, n_boot = n_boot)

  bootstrap_results[[name]] <- list(
    year_season = name,
    year        = year_val,
    season      = season_val,
    rs_mean = res_rs$mean, rs_sd = res_rs$sd, rs_lo = res_rs$lo, rs_hi = res_rs$hi,
    rw_mean = res_rw$mean, rw_sd = res_rw$sd, rw_lo = res_rw$lo, rw_hi = res_rw$hi,
    sw_mean = res_sw$mean, sw_sd = res_sw$sd, sw_lo = res_sw$lo, sw_hi = res_sw$hi
  )
}

plan(sequential)

boot_df <- bind_rows(bootstrap_results) %>% as_tibble()

boot_df %>%
  group_by(season) %>%
  summarise(rw = mean(rw_sd, na.rm = TRUE)) %>%
  mutate(season = factor(season, 
                         levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  arrange(season)

message("Bootstrap complete - rows: ", nrow(boot_df))
summary(boot_df$rs_sd)
summary(boot_df$rw_sd)
summary(boot_df$sw_sd)

# --- 9. Merge bootstrap CIs into plot_df ---------------------

plot_df <- summary_df %>%
  left_join(
    boot_df %>% dplyr::select(year_season, rs_lo, rs_hi, rs_sd, rw_lo, rw_hi, rw_sd, sw_lo, sw_hi, sw_sd),
    by = "year_season"
  ) %>%
  mutate(year = as.numeric(year))

# --- 10. 12-panel Kendall's tau plot (Figure 14 nighttime) ---
# Kendall's tau plot - nighttime (X_r, X_w) only 
# At nighttime X_s = 0 by construction; only the (X_r, X_w)
# pair is informative.

plot_long_night <- plot_df %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  dplyr::select(year, season, tau = rw_tau, lo = rw_lo, hi = rw_hi)

ggplot(plot_long_night %>% mutate(pair = "X[r] - X[w]"),
       aes(x = year, y = tau)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "steelblue") +
  geom_line(color  = "black", size = 0.6) +
  geom_point(color = "black", size = 1.5) +
  facet_grid(
    pair ~ season,
    labeller = labeller(pair = label_parsed)
  ) +
  labs(
    x     = "Year",
    y     = expression(tau),
    title = "Nighttime"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

# --- 11. Uncertainty visualisation ---------------------------

# Mean SD by season and pair (bar chart)
sd_summary <- boot_df %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  group_by(season) %>%
  summarise(
    rs = mean(rs_sd, na.rm = TRUE),
    rw = mean(rw_sd, na.rm = TRUE),
    sw = mean(sw_sd, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(rs, rw, sw), names_to = "pair", values_to = "mean_sd") %>%
  mutate(pair = factor(pair, levels = c("rs", "rw", "sw")))

pair_labels_bar <- c(
  rs = expression(X[r]-X[s]),
  rw = expression(X[r]-X[w]),
  sw = expression(X[s]-X[w])
)

ggplot(sd_summary, aes(x = season, y = mean_sd, fill = pair)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(
    values = c("rs" = "steelblue", "rw" = "darkorange", "sw" = "forestgreen"),
    labels = pair_labels_bar
  ) +
  labs(
    x     = "Season",
    y     = expression("Mean SD of " * tau),
    fill  = "Pair",
    title = "Meteorological-Driven Uncertainty by Season (Nighttime)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

# Top 5 highest uncertainty season-year combinations per pair
boot_df %>%
  dplyr::select(year_season, year, season, rs_sd, rw_sd, sw_sd) %>%
  pivot_longer(cols = c(rs_sd, rw_sd, sw_sd),
               names_to = "pair", values_to = "sd") %>%
  mutate(pair = gsub("_sd", "", pair)) %>%
  group_by(pair) %>%
  slice_max(sd, n = 5) %>%
  arrange(pair, desc(sd))

getwd()
