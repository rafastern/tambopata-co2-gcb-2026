# ───────────────────────────────────────────────────────────────────────────────
# resampling + within-condition ANOVA for ENSO × season subsets
# now reads: dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv
# (season, ENSO, NEE_ok, Reco, GEP are assumed to already exist in the input)
#
# key behavior:
#   - analysis filter uses NEE_ok (not raw NEE)
#   - drops empty ENSO×season groups automatically
#   - resampling:
#       * summary resampling uses a per-group fraction (default 0.8)
#       * single-iteration export uses Fix A: n_samples_safe = min group size
# ───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
  library(ggplot2)
})

# ───────────────────────────────────────────────────────────────────────────────
# paths
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
graphs_path <- file.path(graphs_path, "random_resample_anova")

# 1) load data (enriched)
in_fp <- file.path(input_folder, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
stopifnot(file.exists(in_fp))

df <- readr::read_csv(in_fp, show_col_types = FALSE)

# ───────────────────────────────────────────────────────────────────────────────
# parse DateTime + diagnostics (try BOTH interpretations)
if (!("DateTime" %in% names(df))) stop("missing DateTime column in: ", in_fp)

dt_str <- df$DateTime

# a) interpret strings as local Lima time (no conversion)
dt_local_assuming_local <- lubridate::ymd_hms(dt_str, tz = "America/Lima", quiet = TRUE)
if (all(is.na(dt_local_assuming_local))) {
  dt_local_assuming_local <- as.POSIXct(dt_str, tz = "America/Lima")
}

# b) interpret strings as UTC, then convert to Lima time
dt_utc <- lubridate::ymd_hms(dt_str, tz = "UTC", quiet = TRUE)
if (all(is.na(dt_utc))) {
  dt_utc <- as.POSIXct(dt_str, tz = "UTC")
}
dt_local_from_utc <- lubridate::with_tz(dt_utc, "America/Lima")


# diel PAR helper
stopifnot("PAR" %in% names(df))
if (!"PAR_corrected_SWin" %in% names(df) && "PPFD_IN_1_1_1" %in% names(df)) {
  df$PAR_corrected_SWin <- suppressWarnings(as.numeric(df$PPFD_IN_1_1_1))
}
if (!"PAR_corrected_SWin" %in% names(df)) {
  stop("missing PAR_corrected_SWin in enriched input")
}
par_mismatch <- is.finite(df$PAR) & is.finite(df$PAR_corrected_SWin) & abs(df$PAR - df$PAR_corrected_SWin) > 1e-8
if (any(par_mismatch, na.rm = TRUE)) {
  stop("canonical PAR does not match PAR_corrected_SWin in enriched input")
}

make_diel_par <- function(dt_vec, par_vec, label) {
  tibble(DateTime = dt_vec, PAR = as.numeric(par_vec)) %>%
    filter(!is.na(DateTime), is.finite(PAR)) %>%
    mutate(
      DateTime_30 = floor_date(DateTime, unit = "30 minutes"),
      hour = format(DateTime_30, "%H:%M:%S")
    ) %>%
    group_by(hour) %>%
    summarise(mean_par = mean(PAR, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    mutate(mode = label)
}

diel_local <- make_diel_par(dt_local_assuming_local, df$PAR, "DateTime parsed as America/Lima (no conversion)")
diel_utc   <- make_diel_par(dt_local_from_utc,       df$PAR, "DateTime parsed as UTC -> converted to America/Lima")

diel_all <- bind_rows(diel_local, diel_utc)

halfhour_levels <- sprintf("%02d:%02d:00", rep(0:23, each = 2), rep(c(0, 30), times = 24))
diel_all$hour <- factor(diel_all$hour, levels = halfhour_levels)

p_par <- ggplot(diel_all, aes(x = hour, y = mean_par, group = mode, color = mode)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.4) +
  scale_x_discrete(breaks = c("02:00:00", "08:00:00", "14:00:00", "20:00:00")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "top") +
  labs(title = "corrected PAR diel cycle check after reading enriched CSV",
       y = expression(PAR ~ (mu*mol ~ m^-2 ~ s^-1)),
       color = NULL)

print(p_par)


cat("\ncorrected PAR peak hour check:\n")
print(diel_all %>% group_by(mode) %>% slice_max(mean_par, n = 1, with_ties = FALSE))
print(diel_all %>% group_by(mode) %>% slice_max(mean_par, n = 3, with_ties = FALSE) %>% arrange(mode, desc(mean_par)))

# choose the correct DateTime for the rest of the script:
# if the UTC->local curve looks right, use dt_local_from_utc; otherwise use dt_local_assuming_local.
df$DateTime <- dt_local_from_utc
# ───────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────


# required grouping columns
req_cols <- c("season", "ENSO")
missing_req <- setdiff(req_cols, names(df))
if (length(missing_req)) stop("missing required columns: ", paste(missing_req, collapse = ", "))

# ensure numeric columns are numeric (safe coercion)
num_cols <- intersect(c("NEE","NEE_f","NEE_ok","Reco","GEP"), names(df))
for (cc in num_cols) df[[cc]] <- as.numeric(df[[cc]])

# ───────────────────────────────────────────────────────────────────────────────
# quick diagnostics
cat("\n==== diagnostics ====\n")
cat("input:", in_fp, "\n")
cat("rows:", nrow(df), "\n")
cat("DateTime NA:", sum(is.na(df$DateTime)), "\n")
cat("columns present:", paste(names(df), collapse = ", "), "\n")
for (cc in num_cols) {
  cat(cc, "finite:", sum(is.finite(df[[cc]])), " / ", nrow(df), "\n")
}
cat("=====================\n\n")

# derive year/month just in case (overwrites if already present)
df$year  <- year(df$DateTime)
df$month <- month(df$DateTime)

# normalize labels to prevent mismatches
df <- df %>%
  mutate(
    season = tolower(trimws(season)),
    ENSO   = trimws(ENSO)
  )

# keep chronological order
df <- df %>% arrange(DateTime)

# ───────────────────────────────────────────────────────────────────────────────
# analysis filter: use NEE_ok
df_clean <- df[is.finite(df$NEE_ok), ]

# 2) create seasonal–ENSO subsets
df_wet_el_nino <- df_clean %>% filter(season == "wet", ENSO == "El Nino")
df_dry_el_nino <- df_clean %>% filter(season == "dry", ENSO == "El Nino")
df_wet_la_nina <- df_clean %>% filter(season == "wet", ENSO == "La Nina")
df_dry_la_nina <- df_clean %>% filter(season == "dry", ENSO == "La Nina")
df_wet_neutral <- df_clean %>% filter(season == "wet", ENSO == "neutral")
df_dry_neutral <- df_clean %>% filter(season == "dry", ENSO == "neutral")

# 3) define variables and parameters
iterations <- 100

variables_to_resample <- c("NEE", "NEE_f", "NEE_ok", "Reco", "GEP")
variables_to_resample <- variables_to_resample[variables_to_resample %in% names(df_clean)]

# 4) resampling function (mean distribution)
resample_and_test_multiple_vars <- function(data, condition_cols, value_cols, n_samples, iterations) {
  results <- list()
  
  for (value_col in value_cols) {
    variable_results <- rep(NA_real_, iterations)
    
    for (i in 1:iterations) {
      # if the group itself is too small, this iteration is NA
      if (nrow(data) < n_samples) {
        variable_results[i] <- NA_real_
        next
      }
      
      sampled_data <- data
      
      if (!is.null(condition_cols)) {
        sampled_data <- sampled_data %>%
          group_by(across(all_of(condition_cols))) %>%
          filter(n() >= n_samples) %>%
          sample_n(n_samples, replace = FALSE) %>%
          ungroup()
      } else {
        sampled_data <- sampled_data %>%
          sample_n(n_samples, replace = FALSE)
      }
      
      # extra safety
      if (nrow(sampled_data) < n_samples) {
        variable_results[i] <- NA_real_
        next
      }
      
      vals <- sampled_data[[value_col]]
      if (all(is.na(vals))) {
        variable_results[i] <- NA_real_
      } else {
        variable_results[i] <- mean(vals, na.rm = TRUE)
      }
    }
    
    results[[value_col]] <- variable_results
  }
  
  results
}

# 5) store all subsets in a list
df_list <- list(
  wet_el_nino = df_wet_el_nino,
  dry_el_nino = df_dry_el_nino,
  wet_la_nina = df_wet_la_nina,
  dry_la_nina = df_dry_la_nina,
  wet_neutral = df_wet_neutral,
  dry_neutral = df_dry_neutral
)

# 6) group sizes and drop empty groups
counts <- sapply(df_list, nrow)
cat("\nGroup sizes (rows with finite NEE_ok):\n")
print(counts)

df_list <- df_list[counts > 0]
counts  <- counts[counts > 0]

if (length(df_list) == 0) stop("all ENSO×season groups are empty after filtering by NEE_ok")

# choose resampling size strategy for the mean-distribution plots
sample_frac_group <- 0.8

# 7) apply resampling (mean distribution) to each subset
resampling_results <- list()

for (group_name in names(df_list)) {
  cat("\nProcessing:", group_name, "\n")
  df_group <- df_list[[group_name]]
  
  n_samples_group <- max(1L, floor(sample_frac_group * nrow(df_group)))
  cat("n_samples_group:", n_samples_group, " / group_n:", nrow(df_group), "\n")
  
  resampling_results[[group_name]] <- resample_and_test_multiple_vars(
    data = df_group,
    condition_cols = NULL,
    value_cols = variables_to_resample,
    n_samples = n_samples_group,
    iterations = iterations
  )
}

# check each group-variable combination
for (group_name in names(df_list)) {
  df_group <- df_list[[group_name]]
  cat("\nGroup:", group_name, "\n")
  for (variable in variables_to_resample) {
    non_na_count <- sum(is.finite(df_group[[variable]]))
    cat(variable, ": ", non_na_count, " finite values\n")
  }
}

# 8) summarize and plot resampling distributions
for (group_name in names(resampling_results)) {
  cat("\n===== Group:", group_name, "=====\n")
  
  for (variable in variables_to_resample) {
    cat("toggle variable:", variable, "\n")
    res <- resampling_results[[group_name]][[variable]]
    
    ok <- is.finite(res)
    if (!any(ok)) {
      cat("skipping: no finite resampled means for", variable, "in", group_name, "\n")
      next
    }
    
    print(summary(res[ok]))
    cat("Variance:", var(res[ok]), "\n")
    
    hist(
      res[ok],
      main = paste("Distribution of Means -", group_name, "-", variable),
      xlab = paste("Mean", variable)
    )
  }
}

###################################################################################################
# ANOVA
perform_within_condition_anova <- function(df_list, variables, group_size, n_groups) {
  set.seed(42)
  results <- list()
  
  for (condition_name in names(df_list)) {
    cat("\n===== Condition:", condition_name, "=====\n")
    df0 <- df_list[[condition_name]]
    
    for (variable in variables) {
      cat("\n--- Variable:", variable, "---\n")
      
      if (!(variable %in% names(df0))) {
        cat("missing variable", variable, "in", condition_name, "\n")
        next
      }
      
      dfv <- df0 %>% filter(is.finite(.data[[variable]]))
      
      if (nrow(dfv) < group_size) {
        cat("Not enough finite data in", condition_name, "for variable", variable, "\n")
        next
      }
      
      group_samples <- lapply(1:n_groups, function(i) {
        sampled <- dfv %>% sample_n(group_size, replace = FALSE)
        data.frame(
          value = sampled[[variable]],
          group = as.factor(paste0("G", i))
        )
      })
      
      anova_df <- do.call(rbind, group_samples)
      
      aov_model <- aov(value ~ group, data = anova_df)
      print(summary(aov_model))
      
      results[[paste(condition_name, variable, sep = "_")]] <- summary(aov_model)
    }
  }
  
  results
}

variables_to_test <- c("NEE", "NEE_f", "NEE_ok", "Reco", "GEP")
variables_to_test <- variables_to_test[variables_to_test %in% names(df_clean)]

group_size <- 800
n_groups <- 100

df_list_for_anova <- df_list[names(df_list) != "wet_el_nino"]
anova_within_results <- perform_within_condition_anova(df_list_for_anova, variables_to_test, group_size, n_groups)

######################################################################################################
# saving to load in other scripts

cols_to_keep <- c(
  "DateTime", "season", "ENSO", "year",
  "NEE", "NEE_f", "NEE_ok", "Reco", "GEP",
  "H", "LE",
  "TA_1_1_1", "WS_1_1_1", "WD_1_1_1", "USTAR", "VPD_kPa", "SWC_1_1_1", "TS_3",
  "PAR", "PAR_corrected_SWin", "PAR_source", "PPFD_IN_1_1_1",
  "SW_IN_1_1_1","NETRAD_1_1_1","FC"
)

# show what’s missing per group
print(lapply(df_list, function(x) setdiff(cols_to_keep, names(x))))

# re-select safely
df_dry_el_nino <- df_dry_el_nino %>% select(any_of(cols_to_keep))
df_wet_el_nino <- df_wet_el_nino %>% select(any_of(cols_to_keep))
df_wet_la_nina <- df_wet_la_nina %>% select(any_of(cols_to_keep))
df_dry_la_nina <- df_dry_la_nina %>% select(any_of(cols_to_keep))
df_wet_neutral <- df_wet_neutral %>% select(any_of(cols_to_keep))
df_dry_neutral <- df_dry_neutral %>% select(any_of(cols_to_keep))

save(
  df_dry_el_nino,
  df_dry_la_nina,
  df_dry_neutral,
  df_wet_el_nino,
  df_wet_la_nina,
  df_wet_neutral,
  file = file.path(output_path, "enso_conditions_data.RData")
)

######################################################################################################
# resample and save a single iteration (Fix A: n_samples_safe = min group size)

resample_and_save_single_iteration <- function(data, condition_cols, value_cols, n_samples, iterations, output_file) {
  final_resampled_data <- data.frame()
  
  for (value_col in value_cols) {
    if (!(value_col %in% names(data))) next
    if (all(is.na(data[[value_col]]))) next
    
    selected_iteration <- sample(1:iterations, 1)
    
    sampled_data <- data %>%
      group_by(across(all_of(condition_cols))) %>%
      filter(n() >= n_samples) %>%
      sample_n(n_samples, replace = FALSE) %>%
      ungroup()
    
    sampled_data$iteration <- selected_iteration
    sampled_data$value <- sampled_data[[value_col]]
    sampled_data$variable <- value_col
    
    final_resampled_data <- rbind(final_resampled_data, sampled_data)
  }
  
  tmp_output_file <- tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(output_file)), "_"),
    tmpdir = dirname(output_file),
    fileext = ".csv"
  )
  write.csv(final_resampled_data, file = tmp_output_file, row.names = FALSE)
  if (file.exists(output_file)) unlink(output_file)
  if (!file.rename(tmp_output_file, output_file)) {
    file.copy(tmp_output_file, output_file, overwrite = TRUE)
    unlink(tmp_output_file)
  }
  final_resampled_data
}

# compute safe n_samples for sampling within (season, ENSO)
counts2 <- df_clean %>% count(season, ENSO, name = "n")
cat("\nGroup sizes in df_clean (for Fix A):\n")
print(counts2)

n_samples_safe <- min(counts2$n)
cat("n_samples_safe (min group size):", n_samples_safe, "\n")

output_file <- file.path(output_path, "resampled_data.csv")

final_resampled_data <- resample_and_save_single_iteration(
  df_clean,
  c("season", "ENSO"),
  variables_to_resample,
  n_samples_safe,
  iterations,
  output_file
)

print(head(final_resampled_data))
cat("\nRows per group used in resampling distributions:\n")
print(sapply(df_list, nrow))
