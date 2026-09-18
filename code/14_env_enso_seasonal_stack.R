# seasonal ENSO composites of environmental drivers

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(ggnewscale)
  library(data.table)
})

# paths
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

tz_local <- "America/Lima"
dry_months <- 5:10
cutoff_datetime <- as.POSIXct("2024-10-15 23:59:59", tz = tz_local)

out_csv <- file.path(output_path, "env_ENSO_monthly_composites.csv")
out_png <- file.path(fig_path, "Fig2_env_ENSO_seasonal_stack.png")

flux_fp <- file.path(input_folder, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
if (!file.exists(flux_fp)) {
  stop("missing prepared flux file: ", flux_fp)
}

normalize_enso <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x0 %in% c("el nino", "el niño", "el niÃ±o") ~ "El Niño",
    x0 %in% c("la nina", "la niña", "la niÃ±a") ~ "La Niña",
    x0 %in% c("neutral") ~ "Neutral",
    TRUE ~ NA_character_
  )
}

normalize_season <- function(x, month_num = NULL) {
  x0 <- tolower(trimws(as.character(x)))
  out <- dplyr::case_when(
    startsWith(x0, "d") ~ "Dry",
    startsWith(x0, "w") ~ "Wet",
    TRUE ~ NA_character_
  )
  if (!is.null(month_num)) {
    out <- ifelse(is.na(out), ifelse(month_num %in% dry_months, "Dry", "Wet"), out)
  }
  out
}

parse_local_datetime <- function(x) {
  if (inherits(x, "POSIXct")) {
    return(lubridate::with_tz(x, tz_local))
  }
  out <- suppressWarnings(lubridate::ymd_hms(x, quiet = TRUE, tz = tz_local))
  if (all(is.na(out))) {
    out <- suppressWarnings(lubridate::parse_date_time(
      x,
      orders = c("Y-m-d H:M:S", "Y-m-d H:M", "Ymd HMS", "Ymd HM", "YmdHMS", "YmdHM"),
      tz = tz_local
    ))
  }
  out
}

first_existing <- function(nm, candidates) {
  hit <- intersect(candidates, nm)
  if (length(hit)) hit[1] else NA_character_
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) NA_real_ else stats::sd(x)
}

month_levels <- month.abb
enso_levels <- c("El Niño", "La Niña", "Neutral")
# ENSO colours now come from palette.R (enso_cols: gold / purple / grey)

read_flux_data <- function(path) {
  dat <- data.table::fread(path, data.table = FALSE)
  if (!"PAR_corrected_SWin" %in% names(dat) && "PPFD_IN_1_1_1" %in% names(dat)) {
    dat$PAR_corrected_SWin <- suppressWarnings(as.numeric(dat$PPFD_IN_1_1_1))
  }
  if (!"PAR" %in% names(dat) && "PAR_corrected_SWin" %in% names(dat)) {
    dat$PAR <- suppressWarnings(as.numeric(dat$PAR_corrected_SWin))
  }
  if (!"PAR" %in% names(dat)) {
    stop("missing MATLAB-corrected PAR column in: ", path)
  }
  par_mismatch <- is.finite(dat$PAR) & is.finite(dat$PAR_corrected_SWin) &
    abs(dat$PAR - dat$PAR_corrected_SWin) > 1e-8
  if (any(par_mismatch, na.rm = TRUE)) {
    stop("canonical PAR does not match PAR_corrected_SWin in: ", path)
  }
  dt_col <- first_existing(names(dat), c("DateTime", "tv_dt"))
  if (is.na(dt_col)) stop("could not find DateTime or tv_dt in: ", path)

  dat %>%
    mutate(
      DateTime = parse_local_datetime(.data[[dt_col]]),
      Date = as.Date(DateTime),
      Year = lubridate::year(DateTime),
      Month = lubridate::month(DateTime),
      Month_label = factor(month.abb[Month], levels = month_levels),
      ENSO = normalize_enso(ENSO),
      Season = normalize_season(season, Month)
    ) %>%
    filter(
      !is.na(DateTime),
      DateTime <= cutoff_datetime,
      ENSO %in% enso_levels,
      !is.na(Year),
      !is.na(Month)
    )
}

monthly_flux_variable <- function(dat, variable, source_col, unit, panel_label, value_type = "mean") {
  if (is.na(source_col) || !source_col %in% names(dat)) {
    message("skip ", variable, ": no matching column found")
    return(NULL)
  }

  out <- dat %>%
    transmute(
      Year,
      Month,
      Month_label,
      ENSO,
      value = suppressWarnings(as.numeric(.data[[source_col]]))
    ) %>%
    filter(is.finite(value)) %>%
    group_by(Year, Month, Month_label, ENSO) %>%
    summarise(value = mean_or_na(value), n_obs = dplyr::n(), .groups = "drop") %>%
    mutate(
      variable = variable,
      source_column = source_col,
      unit = unit,
      panel_label = panel_label,
      value_type = value_type
    )

  if (!nrow(out)) {
    message("skip ", variable, ": no finite values in ", source_col)
    return(NULL)
  }

  out
}

read_monthly_precip <- function(flux_dates) {
  imerg_csv <- list.files(input_folder, pattern = "(?i)imerg.*30min.*\\.csv$", full.names = TRUE)
  if (!length(imerg_csv)) {
    message("skip precipitation: no IMERG 30-minute CSV found in data/")
    return(NULL)
  }

  imerg_csv <- imerg_csv[order(file.info(imerg_csv)$mtime, decreasing = TRUE)][1]
  imerg <- data.table::fread(imerg_csv, data.table = FALSE)

  utc_col <- first_existing(names(imerg), c("time_utc_iso", "time_utc", "time", "system:time_start"))
  local_col <- first_existing(names(imerg), c("time_local_minus05", "time_local"))

  if (!is.na(local_col)) {
    DateTime <- suppressWarnings(lubridate::parse_date_time(
      imerg[[local_col]],
      orders = c("Y-m-d H:M:S", "Y-m-d H:M", "Ymd HMS", "Ymd HM", "YmdHMS", "YmdHM"),
      tz = tz_local
    ))
  } else if (!is.na(utc_col)) {
    t_utc <- suppressWarnings(lubridate::parse_date_time(
      imerg[[utc_col]],
      orders = c("Y-m-d H:M:S", "Y-m-d H:M", "Ymd HMS", "Ymd HM", "YmdHMS", "YmdHM"),
      tz = "UTC"
    ))
    DateTime <- lubridate::with_tz(t_utc, tz_local)
  } else {
    message("skip precipitation: no recognized IMERG time column")
    return(NULL)
  }

  if ("precip_mm_30min" %in% names(imerg)) {
    pr30 <- suppressWarnings(as.numeric(imerg$precip_mm_30min))
    precip_source <- "precip_mm_30min"
  } else {
    rate_col <- first_existing(names(imerg), c("precip_mm_per_hr", "precipitation", "precipitationCal"))
    if (is.na(rate_col)) {
      message("skip precipitation: no recognized precipitation column")
      return(NULL)
    }
    pr30 <- suppressWarnings(as.numeric(imerg[[rate_col]])) * 0.5
    precip_source <- paste0(rate_col, " converted from mm hr-1")
  }

  enso_by_date <- flux_dates %>%
    filter(!is.na(Date), !is.na(ENSO), !is.na(Season)) %>%
    distinct(Date, ENSO, Season)

  daily <- tibble(DateTime = DateTime, precip_mm_30min = pr30) %>%
    filter(!is.na(DateTime), DateTime <= cutoff_datetime, is.finite(precip_mm_30min)) %>%
    mutate(Date = as.Date(DateTime)) %>%
    group_by(Date) %>%
    summarise(precip_day_mm = sum(precip_mm_30min, na.rm = TRUE), n_obs = dplyr::n(), .groups = "drop") %>%
    left_join(enso_by_date, by = "Date") %>%
    mutate(
      Year = lubridate::year(Date),
      Month = lubridate::month(Date),
      Month_label = factor(month.abb[Month], levels = month_levels),
      ENSO = normalize_enso(ENSO),
      Season = normalize_season(Season, Month)
    ) %>%
    filter(ENSO %in% enso_levels)

  out <- daily %>%
    group_by(Year, Month, Month_label, ENSO) %>%
    summarise(value = sum(precip_day_mm, na.rm = TRUE), n_obs = sum(n_obs, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      variable = "Precipitation",
      source_column = precip_source,
      unit = "mm month^-1",
      panel_label = "atop(Precipitation, (mm~month^{-1}))",
      value_type = "monthly total"
    )

  if (!nrow(out)) {
    message("skip precipitation: no labeled finite monthly totals")
    return(NULL)
  }

  out
}

make_composite_table <- function(monthly_values) {
  monthly_values %>%
    group_by(variable, source_column, unit, panel_label, value_type, Month, Month_label, ENSO) %>%
    summarise(
      mean = mean_or_na(value),
      sd = sd_or_na(value),
      n_years = dplyr::n_distinct(Year[is.finite(value)]),
      n_monthly_values = sum(is.finite(value)),
      min_year = min(Year[is.finite(value)], na.rm = TRUE),
      max_year = max(Year[is.finite(value)], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = ifelse(n_years >= 2, mean - sd, NA_real_),
      ymax = ifelse(n_years >= 2, mean + sd, NA_real_),
      ENSO = factor(ENSO, levels = enso_levels),
      Month_label = factor(as.character(Month_label), levels = month_levels)
    )
}

# `base_size` defaults to the standalone figure. It is a parameter so the same
# stack can be re-rendered small enough to sit in the two-column Figure 2
# composite at 1:1 scale; rescaling the 7.2-inch raster into a 3.25-inch slot
# would leave the type at about 5 pt. Line and point sizes scale with it.
make_stack_plot <- function(composite, base_size = 11, legend_box = "horizontal") {
  sf <- base_size / 11
  panel_order <- c(
    "atop(italic(T)[plain(a)], (degree*C))",
    "atop(VPD, (kPa))",
    "atop(Precipitation, (mm~month^{-1}))",
    "atop(PAR, (mu*mol~m^{-2}~s^{-1}))",
    "atop(SWC, (cm))",
    "atop(italic(T)[plain(s)], (degree*C))"
  )
  composite <- composite %>%
    mutate(panel_label = factor(panel_label, levels = panel_order[panel_order %in% unique(panel_label)]))

  season_bands <- tibble(
    xmin = c(0.5, 4.5, 10.5),
    xmax = c(4.5, 10.5, 12.5),
    Season = factor(c("Wet", "Dry", "Wet"), levels = c("Wet", "Dry"))
  )

  ggplot(composite, aes(x = Month, y = mean, color = ENSO, linetype = ENSO, group = ENSO)) +
    # wet/dry season background bands (own fill scale -> "Season" legend)
    geom_rect(
      data = season_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Season),
      inherit.aes = FALSE,
      # 0.25, not the 0.6 used until now: season_box_fills are boxplot tints
      # chosen so BLACK lines stay readable inside a box, and at 0.6 they drove
      # the contrast of the coloured ENSO lines drawn on top down to 1.3-2.2,
      # against a minimum of 3.0. At 0.25 the bands are near-white washes and the
      # worst-case line contrast is 3.07. See the note in code/palette.R.
      alpha = 0.25
    ) +
    scale_fill_manual(
      values = season_box_fills,
      name = "Season",
      guide = guide_legend(order = 2, override.aes = list(alpha = 0.9))
    ) +
    ggnewscale::new_scale_fill() +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = ENSO), color = NA, alpha = 0.16, na.rm = TRUE) +
    geom_line(linewidth = 1.2 * sf, na.rm = TRUE) +
    geom_point(size = 2.0 * sf, stroke = 0.2, na.rm = TRUE) +
    facet_grid(panel_label ~ ., scales = "free_y", switch = "y", labeller = label_parsed) +
    scale_x_continuous(breaks = 1:12, labels = month_levels, expand = expansion(mult = c(0.01, 0.01))) +
    scale_color_manual(values = enso_cols, drop = FALSE, guide = guide_legend(order = 1)) +
    scale_fill_manual(values = enso_cols, drop = FALSE, guide = "none") +
    scale_linetype_manual(values = enso_linetypes, drop = FALSE, guide = guide_legend(order = 1)) +
    labs(
      x = NULL,
      y = NULL,
      color = "ENSO",
      linetype = "ENSO"
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = "top",
      # "vertical" stacks the ENSO and Season keys on separate rows; at the
      # narrow panel width the single-row arrangement runs off the canvas
      legend.box = legend_box,
      legend.title = element_text(size = base_size - 2, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(1.2 * sf, "lines"),
      # both reproduce the ggplot defaults at base_size 11 and shrink with sf
      legend.margin = margin(rep(5.5 * sf, 4)),
      legend.box.spacing = unit(11 * sf, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = base_size - 1),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(size = base_size - 2),
      panel.spacing.y = unit(0.45, "lines")
    )
}

flux <- read_flux_data(flux_fp)

variable_specs <- tibble::tribble(
  ~variable, ~candidates, ~unit, ~panel_label,
  "VPD", c("VPD_kPa", "VPD"), "kPa", "atop(VPD, (kPa))",
  "SWC", c("SWC_1_1_1", "SWC"), "cm", "atop(SWC, (cm))",
  "PAR", c("PAR", "PAR_corrected_SWin"), "umol m^-2 s^-1", "atop(PAR, (mu*mol~m^{-2}~s^{-1}))",
  "Ta", c("TA_1_1_1", "Ta", "Tair"), "°C", "atop(italic(T)[plain(a)], (degree*C))",
  "Ts", c("TS_3", "Ts"), "°C", "atop(italic(T)[plain(s)], (degree*C))"
)

flux_monthly <- bind_rows(lapply(seq_len(nrow(variable_specs)), function(i) {
  spec <- variable_specs[i, ]
  source_col <- first_existing(names(flux), spec$candidates[[1]])
  monthly_flux_variable(
    dat = flux,
    variable = spec$variable,
    source_col = source_col,
    unit = spec$unit,
    panel_label = spec$panel_label
  )
}))

precip_monthly <- read_monthly_precip(flux %>% select(Date, ENSO, Season))
monthly_values <- bind_rows(precip_monthly, flux_monthly) %>%
  mutate(
    ENSO = factor(ENSO, levels = enso_levels),
    Month_label = factor(as.character(Month_label), levels = month_levels)
  )

if (!nrow(monthly_values)) {
  stop("no environmental variables were available for plotting")
}

composite <- make_composite_table(monthly_values) %>%
  arrange(variable, ENSO, Month)

write.csv(composite, out_csv, row.names = FALSE)

p_stack <- make_stack_plot(composite)
ggsave(out_png, p_stack, width = 7.2, height = 9.4, dpi = 300)

# small-format copy for the two-column Figure 2 composite (code/15), authored at
# its final size so it can be placed 1:1
out_png_panel <- file.path(fig_path, "Fig2c_env_stack_panel.png")
ggsave(out_png_panel, make_stack_plot(composite, base_size = 7, legend_box = "vertical"),
       width = 3.25, height = 5.92, dpi = 300)
message("saved Figure 2 panel copy: ", out_png_panel)

included <- composite %>%
  distinct(variable, source_column, unit) %>%
  arrange(variable)

message("wrote composite table: ", out_csv)
message("saved figures:")
message("  - ", out_png)
message("included variables:")
for (i in seq_len(nrow(included))) {
  message("  - ", included$variable[i], " from ", included$source_column[i], " (", included$unit[i], ")")
}
