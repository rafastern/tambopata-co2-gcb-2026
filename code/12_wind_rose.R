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

dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

flux_fp <- file.path(input_path, "dataset_from_matlab.csv")

# ───────────────────────────────────────────────────────────────────────────────
# read data
df <- read_csv(flux_fp, show_col_types = FALSE, progress = FALSE)

ws_col <- "WS_1_1_1"
wd_col <- "WD_1_1_1"

# ───────────────────────────────────────────────────────────────────────────────
# prepare data
df2 <- df |>
  mutate(
    tv_dt = dmy_hms(tv_dt, tz = "America/Lima", quiet = TRUE),
    date = as.POSIXct(tv_dt, tz = "America/Lima"),
    hour_local = hour(date)
  ) |>
  filter(
    !is.na(date),
    !is.na(.data[[ws_col]]),
    !is.na(.data[[wd_col]]),
    .data[[ws_col]] >= 0,
    .data[[wd_col]] >= 0,
    .data[[wd_col]] <= 360
  )

# define day/night in local time
df_day <- df2 |>
  filter(hour_local >= 6, hour_local < 18)

df_night <- df2 |>
  filter(hour_local < 6 | hour_local >= 18)

# ───────────────────────────────────────────────────────────────────────────────
wind_speed_breaks <- c(0, 0.5, 1, 2, 3, 4, Inf)
wind_speed_labels <- c("0-0.5", "0.5-1", "1-2", "2-3", "3-4", ">4")
wind_dir_labels <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                     "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")

prepare_windrose_bins <- function(dat, facet_var = NULL) {
  d <- dat |>
    mutate(
      wd_sector = floor(((.data[[wd_col]] + 11.25) %% 360) / 22.5),
      wd_sector = factor(wd_sector, levels = 0:15, labels = wind_dir_labels),
      ws_bin = cut(
        .data[[ws_col]],
        breaks = wind_speed_breaks,
        labels = wind_speed_labels,
        right = FALSE,
        include.lowest = TRUE
      )
    ) |>
    filter(!is.na(wd_sector), !is.na(ws_bin))

  if (!is.null(facet_var)) {
    d |>
      count(.data[[facet_var]], wd_sector, ws_bin, name = "n") |>
      group_by(.data[[facet_var]]) |>
      mutate(percent = 100 * n / sum(n)) |>
      ungroup()
  } else {
    d |>
      count(wd_sector, ws_bin, name = "n") |>
      mutate(percent = 100 * n / sum(n))
  }
}

# base_size defaults to 11, the value every wind-rose figure was built with, so
# adding the argument leaves them byte-identical. The two explicit sizes below
# are expressed as ratios of base_size for the same reason: at 11 they evaluate
# to the original 7 and 9. Only the wet/dry pair (supplement Fig. S5) overrides
# it -- that figure is authored 12 in wide and placed at 180 mm, so its type was
# reaching the page at 4-5 pt.
make_windrose_plot <- function(dat, main_title, facet_var = NULL, base_size = 11) {
  rose <- prepare_windrose_bins(dat, facet_var)

  p <- ggplot(rose, aes(x = wd_sector, y = percent, fill = ws_bin)) +
    geom_col(width = 1, color = "grey35", linewidth = 0.15) +
    coord_polar(start = -pi / 16) +
    scale_fill_brewer(palette = "YlGnBu", direction = 1, name = "WS (m/s)") +
    labs(title = main_title, x = NULL, y = "Frequency (%)") +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = base_size * 7 / 11),
      axis.title.y = element_text(size = base_size * 9 / 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(stats::as.formula(paste("~", facet_var)))
  }

  p
}

# helper to save a wind rose
save_windrose <- function(dat, out_file, main_title) {
  p <- make_windrose_plot(
    dat |> mutate(month_name = factor(month.abb[month(date)], levels = month.abb)),
    main_title,
    facet_var = "month_name"
  )
  ggsave(out_file, p, width = 12, height = 9, dpi = 200, bg = "white")
  message("saved: ", out_file)
}

# ───────────────────────────────────────────────────────────────────────────────
# all data
save_windrose(
  df2,
  file.path(fig_path, "tambopata_windrose_by_month_panels_all.png"),
  "Tambopata - Wind roses by month - all data"
)

# daytime only
save_windrose(
  df_day,
  file.path(fig_path, "tambopata_windrose_by_month_panels_daytime.png"),
  "Tambopata - Wind roses by month - daytime"
)

# nighttime only
save_windrose(
  df_night,
  file.path(fig_path, "tambopata_windrose_by_month_panels_nighttime.png"),
  "Tambopata - Wind roses by month - nighttime"
)


# ───────────────────────────────────────────────────────────────────────────────
# additional wind-rose figures: all data combined, and wet/dry seasons

df2 <- df2 |>
  mutate(
    season = case_when(
      month(date) %in% 5:10 ~ "(b) Dry season",
      TRUE                  ~ "(a) Wet season"
    ),
    season = factor(season, levels = c("(a) Wet season", "(b) Dry season"))
  )

# helper to save a single-panel wind rose
save_windrose_single <- function(dat, out_file, main_title) {
  p <- make_windrose_plot(dat, main_title)
  ggsave(out_file, p, width = 9, height = 8, dpi = 200, bg = "white")
  message("saved: ", out_file)
}

# helper to save wet/dry wind rose panels
save_windrose_season <- function(dat, out_file) {
  # 18 matches supplement Fig. S3, which is also authored 12 in wide, so the two
  # figures carry the same effective type once placed at 180 mm.
  p <- make_windrose_plot(dat, "", facet_var = "season", base_size = 18)
  ggsave(out_file, p, width = 12, height = 7, dpi = 200, bg = "white")
  message("saved: ", out_file)
}

# single-panel wind rose: all years and months together
save_windrose_single(
  df2,
  file.path(fig_path, "tambopata_windrose_all_years_all_months_single_panel.png"),
  "Tambopata - Wind rose - all years and months"
)

# two-panel wind rose: wet versus dry season
save_windrose_season(
  df2,
  file.path(fig_path, "tambopata_windrose_wet_dry_seasons.png")
)
