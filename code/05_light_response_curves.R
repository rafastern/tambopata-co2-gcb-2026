# ───────────────────────────────────────────────────────────────────────────────
# light response curves (NEE_ok vs PAR) by ENSO × season
# - uses tambopata_48points_per_month.csv (monthly × half-hour bins)
# - fits rectangular hyperbola to FC_pos = -NEE_ok
# - colors points by chosen environmental variable (e.g., VPD, Ta, Hour)
# - applies cutoff using ym (since this file has no full DateTime)
# - outputs pooled (all-data) 2×3 grids (5 panels + legend panel)
# - no minimum-point filtering (beyond basic "enough variation" to fit)
# - no per-year logic
#
# requested layout changes:
#   - remove figure title
#   - use a single shared Y axis label (left side of whole figure)
#   - unify X axis label as one shared label at the bottom
#   - show X tick labels only on bottom row panels
# ───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(minpack.lm)
  library(gridExtra)
  library(grid)
  library(viridis)
  library(rlang)
  library(tidyr)
  library(data.table)
  library(readr)
})

# ───────────────────────────────────────────────────────────────────────────────
# paths
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
graphs_path         <- file.path(graphs_path, "NEE_light_response_no_precip")

dir.create(paper_figures_path, showWarnings = FALSE, recursive = TRUE)
dir.create(graphs_path, showWarnings = FALSE, recursive = TRUE)

# ───────────────────────────────────────────────────────────────────────────────
# config
cutoff_ym <- "2024-10"  # this file has ym (YYYY-MM), not DateTime

pretty_names <- c(
  wet_el_nino   = "El Niño – wet",
  dry_el_nino   = "El Niño – dry",
  wet_la_nina   = "La Niña – wet",
  dry_la_nina   = "La Niña – dry",
  wet_neutral   = "Neutral – wet",
  dry_neutral   = "Neutral – dry"
)

# display set (5 plotted panels; legend occupies the 6th slot)
panel_order_display  <- c("dry_el_nino", "dry_la_nina", "dry_neutral", "wet_el_nino", "wet_la_nina")
panel_labels_display <- c("(a)", "(b)", "(c)", "(d)", "(e)")
legend_panel_label   <- NULL

# figure font sizes
fig_base_size <- 16
fig_title_size <- 17
fig_subtitle_size <- 13
fig_axis_text_size <- 14
fig_shared_axis_size <- 18
fig_legend_title_size <- 18
fig_legend_text_size <- 15

# ───────────────────────────────────────────────────────────────────────────────
# load csv dataset
# nee_source selects the response: "measured" = non-gap-filled raw NEE (primary; matches the
#   non-gap-filled-primary contract and recovers the 15 months where the RP gap-fill is empty),
#   "gapfilled" = NEE_ok (= NEE_f), kept as a sensitivity option. gap-filled-vs-measured
#   comparison was run separately as a sensitivity check.
nee_source <- "measured"
if (identical(nee_source, "measured")) {
  csv_fp <- file.path(input_folder, "tambopata_48points_per_month_measured.csv")
  df_all <- readr::read_csv(csv_fp, show_col_types = FALSE)
  # internal response column NEE_ok holds the non-gap-filled measured NEE in this mode
  df_all$NEE_ok <- df_all$NEE_meas
} else {
  csv_fp <- file.path(input_folder, "tambopata_48points_per_month.csv")
  df_all <- readr::read_csv(csv_fp, show_col_types = FALSE)
}

# rename mapping (what exists in this csv)
rename_map <- c(
  "TA_1_1_1"      = "Ta",
  "VPD_kPa"       = "VPD",
  "SWC_1_1_1"     = "SWC",
  "TS_3"          = "Ts"
)

# ───────────────────────────────────────────────────────────────────────────────
# helpers

rename_with_map <- function(df, rename_map) {
  old <- names(rename_map)
  present <- old[old %in% names(df)]
  if (length(present)) {
    for (old_name in present) {
      new_name <- unname(rename_map[[old_name]])
      if (!new_name %in% names(df)) {
        names(df)[match(old_name, names(df))] <- new_name
      }
    }
  }
  df
}

# normalize text with transliteration (niño->nino, niña->nina)
norm_text <- function(x) {
  xx <- trimws(as.character(x))
  xx <- iconv(xx, from = "", to = "ASCII//TRANSLIT")
  tolower(xx)
}

norm_season <- function(x) {
  xx <- norm_text(x)
  dplyr::case_when(
    xx %in% c("wet","rainy","rain","w") ~ "wet",
    xx %in% c("dry","drought","d")     ~ "dry",
    grepl("wet", xx)                   ~ "wet",
    grepl("dry", xx)                   ~ "dry",
    TRUE                               ~ NA_character_
  )
}

norm_enso <- function(x) {
  xx <- norm_text(x)
  dplyr::case_when(
    grepl("la", xx) & grepl("nina", xx) ~ "la_nina",
    grepl("nina", xx)                  ~ "la_nina",
    grepl("el", xx) & grepl("nino", xx) ~ "el_nino",
    grepl("nino", xx)                  ~ "el_nino",
    grepl("neutral", xx)               ~ "neutral",
    xx %in% c("neu","n")               ~ "neutral",
    TRUE                               ~ NA_character_
  )
}

# build numeric hour from "HH:MM:SS" (or accept numeric hour if already numeric)
make_hour_numeric <- function(df) {
  if ("Hour" %in% names(df)) return(df)
  
  if ("hour" %in% names(df)) {
    h <- df$hour
    if (is.numeric(h)) {
      df$Hour <- as.numeric(h)
      return(df)
    }
    hh <- suppressWarnings(as.integer(substr(h, 1, 2)))
    mm <- suppressWarnings(as.integer(substr(h, 4, 5)))
    ok <- is.finite(hh) & is.finite(mm)
    df$Hour <- NA_real_
    df$Hour[ok] <- hh[ok] + mm[ok] / 60
  }
  df
}

apply_cutoff_ym <- function(df, cutoff_ym) {
  if (!"ym" %in% names(df)) return(df)
  df %>% filter(as.character(ym) <= cutoff_ym)
}

blank_panel <- function(title_text) {
  ggplot() +
    theme_void() +
    labs(title = title_text) +
    theme(
      plot.title = element_text(hjust = 0, size = 12, face = "bold"),
      plot.subtitle = element_blank()
    )
}

light_response <- function(Q, phi0, Pmax, Rd) {
  ((phi0 * Pmax * Q) / (phi0 * Q + Pmax)) - Rd
}

# high-light reference PAR used instead of the asymptotic Pmax
q_ref <- 2000

compute_lcp <- function(phi0, Pmax, Rd) {
  if (is.finite(phi0) && is.finite(Pmax) && is.finite(Rd) && (Pmax - Rd) > 0) {
    Rd * Pmax / (phi0 * (Pmax - Rd))
  } else {
    NA_real_
  }
}

compute_fit_summary <- function(fit, dff, q_ref = 2000) {
  cf <- coef(fit)
  phi0 <- cf["phi0"]
  Pmax <- cf["Pmax"]
  Rd   <- cf["Rd"]
  
  data.frame(
    phi0  = as.numeric(phi0),
    Pmax  = as.numeric(Pmax),
    Rd    = as.numeric(Rd),
    LCP   = as.numeric(compute_lcp(phi0, Pmax, Rd)),
    P2000 = as.numeric(light_response(q_ref, phi0, Pmax, Rd)),
    N     = nrow(dff)
  )
}

bootstrap_light_ci <- function(df, n_boot = 500, q_ref = 2000,
                               par_grid = seq(20, 2200, length.out = 60)) {
  dff <- filter_for_fit(df)
  if (nrow(dff) < 10) return(NULL)
  
  boot_rows <- vector("list", n_boot)
  
  for (i in seq_len(n_boot)) {
    sampled <- dff[sample(seq_len(nrow(dff)), replace = TRUE), , drop = FALSE]
    res <- fit_light_curve(sampled)
    
    if (!is.null(res)) {
      boot_rows[[i]] <- compute_fit_summary(res$fit, res$df, q_ref = q_ref)
    }
  }
  
  boot_tbl <- dplyr::bind_rows(boot_rows)
  if (!nrow(boot_tbl)) return(NULL)

  ci <- boot_tbl %>%
    summarise(
      phi0_low  = quantile(phi0,  0.025, na.rm = TRUE),
      phi0_high = quantile(phi0,  0.975, na.rm = TRUE),
      Rd_low    = quantile(Rd,    0.025, na.rm = TRUE),
      Rd_high   = quantile(Rd,    0.975, na.rm = TRUE),
      LCP_low   = quantile(LCP,   0.025, na.rm = TRUE),
      LCP_high  = quantile(LCP,   0.975, na.rm = TRUE),
      P2000_low = quantile(P2000, 0.025, na.rm = TRUE),
      P2000_high = quantile(P2000, 0.975, na.rm = TRUE),
      n_boot_success = n()
    )

  # Confidence band for the CURVE, from the draws already taken above. A
  # co-author could not judge whether the four fitted lines differ, because the
  # figure drew no uncertainty at all; this supplies it without consuming any
  # extra random numbers, so the parameter CIs above are unchanged.
  band <- do.call(rbind, lapply(par_grid, function(Q) {
    y <- light_response(Q, boot_tbl$phi0, boot_tbl$Pmax, boot_tbl$Rd)   # positive = uptake
    data.frame(PAR = Q,
               nee_med = -median(y, na.rm = TRUE),                      # NEE sign for plotting
               nee_low = -quantile(y, 0.975, na.rm = TRUE),
               nee_high = -quantile(y, 0.025, na.rm = TRUE))
  }))
  rownames(band) <- NULL

  list(summary = ci, band = band, draws = boot_tbl[, c("phi0", "Pmax", "Rd")])
}

filter_for_fit <- function(df) {
  if (!all(c("PAR", "NEE_ok") %in% names(df))) return(df[0, , drop = FALSE])
  out <- df[df$PAR > 20 & is.finite(df$PAR) & is.finite(df$NEE_ok), , drop = FALSE]
  out$FC_pos <- -out$NEE_ok
  out
}

fit_light_curve <- function(df) {
  df <- filter_for_fit(df)
  
  # basic sanity so nlsLM doesn't crash immediately
  if (is.null(df) || nrow(df) < 5) return(NULL)
  if (length(unique(df$PAR)) < 3 || length(unique(df$FC_pos)) < 3) return(NULL)
  
  tryCatch({
    fit <- nlsLM(
      FC_pos ~ ((phi0 * Pmax * PAR) / (phi0 * PAR + Pmax)) - Rd,
      data = df,
      start = list(phi0 = 0.03, Pmax = 30, Rd = 5),
      lower = c(phi0 = 0.005, Pmax = 5, Rd = 0),
      upper = c(phi0 = 0.30,  Pmax = 100, Rd = 20),
      control = nls.lm.control(maxiter = 500)
    )
    list(fit = fit, df = df)
  }, error = function(e) NULL)
}

extract_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(idx)) g$grobs[[idx[1]]] else NULL
}

legend_panel_grob <- function(legend_grob, label = NULL) {
  if (is.null(legend_grob)) return(ggplotGrob(blank_panel("")))
  
  if (is.null(label) || identical(label, "")) {
    return(
      arrangeGrob(
        grobs = list(legend_grob),
        layout_matrix = matrix(1),
        vp = viewport(
          x = unit(0.02, "npc"),
          y = unit(0.98, "npc"),
          just = c("left", "top")
        )
      )
    )
  }
  
  arrangeGrob(
    grobs = list(
      textGrob(
        label,
        x = unit(0.02, "npc"),
        y = unit(0.98, "npc"),
        just = c("left", "top"),
        gp = gpar(fontsize = 12, fontface = "bold")
      ),
      legend_grob
    ),
    layout_matrix = rbind(c(1), c(2)),
    heights = unit.c(unit(1.2, "lines"), unit(1, "null")),
    vp = viewport(
      x = unit(0.02, "npc"),
      y = unit(0.98, "npc"),
      just = c("left", "top")
    )
  )
}


legend_title_for <- function(color_var) {
  if (identical(color_var, "Ta"))  return(expression(italic(T)[plain(a)] ~ "(" * degree*C * ")"))
  if (identical(color_var, "VPD")) return("VPD (kPa)")
  return(color_var)
}


# ───────────────────────────────────────────────────────────────────────────────
# prep data (this csv has: ym, season, ENSO, hour, ...)
df_all <- rename_with_map(df_all, rename_map)
if (!"PAR_corrected_SWin" %in% names(df_all) && "PPFD_IN_1_1_1" %in% names(df_all)) {
  df_all$PAR_corrected_SWin <- suppressWarnings(as.numeric(df_all$PPFD_IN_1_1_1))
}
if (!"PAR" %in% names(df_all) && "PAR_corrected_SWin" %in% names(df_all)) {
  df_all$PAR <- suppressWarnings(as.numeric(df_all$PAR_corrected_SWin))
}
if (!"PAR" %in% names(df_all)) {
  stop("MATLAB-corrected PAR column is missing from light-response input")
}
par_mismatch <- is.finite(df_all$PAR) & is.finite(df_all$PAR_corrected_SWin) & abs(df_all$PAR - df_all$PAR_corrected_SWin) > 1e-8
if (any(par_mismatch, na.rm = TRUE)) {
  stop("canonical PAR does not match PAR_corrected_SWin in light-response input")
}
df_all <- make_hour_numeric(df_all)
df_all <- apply_cutoff_ym(df_all, cutoff_ym)

# detect season + enso columns (this file uses season + ENSO)
season_col_candidates <- c("season","Season","SEASON")
enso_col_candidates   <- c("ENSO","enso","enso_phase","ENSO_phase","phase","Phase")

season_col <- season_col_candidates[season_col_candidates %in% names(df_all)][1]
enso_col   <- enso_col_candidates[enso_col_candidates %in% names(df_all)][1]

if (is.na(season_col) || is.na(enso_col)) {
  stop(
    paste0(
      "could not find season/ENSO columns in ", basename(csv_fp), ".\n",
      "edit season_col_candidates / enso_col_candidates."
    )
  )
}

if (!all(c("PAR", "NEE_ok") %in% names(df_all))) {
  stop("dataset must include columns PAR and NEE_ok (and optionally Ta/VPD/SWC/etc).")
}

df_all <- df_all %>%
  mutate(
    season_norm = norm_season(.data[[season_col]]),
    enso_norm   = norm_enso(.data[[enso_col]])
  ) %>%
  filter(!is.na(season_norm), !is.na(enso_norm))

cat("\nseason×enso counts (sanity check):\n")
print(table(df_all$season_norm, df_all$enso_norm, useNA = "ifany"))

# build groups (note: neutral-wet may truly be 0 in your csv)
make_group <- function(season_val, enso_val) {
  df_all %>%
    filter(
      .data$season_norm == .env$season_val,
      .data$enso_norm   == .env$enso_val
    ) %>%
    select(-season_norm, -enso_norm)
}

df_list <- list(
  wet_el_nino   = make_group("wet", "el_nino"),
  dry_el_nino   = make_group("dry", "el_nino"),
  wet_la_nina   = make_group("wet", "la_nina"),
  dry_la_nina   = make_group("dry", "la_nina"),
  wet_neutral   = make_group("wet", "neutral"),
  dry_neutral   = make_group("dry", "neutral")
)

cat("\nrows per group:\n")
print(sapply(df_list, nrow))

# ───────────────────────────────────────────────────────────────────────────────
# plotting engines (pooled only)

make_panel <- function(df,
                       title_text,
                       color_var,
                       global_x_range,
                       global_y_range,
                       color_levels,
                       color_range,
                       band = NULL,
                       show_legend = FALSE,
                       show_y_ticks = TRUE,
                       show_x_ticks = TRUE) {
  df <- filter_for_fit(df)
  if (!nrow(df)) return(blank_panel(title_text))
  
  res <- fit_light_curve(df)
  if (is.null(res)) return(blank_panel(title_text))
  
  fit <- res$fit
  dff <- res$df
  
  params <- coef(fit)
  phi0 <- params["phi0"]
  Pmax <- params["Pmax"]
  Rd   <- params["Rd"]
  N    <- nrow(dff)
  
  LCP <- if (is.finite(phi0) && is.finite(Pmax) && is.finite(Rd) && (Pmax - Rd) > 0) {
    Rd * Pmax / (phi0 * (Pmax - Rd))
  } else {
    NA_real_
  }
  
  P2000 <- light_response(q_ref, phi0, Pmax, Rd)
  
  Q_vals <- seq(global_x_range[1], global_x_range[2], length.out = 200)
  fit_df <- data.frame(PAR = Q_vals, FC_pos = light_response(Q_vals, phi0, Pmax, Rd))

  # Display the paper's NEE sign convention (NEE > 0 = release, NEE < 0 = uptake).
  # The curve is fitted on uptake (FC_pos = -NEE_ok); negate only for plotting so
  # the y-axis shows NEE. Fit parameters (phi0, P2000, Rd, LCP) are unchanged.
  dff$FC_pos    <- -dff$FC_pos
  fit_df$FC_pos <- -fit_df$FC_pos

  # 95% bootstrap band for the fitted curve, already on the NEE sign convention.
  # Drawn under the points so it does not hide data.
  band_layer <- if (!is.null(band) && nrow(band)) {
    list(geom_ribbon(data = band,
                     aes(x = PAR, ymin = nee_low, ymax = nee_high),
                     inherit.aes = FALSE, fill = "grey35", alpha = 0.22))
  } else NULL

  if (!color_var %in% names(dff)) {
    p <- ggplot(dff, aes(x = PAR, y = FC_pos)) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
      band_layer +
      geom_point(alpha = 0.6, na.rm = TRUE) +
      geom_line(data = fit_df, aes(x = PAR, y = FC_pos), color = "black", linewidth = 1.1, inherit.aes = FALSE) +
      guides(color = "none")
  } else if (identical(color_var, "friagem")) {
    dff$friagem <- as.logical(dff$friagem)
    p <- ggplot(dff, aes(x = PAR, y = FC_pos, color = friagem)) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
      band_layer +
      geom_point(alpha = 0.6, na.rm = TRUE) +
      geom_line(data = fit_df, aes(x = PAR, y = FC_pos), color = "black", linewidth = 1.1, inherit.aes = FALSE) +
      scale_color_manual(name = "friagem", values = c("TRUE" = "blue", "FALSE" = "red"))
    if (!show_legend) p <- p + guides(color = "none")
  } else {
    p <- ggplot(dff, aes(x = PAR, y = FC_pos, color = .data[[color_var]])) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
      band_layer +
      geom_point(alpha = 0.6, na.rm = TRUE) +
      geom_line(data = fit_df, aes(x = PAR, y = FC_pos), color = "black", linewidth = 1.1, inherit.aes = FALSE)
    
    cd <- dff[[color_var]]
    is_discrete <- is.factor(cd) || is.character(cd) || is.logical(cd)
    
    if (is_discrete) {
      if (!is.null(color_levels)) {
        p <- p + scale_color_viridis_d(option = "C", name = legend_title_for(color_var), limits = color_levels)
      } else {
        p <- p + scale_color_viridis_d(option = "C", name = legend_title_for(color_var))
      }
    } else {
      if (!is.null(color_range) && all(is.finite(color_range))) {
        p <- p + scale_color_viridis_c(option = "C", name = legend_title_for(color_var), limits = color_range)
      } else {
        p <- p + scale_color_viridis_c(option = "C", name = legend_title_for(color_var))
      }
    }
    
    if (!show_legend) p <- p + guides(color = "none")
  }
  
  if (is.finite(LCP)) p <- p + geom_vline(xintercept = LCP, linetype = "dotted")
  if (all(is.finite(global_x_range)) && all(is.finite(global_y_range))) {
    p <- p + coord_cartesian(
      xlim = global_x_range,
      ylim = global_y_range
    )
  }
  
  # parameter subtitle as plotmath: italic symbols with roman subscripts (phi0, P2000, Rd; LCP roman)
  line1 <- bquote(
    italic(phi)[0] == .(sprintf("%.3f", phi0)) * "," ~
      italic(P)[2000] == .(sprintf("%.1f", P2000)) * "," ~
      italic(R)[d] == .(sprintf("%.2f", Rd))
  )
  line2 <- if (is.finite(LCP)) {
    bquote("LCP" %~~% .(sprintf("%.0f", LCP)) * "," ~ italic(N) == .(N))
  } else {
    bquote(italic(N) == .(N))
  }
  # single line: the four parameters plus LCP and N all fit on one row
  subtitle_expr <- bquote(.(line1) * "," ~ .(line2))

  p <- p + labs(
    title = title_text,
    subtitle = subtitle_expr,
    x = NULL,
    y = NULL
  ) +
    theme_minimal(base_size = fig_base_size) +
    theme(
      plot.title = element_text(hjust = 0, size = fig_title_size, face = "bold"),
      plot.subtitle = element_text(
        hjust = 0,
        size = fig_subtitle_size,
        lineheight = 1.05
      ),
      axis.text = element_text(size = fig_axis_text_size),
      legend.position = if (show_legend) "right" else "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  
  if (!show_y_ticks) {
    p <- p + theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  if (!show_x_ticks) {
    p <- p + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  p
}

# bootstrap curve bands, filled in by create_plot_list_all_data() when do_ci is
# on. Kept at script scope so the overlay panel built after the grid can reuse
# the same draws rather than bootstrapping a second time.
curve_bands <- list()
curve_draws <- list()

create_plot_list_all_data <- function(color_var, do_ci = FALSE) {
  if (do_ci) { curve_bands <<- list(); curve_draws <<- list() }
  pooled <- bind_rows(lapply(df_list, filter_for_fit), .id = "GroupID")

  global_x_range <- c(0, 2500)
  
  global_y_range <- if (nrow(pooled)) {
    range(-pooled$FC_pos, na.rm = TRUE)  # NEE convention: uptake plotted negative
  } else {
    c(NA_real_, NA_real_)
  }
  
  color_levels <- NULL
  color_range  <- NULL
  if (nrow(pooled) && color_var %in% names(pooled)) {
    cd <- pooled[[color_var]]
    is_discrete <- is.factor(cd) || is.character(cd) || is.logical(cd)
    if (is_discrete) {
      vals <- unique(cd[!is.na(cd)])
      if (length(vals)) color_levels <- sort(vals)
    } else {
      vals <- cd[is.finite(cd)]
      if (length(vals)) color_range <- range(vals, na.rm = TRUE)
    }
  }
  
  plot_list <- list()
  param_table <- data.frame(
    Group = character(),
    phi0 = numeric(),
    Pmax = numeric(),
    Rd = numeric(),
    LCP = numeric(),
    P2000 = numeric(),
    N = integer(),
    phi0_low = numeric(),
    phi0_high = numeric(),
    Rd_low = numeric(),
    Rd_high = numeric(),
    LCP_low = numeric(),
    LCP_high = numeric(),
    P2000_low = numeric(),
    P2000_high = numeric(),
    n_boot_success = integer()
  )
  
  i <- 1
  for (group_name in panel_order_display) {
    d <- df_list[[group_name]]
    title_text <- paste(panel_labels_display[i], pretty_names[[group_name]])
    i <- i + 1
    
    # y ticks only on left column; x ticks only on bottom row
    is_left_col   <- group_name %in% c("dry_el_nino", "wet_el_nino")
    is_bottom_row <- group_name %in% c("wet_el_nino", "wet_la_nina")

    # bootstrap first, so the panel can draw its confidence ribbon. Moving this
    # ahead of make_panel is RNG-neutral: make_panel only refits deterministically.
    if (do_ci) {
      boot_out   <- bootstrap_light_ci(d, n_boot = 500, q_ref = q_ref)
      ci_summary <- boot_out$summary
      if (!is.null(boot_out$band)) {
        curve_bands[[pretty_names[[group_name]]]] <<- boot_out$band
        curve_draws[[pretty_names[[group_name]]]] <<- boot_out$draws
      }
    } else {
      boot_out <- NULL; ci_summary <- NULL
    }

    p <- make_panel(
      band = boot_out$band,
      df = d,
      title_text = title_text,
      color_var = color_var,
      global_x_range = global_x_range,
      global_y_range = global_y_range,
      color_levels = color_levels,
      color_range = color_range,
      show_legend = FALSE,
      show_y_ticks = is_left_col,
      show_x_ticks = is_bottom_row
    )
    
    plot_list[[group_name]] <- p
    
    res <- fit_light_curve(d)
    if (!is.null(res)) {
      cf <- coef(res$fit)
      phi0 <- cf["phi0"]
      Pmax <- cf["Pmax"]
      Rd   <- cf["Rd"]
      N    <- nrow(res$df)
      LCP  <- if (is.finite(phi0) && is.finite(Pmax) && is.finite(Rd) && (Pmax - Rd) > 0) {
        Rd * Pmax / (phi0 * (Pmax - Rd))
      } else {
        NA_real_
      }
      
      fit_summary <- compute_fit_summary(res$fit, res$df, q_ref = q_ref)
      
      # ci_summary was computed above, before make_panel, so the panel could draw
      # its ribbon from the same draws

      if (is.null(ci_summary)) {
        ci_summary <- data.frame(
          phi0_low = NA_real_,
          phi0_high = NA_real_,
          Rd_low = NA_real_,
          Rd_high = NA_real_,
          LCP_low = NA_real_,
          LCP_high = NA_real_,
          P2000_low = NA_real_,
          P2000_high = NA_real_,
          n_boot_success = NA_integer_
        )
      }
      
      param_table <- rbind(
        param_table,
        cbind(
          data.frame(Group = pretty_names[[group_name]]),
          fit_summary,
          ci_summary
        )
      )
    }
  }
  
  # pick a legend donor that exists and has rows
  legend_donor <- NULL
  for (k in panel_order_display) {
    if (!is.null(df_list[[k]]) && nrow(df_list[[k]]) > 0) { legend_donor <- k; break }
  }
  if (is.null(legend_donor)) legend_donor <- panel_order_display[1]
  
  meta <- list(
    color_var = color_var,
    global_x_range = global_x_range,
    global_y_range = global_y_range,
    color_levels = color_levels,
    color_range = color_range,
    df_for_legend = df_list[[legend_donor]]
  )
  
  list(plots = plot_list, params = param_table, meta = meta)
}

# ───────────────────────────────────────────────────────────────────────────────
# save helper (2×3, legend in bottom-right)
# - no title
# - shared X label at bottom
# - shared Y label at left (rotated)

save_grid <- function(results, filename_prefix, folder = graphs_path, paper = FALSE) {
  m <- results$meta
  
  legend_plot <- make_panel(
    df = m$df_for_legend,
    title_text = legend_panel_label,
    color_var = m$color_var,
    global_x_range = m$global_x_range,
    global_y_range = m$global_y_range,
    color_levels = m$color_levels,
    color_range = m$color_range,
    show_legend = TRUE,
    show_y_ticks = FALSE,
    show_x_ticks = TRUE
  ) +
    theme(
      plot.subtitle = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    guides(
      color = guide_colorbar(
        barheight = unit(6.5, "cm"),
        barwidth  = unit(0.8, "cm"),
        title.position = "top",
        title.hjust = 0
      )
    ) +
    theme(
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 13)
    )
  
  lg <- extract_legend(legend_plot)
  lg_panel <- legend_panel_grob(lg, label = legend_panel_label)
  
  grobs_in_order <- list(
    results$plots[["dry_el_nino"]],
    results$plots[["dry_la_nina"]],
    results$plots[["dry_neutral"]],
    results$plots[["wet_el_nino"]],
    results$plots[["wet_la_nina"]],
    lg_panel
  )
  
  shared_x <- textGrob(
    "PAR (µmol m⁻² s⁻¹)", 
    gp = gpar(fontsize = fig_shared_axis_size), 
    x = unit(1/3, "npc"), just = "center" 
  )
  
  shared_y <- textGrob(
    expression(NEE ~ (mu*mol ~ CO[2] ~ m^{-2} ~ s^{-1})),
    gp = gpar(fontsize = fig_shared_axis_size), rot = 90, 
    x = unit(0.5, "npc"), just = "center"
  )
  
  core_grid <- arrangeGrob(
    grobs = grobs_in_order,
    ncol = 3,
    bottom = shared_x
  )
  
  grid_plot <- arrangeGrob(
    shared_y,
    core_grid,
    ncol = 2,
    widths = unit.c(unit(2.6, "lines"), unit(1, "null"))
  )
  
  fname <- if (paper) {
    file.path(paper_figures_path, paste0("FigS_", filename_prefix, "_all_data.png"))
  } else {
    file.path(folder, paste0("light_response_", filename_prefix, ".png"))
  }
  
  ggsave(filename = fname, plot = grid_plot, width = 15, height = 13.5, dpi = 300, bg = "white")
  cat("\nsaved:", fname, "\n")
}

# ───────────────────────────────────────────────────────────────────────────────
# save helper for combined Ta / VPD grid (2 rows × 3 columns)
# Row 1: Ta (El Niño dry, La Niña dry, Legend)
# Row 2: VPD (El Niño wet, La Niña wet, Legend)
# ───────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────
# overlay panel: every fitted curve on shared axes
#
# A co-author could not compare the fitted lines because each sits in its own
# subplot and none carried a confidence band. This puts all five regimes on one
# set of axes with their bootstrap bands. ENSO takes colour and season takes
# linetype, the same encoding as Figure 2b.
#
# The inset matters: the within-season ENSO pairs separate only at LOW light,
# where Rd dominates, while the seasonal pairs separate only at high light. On
# the full-range axes the low-light end is compressed into almost nothing.

make_overlay_panel <- function(bands, x_range = c(0, 2500)) {
  if (!length(bands)) return(blank_panel(""))
  bd <- dplyr::bind_rows(lapply(names(bands), function(nm) {
    b <- bands[[nm]]; b$Group <- nm; b
  }))
  bd <- bd %>%
    mutate(
      ENSO = dplyr::case_when(grepl("La Ni", Group) ~ "La Nina",
                              grepl("El Ni", Group) ~ "El Nino",
                              TRUE ~ "Neutral"),
      Season = ifelse(grepl("wet", Group), "Wet", "Dry")
    )

  base <- function(dat, xlim, show_key) {
    ggplot(dat, aes(x = PAR, group = Group)) +
      geom_hline(yintercept = 0, colour = "black", linetype = "dashed", linewidth = 0.9) +
      geom_ribbon(aes(ymin = nee_low, ymax = nee_high, fill = ENSO),
                  alpha = 0.16, colour = NA) +
      geom_line(aes(y = nee_med, colour = ENSO, linetype = Season), linewidth = 1.15) +
      scale_colour_manual(values = enso_cols, name = "ENSO",
                          labels = c("El Nino" = "El Niño", "La Nina" = "La Niña")) +
      scale_fill_manual(values = enso_cols, guide = "none") +
      scale_linetype_manual(values = c("Wet" = "solid", "Dry" = "dashed"), name = "Season") +
      coord_cartesian(xlim = xlim) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = fig_base_size) +
      theme(panel.grid.minor = element_blank(),
            legend.position = if (show_key) c(0.985, 0.03) else "none",
            legend.justification = c(1, 0), legend.box = "horizontal",
            legend.background = element_rect(fill = scales::alpha("white", 0.75), colour = "grey60"),
            legend.margin = margin(2, 4, 2, 4),
            legend.title = element_text(size = fig_legend_title_size - 4, face = "bold"),
            legend.text = element_text(size = fig_legend_text_size - 3))
  }

  main <- base(bd, x_range, TRUE)
  inset <- base(bd %>% filter(PAR <= 500), c(0, 500), FALSE) +
    labs(title = "low light") +
    theme(plot.title = element_text(size = fig_legend_text_size - 3, hjust = 0.02,
                                    margin = margin(b = 1)),
          axis.text = element_text(size = fig_axis_text_size - 5),
          plot.background = element_rect(fill = "white", colour = "grey60"),
          plot.margin = margin(1, 2, 1, 1))

  # the upper-right quadrant is empty (all curves are well below zero above
  # ~1200 umol), and the key sits bottom-right, so the inset goes there. Keep
  # ymax inside the data range or annotation_custom is clipped by the panel.
  yr <- range(c(bd$nee_low, bd$nee_high), na.rm = TRUE)
  main + annotation_custom(ggplotGrob(inset),
                           xmin = 1120, xmax = 2480,
                           ymin = yr[1] + 0.60 * diff(yr),
                           ymax = yr[2] - 0.01 * diff(yr))
}

save_combined_grid <- function(res_ta, res_vpd, filename_prefix, folder = graphs_path, paper = FALSE) {
  
  # --- 1. Create Legend Panel for Ta (Row 1, Col 3) ---
  m_ta <- res_ta$meta
  legend_plot_ta <- make_panel(
    df = m_ta$df_for_legend,
    title_text = "",
    color_var = m_ta$color_var,
    global_x_range = m_ta$global_x_range,
    global_y_range = m_ta$global_y_range,
    color_levels = m_ta$color_levels,
    color_range = m_ta$color_range,
    show_legend = TRUE,
    show_y_ticks = FALSE,
    show_x_ticks = FALSE # No x-ticks for top row legend
  ) +
    theme(
      plot.subtitle = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(), panel.grid = element_blank()
    ) +
    guides(
      color = guide_colorbar(
        barheight = unit(5, "cm"), barwidth = unit(0.8, "cm"), 
        title.position = "top", title.hjust = 0
      )
    ) +
    theme(
      legend.title = element_text(size = fig_legend_title_size, face = "bold"),
      legend.text  = element_text(size = fig_legend_text_size)
    )
  
  lg_ta <- extract_legend(legend_plot_ta)
  lg_panel_ta <- legend_panel_grob(lg_ta)
  
  # --- 2. Create Legend Panel for VPD (Row 2, Col 3) ---
  m_vpd <- res_vpd$meta
  legend_plot_vpd <- make_panel(
    df = m_vpd$df_for_legend,
    title_text = "",
    color_var = m_vpd$color_var,
    global_x_range = m_vpd$global_x_range,
    global_y_range = m_vpd$global_y_range,
    color_levels = m_vpd$color_levels,
    color_range = m_vpd$color_range,
    show_legend = TRUE,
    show_y_ticks = FALSE,
    show_x_ticks = TRUE # Keep alignment structure for bottom row
  ) +
    theme(
      plot.subtitle = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(), panel.grid = element_blank()
    ) +
    guides(
      color = guide_colorbar(
        barheight = unit(5, "cm"), barwidth = unit(0.8, "cm"), 
        title.position = "top", title.hjust = 0
      )
    ) +
    theme(
      legend.title = element_text(size = fig_legend_title_size, face = "bold"),
      legend.text  = element_text(size = fig_legend_text_size)
    )
  
  lg_vpd <- extract_legend(legend_plot_vpd)
  lg_panel_vpd <- legend_panel_grob(lg_vpd)
  
  # --- 3. Arrange the panels ---
  # (a) is the overlay comparing every fitted curve; the four scatter panels
  # follow as (b)-(e)
  p_overlay <- make_overlay_panel(curve_bands) +
    labs(title = "(a) all regimes, fitted curves with 95% bootstrap bands")

  p_a <- res_ta$plots[["dry_el_nino"]] +
    labs(title = paste("(b)", pretty_names[["dry_el_nino"]]))

  p_b <- res_ta$plots[["dry_la_nina"]] +
    labs(title = paste("(c)", pretty_names[["dry_la_nina"]])) +
    # inset the single Ta colourbar into panel (b)'s empty top-right corner
    guides(color = guide_colorbar(
      barheight = unit(2.1, "cm"), barwidth = unit(0.33, "cm"),
      title.position = "top", title.hjust = 0
    )) +
    theme(
      legend.position = c(0.99, 0.97),
      legend.justification = c(1, 1),
      legend.title = element_text(size = fig_legend_title_size - 1, face = "bold"),
      legend.text  = element_text(size = fig_legend_text_size - 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.6), color = NA),
      legend.margin = margin(2, 3, 2, 3),
      legend.key = element_blank()
    )
  
  p_c <- res_ta$plots[["wet_el_nino"]] +
    labs(title = paste("(d)", pretty_names[["wet_el_nino"]]))

  p_d <- res_ta$plots[["wet_la_nina"]] +
    labs(title = paste("(e)", pretty_names[["wet_la_nina"]]))

  # the four scatter panels are Ta-coloured; the Ta colourbar is inset in panel
  # (c), so they fill a 2x2 grid with no separate legend column
  grobs_in_order <- list(p_overlay, p_a, p_b, p_c, p_d)
  
  # centered under the two data-panel columns
  shared_x <- textGrob(
    "PAR (µmol m⁻² s⁻¹)",
    gp = gpar(fontsize = fig_shared_axis_size),
    x = unit(0.5, "npc"), just = "center"
  )
  
  shared_y <- textGrob(
    expression(NEE ~ (mu*mol ~ CO[2] ~ m^{-2} ~ s^{-1})),
    gp = gpar(fontsize = fig_shared_axis_size), rot = 90, 
    x = unit(0.5, "npc"), just = "center"
  )
  
  # overlay spans the full width on top, scatter panels in a 2x2 beneath
  core_grid <- arrangeGrob(
    grobs = grobs_in_order,
    layout_matrix = rbind(c(1, 1),
                          c(2, 3),
                          c(4, 5)),
    heights = unit(c(1.15, 1, 1), "null"),
    bottom = shared_x
  )
  
  grid_plot <- arrangeGrob(
    shared_y, core_grid,
    ncol = 2,
    widths = unit.c(unit(2.6, "lines"), unit(1, "null"))
  )
  
  # --- 4. Save Output ---
  fname <- if (paper) {
    file.path(paper_figures_path, paste0("Fig4_", filename_prefix, "_combined.jpg"))
  } else {
    file.path(folder, paste0("light_response_", filename_prefix, "_combined.jpg"))
  }
  
  ggsave(filename = fname, plot = grid_plot, width = 15, height = 13.5, dpi = 300, bg = "white")
  cat("\nsaved combined figure:", fname, "\n")
}

# ───────────────────────────────────────────────────────────────────────────────
# run: pooled (all-data) figures
set.seed(123)
res_hour <- create_plot_list_all_data("Hour", do_ci = TRUE)
save_grid(res_hour, "light_response_Hour", paper = TRUE)

# do_ci = TRUE so Figure 4 can draw confidence ribbons and the overlay panel.
# res_hour above is already computed, so its draws -- and therefore Table S6 --
# are untouched by the extra resampling here.
res_ta <- create_plot_list_all_data("Ta", do_ci = TRUE)
save_grid(res_ta, "light_response_Ta", paper = TRUE)

res_vpd <- create_plot_list_all_data("VPD")
save_grid(res_vpd, "light_response_VPD", paper = TRUE)

save_combined_grid(res_ta, res_vpd, "Ta_VPD", paper = TRUE)

# ───────────────────────────────────────────────────────────────────────────────
# where do the fitted curves actually separate?
#
# The co-author's second point -- "they look very similar to me, especially
# considering the CI" -- is a claim that can be tested rather than argued with.
# For each pair, the fraction of bootstrap draws in which one curve sits below
# the other (more uptake) at each light level.

if (length(curve_draws) >= 2) {
  par_grid <- seq(20, 2200, length.out = 60)
  curve_of <- function(nm) {
    d <- curve_draws[[nm]]
    t(vapply(seq_len(nrow(d)),
             function(i) -light_response(par_grid, d$phi0[i], d$Pmax[i], d$Rd[i]),
             numeric(length(par_grid))))          # rows = draws, cols = PAR; NEE sign
  }
  mats <- lapply(setNames(names(curve_draws), names(curve_draws)), curve_of)

  bands_out <- dplyr::bind_rows(lapply(names(curve_bands), function(nm) {
    b <- curve_bands[[nm]]; b$Group <- nm; b
  }))

  nm <- names(mats)
  sep_rows <- list()
  for (a in seq_along(nm)) for (b in seq_len(a - 1)) {
    A <- mats[[nm[a]]]; B <- mats[[nm[b]]]
    n <- min(nrow(A), nrow(B))
    p <- colMeans(A[seq_len(n), , drop = FALSE] < B[seq_len(n), , drop = FALSE])
    sep <- p >= 0.975 | p <= 0.025
    sep_rows[[length(sep_rows) + 1]] <- data.frame(
      A = nm[a], B = nm[b],
      p_A_more_uptake_at_1500 = p[which.min(abs(par_grid - 1500))],
      separated_from = if (any(sep)) min(par_grid[sep]) else NA_real_,
      separated_to   = if (any(sep)) max(par_grid[sep]) else NA_real_,
      n_grid_separated = sum(sep)
    )
  }
  sep_tbl <- dplyr::bind_rows(sep_rows)

  readr::write_csv(bands_out, file.path(output_path, "table_light_response_curve_bands.csv"))
  readr::write_csv(sep_tbl,  file.path(output_path, "table_light_response_curve_separation.csv"))
  cat("\nsaved curve bands and pairwise separation tables\n")

  cat("\n=== NEE at PAR 1500 (umol CO2 m-2 s-1), 95% bootstrap interval ===\n")
  j <- which.min(abs(par_grid - 1500))
  print(as.data.frame(bands_out %>%
    group_by(Group) %>% slice(which.min(abs(PAR - 1500))) %>% ungroup() %>%
    transmute(Group, NEE = sprintf("%.1f [%.1f, %.1f]", nee_med, nee_low, nee_high))),
    row.names = FALSE)

  cat("\n=== where each pair separates at 95% (umol m-2 s-1) ===\n")
  print(as.data.frame(sep_tbl %>% transmute(
    A, B, p_at_1500 = round(p_A_more_uptake_at_1500, 2),
    separated = ifelse(is.na(separated_from), "nowhere",
                       sprintf("%.0f-%.0f", separated_from, separated_to)))),
    row.names = FALSE)
}

# pooled summary table: mean environment + fit parameters, now including P2000
env_means <- lapply(df_list, function(df) {
  df_filt <- df[df$PAR > 20 & is.finite(df$PAR), , drop = FALSE]
  
  data.frame(
    Ta    = if ("Ta"    %in% names(df_filt)) mean(df_filt$Ta,    na.rm = TRUE) else NA_real_,
    VPD   = if ("VPD"   %in% names(df_filt)) mean(df_filt$VPD,   na.rm = TRUE) else NA_real_,
    SWC   = if ("SWC"   %in% names(df_filt)) mean(df_filt$SWC,   na.rm = TRUE) else NA_real_,
    Ustar = if ("Ustar" %in% names(df_filt)) mean(df_filt$Ustar, na.rm = TRUE) else NA_real_,
    PAR   = if ("PAR"   %in% names(df_filt)) mean(df_filt$PAR,   na.rm = TRUE) else NA_real_
  )
})

env_means_df <- bind_rows(env_means, .id = "GroupKey")
env_means_df$Group <- pretty_names[env_means_df$GroupKey]

param_df <- res_hour$params

summary_table <- left_join(
  env_means_df %>% select(-GroupKey),
  param_df,
  by = "Group"
)

summary_csv <- file.path(output_path, "light_response_mean_environmental_and_fit_params.csv")
write.csv(summary_table, file = summary_csv, row.names = FALSE)

cat("\nsaved pooled env+params table with P2000 to:", summary_csv, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# pooled params csv (figure-matching) using hour-colored pooled params
params_fig_exact <- res_hour$params %>%
  mutate(
    Season = ifelse(grepl("wet", Group, ignore.case = TRUE), "Wet", "Dry"),
    ENSO = dplyr::case_when(
      grepl("La Niña", Group, fixed = TRUE) ~ "La Niña",
      grepl("El Niño", Group, fixed = TRUE) ~ "El Niño",
      grepl("Neutral", Group) ~ "Neutral",
      TRUE ~ NA_character_
    ),
    phi0 = round(phi0, 3),
    Pmax = round(Pmax, 1),
    Rd   = round(Rd, 2),
    LCP  = round(LCP, 0)
  ) %>%
  select(Season, ENSO, phi0, Pmax, Rd, LCP, N)

params_csv <- file.path(output_path, "season_enso_means_phi0_Pmax_Rd_LCP.csv")
write.csv(params_fig_exact, params_csv, row.names = FALSE)
cat("\nsaved pooled params to:", params_csv, "\n")

# ───────────────────────────────────────────────────────────────────────────────
format_est_ci <- function(est, low, high, digits = 1) {
  ifelse(
    is.na(low) | is.na(high),
    sprintf(paste0("%.", digits, "f"), est),
    paste0(
      sprintf(paste0("%.", digits, "f"), est),
      " [",
      sprintf(paste0("%.", digits, "f"), low),
      ", ",
      sprintf(paste0("%.", digits, "f"), high),
      "]"
    )
  )
}

params_fig_exact <- res_hour$params %>%
  mutate(
    Season = ifelse(grepl("wet", Group, ignore.case = TRUE), "Wet", "Dry"),
    ENSO = dplyr::case_when(
      grepl("La Niña", Group, fixed = TRUE) ~ "La Niña",
      grepl("El Niño", Group, fixed = TRUE) ~ "El Niño",
      grepl("Neutral", Group) ~ "Neutral",
      TRUE ~ NA_character_
    ),
    P2000_CI = format_est_ci(P2000, P2000_low, P2000_high, digits = 1),
    Rd_CI    = format_est_ci(Rd,    Rd_low,    Rd_high,    digits = 2),
    phi0_CI  = format_est_ci(phi0,  phi0_low,  phi0_high,  digits = 3),
    LCP_CI   = format_est_ci(LCP,   LCP_low,   LCP_high,   digits = 0)
  ) %>%
  select(
    Season,
    ENSO,
    P2000_CI,
    Rd_CI,
    phi0_CI,
    LCP_CI,
    N,
    n_boot_success
  )

params_csv <- file.path(output_path, "season_enso_light_response_P2000_bootstrap_CI.csv")
write.csv(params_fig_exact, params_csv, row.names = FALSE)
cat("\nsaved pooled params with P2000 and bootstrap CI to:", params_csv, "\n")
