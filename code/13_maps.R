# ───────────────────────────────────────────────────────────────────────────────
# worldclim maps + in-map inset of Tambopata seasonal means (Wet vs Dry)
# outputs 2 separate figures:
#   (A) annual precipitation map (capped at ≥ 3000 mm/yr) + inset
#   (B) dry-season length map (months P < 100 mm/month; capped at ≥ 6) + inset
#
# updates in this version:
#   - legend is horizontal and placed at top-right
#   - ocean is blank/white (mask raster to South America after derivation)
# ───────────────────────────────────────────────────────────────────────────────

# packages
library(terra)
library(sf)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(lubridate)
library(data.table)
library(tidyr)
library(grid)   # unit()

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
dir.create(paper_figures_path, recursive = TRUE, showWarnings = FALSE)
sf::sf_use_s2(FALSE)

# ───────────────────────────────────────────────────────────────────────────────
# settings
# ───────────────────────────────────────────────────────────────────────────────

amazon_bbox <- ext(-82, -35, -20, 13)  # xmin, xmax, ymin, ymax

dry_thresh_mm <- 100
cap_prec_mm_yr <- 3000
cap_drylen_mo  <- 6

# dry -> wet palette:
dry_to_wet_cols <- c(
  "#B10026", "#E31A1C", "#FD8D3C", "#FED976",
  "#C7E9F1", "#41B6C4", "#2C7FB8", "#253494", "#4A1486"
)

# site point
sites <- data.frame(
  name = "Tambopata",
  lon  = -70.25,
  lat  = -12.83
)
sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

# inset placement (relative to panel; tweak if you want)
inset_left   <- 0.58
inset_bottom <- 0.05
inset_right  <- 0.99
inset_top    <- 0.39

# legend placement/size
legend_pos <- c(0.98, 0.98)
legend_barwidth  <- unit(1, "cm")
legend_barheight <- unit(3, "cm")

# borders + land mask
# use all countries intersecting the amazon bbox so french guiana is included

if (requireNamespace("rnaturalearth", quietly = TRUE)) {
  world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
} else {
  message("optional: install.packages('rnaturalearth') for Natural Earth country outlines; using maps::map fallback")
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop("missing country-outline source: install rnaturalearth or maps")
  }
  world_sf <- sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
  world_sf <- sf::st_transform(world_sf, 4326)
}

bbox_sf <- st_as_sfc(
  st_bbox(
    c(
      xmin = xmin(amazon_bbox),
      xmax = xmax(amazon_bbox),
      ymin = ymin(amazon_bbox),
      ymax = ymax(amazon_bbox)
    ),
    crs = st_crs(4326)
  )
)

sa_sf <- st_intersection(world_sf, bbox_sf)
sa_vect <- vect(sa_sf)

# ───────────────────────────────────────────────────────────────────────────────
# worldclim monthly precipitation climatology (12 layers)
# ───────────────────────────────────────────────────────────────────────────────

worldclim_path <- file.path(output_path, "worldclim")
dir.create(worldclim_path, showWarnings = FALSE, recursive = TRUE)

# res=10 is ~10 arc-min (~18 km); try res=5 for finer (bigger download)
prec_files <- list.files(
  file.path(worldclim_path, "climate", "wc2.1_10m"),
  pattern = "^wc2\\.1_10m_prec_[0-9]{2}\\.tif$",
  full.names = TRUE
)
prec_files <- sort(prec_files)
if (length(prec_files) == 12) {
  prec12 <- terra::rast(prec_files)
} else {
  if (!requireNamespace("geodata", quietly = TRUE)) {
    stop(
      "missing cached WorldClim precipitation rasters in: ",
      file.path(worldclim_path, "climate", "wc2.1_10m"),
      "\nInstall geodata or restore the 12 wc2.1_10m_prec_*.tif files."
    )
  }
  prec12 <- geodata::worldclim_global(var = "prec", res = 10, path = worldclim_path)
}
prec12 <- crop(prec12, amazon_bbox)
prec12 <- mask(prec12, sa_vect)  # mask to SA early to reduce junk

# ───────────────────────────────────────────────────────────────────────────────
# derived rasters
# ───────────────────────────────────────────────────────────────────────────────

# annual precipitation (mm/yr)
ann_prec <- app(prec12, fun = sum, na.rm = TRUE)
names(ann_prec) <- "annual_prec_mm"

# cap at 3000 for visualization (values > 3000 shown as top color)
ann_prec_cap <- clamp(ann_prec, lower = 0, upper = cap_prec_mm_yr, values = TRUE)
ann_prec_cap <- mask(ann_prec_cap, sa_vect)  # make ocean NA -> blank/white
names(ann_prec_cap) <- "annual_prec_mm_cap"
ann_df <- as.data.frame(ann_prec_cap, xy = TRUE, na.rm = FALSE)

# dry-season length (# months with P < threshold)
dry_len <- app(prec12, fun = function(x) sum(x < dry_thresh_mm, na.rm = TRUE))
names(dry_len) <- "dry_season_len_months"

# cap at 6 for visualization
dry_len_cap <- clamp(dry_len, lower = 0, upper = cap_drylen_mo, values = TRUE)
dry_len_cap <- mask(dry_len_cap, sa_vect)    # make ocean NA -> blank/white
names(dry_len_cap) <- "dry_season_len_months_cap"
dry_df <- as.data.frame(dry_len_cap, xy = TRUE, na.rm = FALSE)

# ───────────────────────────────────────────────────────────────────────────────
# build the Tambopata inset panel (Wet vs Dry means)
#   uses your "enso_conditions_data.RData" (flux env vars) + IMERG precip
# ───────────────────────────────────────────────────────────────────────────────

# load the seasonal/ENSO-split flux dataframes
# expects objects like: df_wet_el_nino, df_dry_el_nino, ... etc
load(file.path(output_path, "enso_conditions_data.RData"))

# ───────────────────────────────────────────────────────────────────────────────
# cutoff (same as diel script)
cutoff_datetime <- as.POSIXct("2024-10-15 23:59:59", tz = "UTC")
cutoff_date_local <- as.Date("2024-10-15")  # for IMERG daily (local calendar)

apply_cutoff_flux <- function(d, cutoff_dt) {
  if (is.null(d) || !nrow(d) || !"DateTime" %in% names(d)) return(d)
  
  if (!inherits(d$DateTime, "POSIXct")) {
    d$DateTime <- lubridate::ymd_hms(d$DateTime, quiet = TRUE, tz = "UTC")
    if (all(is.na(d$DateTime))) d$DateTime <- as.POSIXct(d$DateTime, tz = "UTC")
  } else {
    attr(d$DateTime, "tzone") <- "UTC"
  }
  
  d %>% dplyr::filter(!is.na(DateTime), DateTime <= cutoff_dt)
}


# helper: safe mean
mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# rename mapping
rename_map <- c(
  "TA_1_1_1"       = "Ta",
  "VPD_kPa"      = "VPD",
  "SWC_1_1_1"      = "SWC",
  "TS_3"           = "Ts"
)

# gather all 6 groups into 1 table, then collapse to Wet vs Dry (multi-annual mean)
df_list <- list(
  wet_el_nino   = df_wet_el_nino,
  dry_el_nino   = df_dry_el_nino,
  wet_la_nina   = df_wet_la_nina,
  dry_la_nina   = df_dry_la_nina,
  wet_neutral   = df_wet_neutral,
  dry_neutral   = df_dry_neutral
)

df_list <- lapply(df_list, function(d) {
  d <- apply_cutoff_flux(d, cutoff_datetime)
  if (!nrow(d)) return(d)
  if (!"PAR_corrected_SWin" %in% names(d) && "PPFD_IN_1_1_1" %in% names(d)) {
    d$PAR_corrected_SWin <- suppressWarnings(as.numeric(d$PPFD_IN_1_1_1))
  }
  if (!"PAR" %in% names(d) && "PAR_corrected_SWin" %in% names(d)) {
    d$PAR <- suppressWarnings(as.numeric(d$PAR_corrected_SWin))
  }
  if (!"PAR" %in% names(d)) {
    stop("MATLAB-corrected PAR column is missing from enso_conditions_data.RData")
  }
  par_mismatch <- is.finite(d$PAR) & is.finite(d$PAR_corrected_SWin) &
    abs(d$PAR - d$PAR_corrected_SWin) > 1e-8
  if (any(par_mismatch, na.rm = TRUE)) {
    stop("canonical PAR does not match PAR_corrected_SWin in enso_conditions_data.RData")
  }
  for (old in names(rename_map)) {
    new <- rename_map[[old]]
    if (old %in% names(d) && !new %in% names(d)) d <- dplyr::rename(d, !!new := !!rlang::sym(old))
  }
  d
})

flux_all <- bind_rows(lapply(names(df_list), function(k) {
  d <- df_list[[k]]
  if (is.null(d) || !nrow(d)) return(NULL)
  d %>%
    mutate(
      Season = ifelse(grepl("^wet", k), "Wet", "Dry"),
      Season = factor(Season, levels = c("Wet","Dry"))
    )
}))

# means from flux variables
season_means_flux <- flux_all %>%
  summarise(
    Ta  = mean_or_na(Ta),
    VPD = mean_or_na(VPD),
    PAR = mean_or_na(PAR),
    SWC = mean_or_na(SWC),
    Ts  = mean_or_na(Ts),
    .by = Season
  )

# precipitation from IMERG (30-min) -> daily totals -> mean mm/day by Season (Wet/Dry)
tz_local <- "America/Lima"

imerg_csv <- list.files(input_folder, pattern = "(?i)imerg.*30min.*\\.csv$", full.names = TRUE)
if (!length(imerg_csv)) stop("no IMERG 30-min CSV found in input_folder")
imerg_csv <- imerg_csv[order(file.info(imerg_csv)$mtime, decreasing = TRUE)][1]
imerg <- data.table::fread(imerg_csv)

utc_col   <- intersect(names(imerg), c("time_utc_iso","time_utc","time","system:time_start"))[1]
local_col <- intersect(names(imerg), c("time_local_minus05","time_local"))[1]

if (!is.na(local_col)) {
  DateTime <- lubridate::parse_date_time(
    imerg[[local_col]],
    orders = c("Y-m-d H:M:S","Y-m-d H:M","Ymd HMS","Ymd HM"),
    tz = tz_local
  )
} else if (!is.na(utc_col)) {
  t_utc <- lubridate::parse_date_time(
    imerg[[utc_col]],
    orders = c("Y-m-d H:M:S","Y-m-d H:M","Ymd HMS","Ymd HM","YmdHMS","YmdHM"),
    tz = "UTC"
  )
  DateTime <- lubridate::with_tz(t_utc, tz_local)
} else stop("could not find a time column (time_utc*/time_local*).")

if ("precip_mm_30min" %in% names(imerg)) {
  pr30 <- as.numeric(imerg$precip_mm_30min)
} else {
  rate_col <- intersect(names(imerg), c("precip_mm_per_hr","precipitation","precipitationCal"))[1]
  if (is.na(rate_col)) stop("no precipitation column found (precip_mm_30min or rate in mm/hr)")
  pr30 <- as.numeric(imerg[[rate_col]]) * 0.5
}

imerg30 <- data.frame(DateTime = DateTime, precip_mm_30min = pr30) %>%
  filter(DateTime >= as.POSIXct("2017-01-01 00:00:00", tz = tz_local),
         DateTime <= as.POSIXct("2024-12-31 23:59:59", tz = tz_local))

pr_daily <- imerg30 %>%
  mutate(Date = as.Date(DateTime)) %>%
  group_by(Date) %>%
  summarise(precip_day_mm = sum(precip_mm_30min, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(Date <= cutoff_date_local) %>%
  mutate(
    Season = ifelse(month(Date) %in% 5:10, "Dry", "Wet"),
    Season = factor(Season, levels = c("Wet","Dry"))
  )



season_means_precip <- pr_daily %>%
  summarise(Precip = mean_or_na(precip_day_mm), .by = Season)

# merge for inset
season_means <- season_means_flux %>%
  left_join(season_means_precip, by = "Season")

# ── interannual variability (SD across years, 2017–2024) for inset error bars ──
# compute per-year seasonal means first, then the SD across years, per Season × variable
flux_year_means <- flux_all %>%
  mutate(year = lubridate::year(DateTime)) %>%
  summarise(
    Ta  = mean_or_na(Ta),
    VPD = mean_or_na(VPD),
    PAR = mean_or_na(PAR),
    SWC = mean_or_na(SWC),
    Ts  = mean_or_na(Ts),
    .by = c(Season, year)
  )
season_sd_flux <- flux_year_means %>%
  summarise(
    Ta  = sd(Ta,  na.rm = TRUE),
    VPD = sd(VPD, na.rm = TRUE),
    PAR = sd(PAR, na.rm = TRUE),
    SWC = sd(SWC, na.rm = TRUE),
    Ts  = sd(Ts,  na.rm = TRUE),
    .by = Season
  )
precip_year_means <- pr_daily %>%
  mutate(year = lubridate::year(Date)) %>%
  summarise(Precip = mean_or_na(precip_day_mm), .by = c(Season, year))
season_sd_precip <- precip_year_means %>%
  summarise(Precip = sd(Precip, na.rm = TRUE), .by = Season)
season_sd <- season_sd_flux %>%
  left_join(season_sd_precip, by = "Season")

# make long tables (mean + interannual SD) for the faceted mini-bars
mean_long <- season_means %>%
  select(Season, Ta, VPD, PAR, SWC, Precip, Ts) %>%
  tidyr::pivot_longer(-Season, names_to = "var_raw", values_to = "value")
sd_long <- season_sd %>%
  select(Season, Ta, VPD, PAR, SWC, Precip, Ts) %>%
  tidyr::pivot_longer(-Season, names_to = "var_raw", values_to = "sd")

season_long <- mean_long %>%
  dplyr::left_join(sd_long, by = c("Season", "var_raw")) %>%
  mutate(
    ymin = value - sd,
    ymax = value + sd,
    # physically non-negative variables: don't let the lower whisker dip below 0
    ymin = ifelse(var_raw %in% c("Precip","PAR","VPD","SWC") & is.finite(ymin) & ymin < 0, 0, ymin),
    var = factor(
      var_raw,
      levels = c("Precip","VPD","SWC","PAR","Ta","Ts"),
      labels = c(
        "Precip~(mm~d^{-1})",
        "VPD~(kPa)",
        "SWC~(cm)",
        "PAR~(mu*mol~m^{-2}~s^{-1})",
        "italic(T)[plain(a)]~(degree*C)",
        "italic(T)[plain(s)]~(degree*C)"
      )
    )
  )


inset_panel <- ggplot(season_long, aes(x = Season, y = value, fill = Season)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.25, linewidth = 0.3,
                color = "black", na.rm = TRUE) +
  facet_wrap(~ var, nrow = 2, scales = "free_y", dir = "h", labeller = label_parsed)+
  scale_fill_manual(values = c(Wet = "gray35", Dry = "gray70")) +
  labs(title = "Tambopata seasonal means (2017–2024)") +
  theme_minimal(base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    plot.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(4, 4, 4, 4)
  )

# ───────────────────────────────────────────────────────────────────────────────
# base theme for maps
# ───────────────────────────────────────────────────────────────────────────────

base_theme_map <- theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    
    panel.grid.minor = element_blank(),
    
    # remove axes completely
    axis.title = element_blank(),
    #axis.text  = element_blank(),
    #axis.ticks = element_blank(),
    
    legend.position = legend_pos,
    legend.justification = c(1, 1),
    legend.direction = "horizontal",
    
    legend.background = element_rect(
      fill = "white",
      color = "grey40",
      linewidth = 0.4
    ),
    legend.key = element_rect(fill = "white"),
    
    plot.title = element_text(face = "bold", margin = margin(b = 6))
  )

# station marker: black star with white halo
layer_station <- list(
  geom_sf(data = sites_sf, shape = 8, size = 4.6, color = "white", stroke = 1.6),
  geom_sf(data = sites_sf, shape = 8, size = 3.6, color = "black", stroke = 0.9)
)

# ───────────────────────────────────────────────────────────────────────────────
# (A) annual precipitation figure + inset
# ───────────────────────────────────────────────────────────────────────────────

p_prec_map <- ggplot() +
  geom_raster(data = ann_df, aes(x = x, y = y, fill = annual_prec_mm_cap)) +
  geom_sf(data = sa_sf, fill = NA, linewidth = 0.6, color = "grey20")+
  layer_station +
  coord_sf(
    xlim = c(xmin(amazon_bbox), xmax(amazon_bbox)),
    ylim = c(ymin(amazon_bbox), ymax(amazon_bbox)),
    expand = FALSE
  ) +
  scale_fill_gradientn(
    name = "Annual precip\n(mm yr⁻¹)",
    colours = dry_to_wet_cols,
    limits = c(0, cap_prec_mm_yr),
    oob = scales::squish,
    breaks = c(0, 500, 1000, 1500, 2000, 2500, cap_prec_mm_yr),
    labels = c("0","500","1,000","1,500","2,000","2,500", paste0(">= ", cap_prec_mm_yr)),
    na.value = "white",
    guide = guide_colorbar(
      direction = "vertical",
      barwidth  = legend_barwidth,
      barheight = legend_barheight,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(title = paste0("Annual precipitation")) +
  base_theme_map

fig_prec <- p_prec_map +
  patchwork::inset_element(
    inset_panel,
    left = inset_left, bottom = inset_bottom, right = inset_right, top = inset_top,
    align_to = "full"
  )

ggsave(
  filename = file.path(paper_figures_path, "amazon_annual_precip_with_inset.png"),
  plot     = fig_prec,
  width    = 8.5,
  height   = 6.0,
  dpi      = 300
)

# ───────────────────────────────────────────────────────────────────────────────
# (B) dry-season length figure + inset
# ───────────────────────────────────────────────────────────────────────────────

p_drylen_map <- ggplot() +
  geom_raster(data = dry_df, aes(x = x, y = y, fill = dry_season_len_months_cap)) +
  geom_sf(data = sa_sf, fill = NA, linewidth = 0.6, color = "grey20")+
  layer_station +
  coord_sf(
    xlim = c(xmin(amazon_bbox), xmax(amazon_bbox)),
    ylim = c(ymin(amazon_bbox), ymax(amazon_bbox)),
    expand = FALSE
  ) +
  scale_fill_gradientn(
    name = "Dry-season length\n(months)",
    colours = rev(dry_to_wet_cols),
    limits  = c(0, cap_drylen_mo),
    oob     = scales::squish,
    breaks  = 0:cap_drylen_mo,
    labels  = c(as.character(0:(cap_drylen_mo - 1)), paste0(">= ", cap_drylen_mo)),
    na.value = "white",
    guide = guide_colorbar(
      direction = "vertical",
      barwidth  = legend_barwidth,
      barheight = legend_barheight,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(title = paste0("Dry-season length")) +
  base_theme_map

fig_drylen <- p_drylen_map +
  patchwork::inset_element(
    inset_panel,
    left = inset_left, bottom = inset_bottom, right = inset_right, top = inset_top,
    align_to = "full"
  )

ggsave(
  filename = file.path(paper_figures_path, "amazon_dryseason_length_with_inset.png"),
  plot     = fig_drylen,
  width    = 8.5,
  height   = 6.0,
  dpi      = 300
)


range(flux_all$DateTime, na.rm = TRUE)
range(pr_daily$Date, na.rm = TRUE)

message("done. wrote:\n  - amazon_annual_precip_with_inset.(png|pdf)\n  - amazon_dryseason_length_with_inset.(png|pdf)")



# ───────────────────────────────────────────────────────────────────────────────
# (C) PanAmazonFlux tower sites + inset
# ───────────────────────────────────────────────────────────────────────────────

panamazon_sites_csv <- file.path(input_folder, "panamazon_flux_sites.csv")
if (!file.exists(panamazon_sites_csv)) {
  stop(
    "missing PanAmazonFlux site snapshot: ", panamazon_sites_csv,
    "\nRun code/update_panamazon_sites.R to regenerate it."
  )
}

# read the local PanAmazonFlux snapshot
flux_sites <- data.table::fread(panamazon_sites_csv, encoding = "UTF-8", data.table = FALSE) %>%
  transmute(
    site_id = as.character(acronym),
    site_name = as.character(name),
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    country = as.character(country),
    status = tolower(as.character(status))
  ) %>%
  mutate(
    is_tambopata = site_id == "TNR",
    status_group = if_else(status == "active", "Active", "Building / reactivating")
  )

if (nrow(flux_sites) < 29) {
  stop("expected at least 29 PanAmazonFlux sites, found ", nrow(flux_sites), ".")
}
if (any(is.na(flux_sites$lon) | is.na(flux_sites$lat))) {
  stop("PanAmazonFlux site snapshot has missing coordinates.")
}
if (any(duplicated(flux_sites$site_id))) {
  stop(
    "PanAmazonFlux site snapshot has duplicate acronyms: ",
    paste(flux_sites$site_id[duplicated(flux_sites$site_id)], collapse = ", ")
  )
}

message(
  "PanAmazonFlux tower sites: ", nrow(flux_sites),
  " total; ", sum(flux_sites$status_group == "Active"), " active plotted; ",
  sum(flux_sites$status_group != "Active"), " not plotted; countries: ",
  paste(sort(unique(flux_sites$country)), collapse = ", ")
)

flux_sites_sf <- st_as_sf(
  flux_sites %>% filter(status_group == "Active"),
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

p_flux_sites_map <- ggplot() +
  geom_raster(data = ann_df, aes(x = x, y = y, fill = annual_prec_mm_cap)) +
  geom_sf(data = sa_sf, fill = NA, linewidth = 0.6, color = "grey20") +
  
  # active PanAmazonFlux sites
  geom_sf(
    data = flux_sites_sf %>% filter(!is_tambopata),
    shape = 21,
    size = 2.2,
    fill = "white",
    color = "black",
    stroke = 0.5
  ) +
  
  # tambopata highlighted: white halo
  geom_sf(
    data = flux_sites_sf %>% filter(is_tambopata),
    shape = 8,
    size = 6.2,
    color = "white",
    stroke = 2.2
  ) +

  # tambopata highlighted: red star (stands out from the white circles / precip map)
  geom_sf(
    data = flux_sites_sf %>% filter(is_tambopata),
    shape = 8,
    size = 4.8,
    color = "red2",
    stroke = 1.4
  ) +

  # tambopata label
  geom_label(
    data = flux_sites_sf %>% filter(is_tambopata),
    aes(x = lon, y = lat),
    label = "Tambopata (PE-TNR)",
    inherit.aes = FALSE,
    nudge_y = -2.2,
    hjust = 0.5,
    size = 2.6,
    fontface = "bold",
    color = "black",
    fill = "white",
    alpha = 0.85,
    label.size = 0.25,
    label.padding = unit(0.12, "lines")
  ) +
  
  coord_sf(
    xlim = c(xmin(amazon_bbox), xmax(amazon_bbox)),
    ylim = c(ymin(amazon_bbox), ymax(amazon_bbox)),
    expand = FALSE
  ) +
  scale_fill_gradientn(
    name = "Annual precip\n(mm yr⁻¹)",
    colours = dry_to_wet_cols,
    limits = c(0, cap_prec_mm_yr),
    oob = scales::squish,
    breaks = c(0, 500, 1000, 1500, 2000, 2500, cap_prec_mm_yr),
    labels = c("0", "500", "1,000", "1,500", "2,000", "2,500", paste0(">= ", cap_prec_mm_yr)),
    na.value = "white",
    guide = guide_colorbar(
      direction = "vertical",
      barwidth = legend_barwidth,
      barheight = legend_barheight,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(title = NULL) +
  base_theme_map

# add the meteorological inset
fig_flux_sites <- p_flux_sites_map +
  patchwork::inset_element(
    inset_panel,
    left = inset_left,
    bottom = inset_bottom,
    right = inset_right,
    top = inset_top,
    align_to = "full"
  )

ggsave(
  filename = file.path(paper_figures_path, "amazon_flux_tower_sites_with_inset.png"),
  plot = fig_flux_sites,
  width = 8.5,
  height = 6.0,
  dpi = 300
)
