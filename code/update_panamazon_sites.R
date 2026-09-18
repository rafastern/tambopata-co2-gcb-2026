#!/usr/bin/env Rscript

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

sites_url <- "https://panamazonflux.org/sites"
out_file <- file.path(input_folder, "panamazon_flux_sites.csv")

fetch_text <- function(url) {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
  paste(readLines(tmp, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

extract_first_match <- function(x, pattern) {
  hit <- regexpr(pattern, x, perl = TRUE)
  if (hit[1] == -1) return(NA_character_)
  regmatches(x, hit)
}

extract_field <- function(object_text, field, numeric = FALSE) {
  pattern <- if (numeric) {
    paste0(field, ":([-]?[0-9]+(?:\\.[0-9]+)?)")
  } else {
    paste0(field, ':("(?:[^"\\\\]|\\\\.)*")')
  }

  hit <- regexpr(pattern, object_text, perl = TRUE)
  if (hit[1] == -1) {
    return(if (numeric) NA_real_ else NA_character_)
  }

  value <- sub(paste0("^", field, ":"), "", regmatches(object_text, hit))
  if (numeric) {
    as.numeric(value)
  } else {
    gsub('^"|"$', "", value)
  }
}

page <- fetch_text(sites_url)
bundle_path <- extract_first_match(page, '/assets/index-[^"]+\\.js')
if (is.na(bundle_path)) {
  stop("Could not find the PanAmazonFlux JavaScript bundle in ", sites_url)
}

bundle_url <- paste0("https://panamazonflux.org", bundle_path)
bundle <- fetch_text(bundle_url)

array_hit <- regexpr("HT=\\[\\{.*?\\}\\],mc=", bundle, perl = TRUE)
if (array_hit[1] == -1) {
  stop("Could not find the embedded PanAmazonFlux site array in ", bundle_url)
}

array_text <- regmatches(bundle, array_hit)
array_text <- sub("^HT=\\[", "", array_text)
array_text <- sub("\\],mc=$", "", array_text)
objects <- strsplit(array_text, "\\},\\{", perl = TRUE)[[1]]
objects <- paste0(
  ifelse(grepl("^\\{", objects), "", "{"),
  objects,
  ifelse(grepl("\\}$", objects), "", "}")
)

sites <- data.frame(
  id = vapply(objects, extract_field, numeric(1), field = "id", numeric = TRUE),
  name = vapply(objects, extract_field, character(1), field = "name"),
  acronym = vapply(objects, extract_field, character(1), field = "acronym"),
  lat = vapply(objects, extract_field, numeric(1), field = "lat", numeric = TRUE),
  lon = vapply(objects, extract_field, numeric(1), field = "lon", numeric = TRUE),
  country = vapply(objects, extract_field, character(1), field = "country"),
  state = vapply(objects, extract_field, character(1), field = "state"),
  vegetation = vapply(objects, extract_field, character(1), field = "vegetation"),
  status = vapply(objects, extract_field, character(1), field = "status"),
  year = vapply(objects, extract_field, character(1), field = "year"),
  pi = vapply(objects, extract_field, character(1), field = "pi"),
  stringsAsFactors = FALSE
)

required <- c("id", "name", "acronym", "lat", "lon", "country", "state", "vegetation", "status")
if (nrow(sites) != 29) {
  stop("Expected 29 PanAmazonFlux sites, found ", nrow(sites), ".")
}
if (anyNA(sites[required])) {
  stop("Missing required PanAmazonFlux site fields.")
}
if (any(duplicated(sites$acronym))) {
  stop("Duplicate PanAmazonFlux acronyms: ", paste(sites$acronym[duplicated(sites$acronym)], collapse = ", "))
}
if (any(!is.finite(sites$lat) | !is.finite(sites$lon) | sites$lat < -90 | sites$lat > 90 | sites$lon < -180 | sites$lon > 180)) {
  stop("Invalid PanAmazonFlux coordinates.")
}

sites <- sites[order(sites$id), ]
utils::write.csv(sites, out_file, row.names = FALSE, fileEncoding = "UTF-8")

message("wrote ", nrow(sites), " PanAmazonFlux sites to ", out_file)
