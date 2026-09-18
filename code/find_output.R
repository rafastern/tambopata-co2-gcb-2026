#!/usr/bin/env Rscript

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

args <- commandArgs(trailingOnly = TRUE)
query <- if (length(args)) paste(args, collapse = " ") else ""

script_files <- list.files(code_path, pattern = "\\.R$", full.names = TRUE)
script_files <- script_files[basename(script_files) != "find_output.R"]

output_ext <- "\\.(png|pdf|csv|RData|rds|jpg|jpeg|tif|tiff)$"

extract_quoted_strings <- function(line) {
  matches <- gregexpr("\"([^\"\\\\]|\\\\.)*\"|'([^'\\\\]|\\\\.)*'", line, perl = TRUE)[[1]]
  if (matches[1] < 0) return(character())

  raw <- regmatches(line, list(matches))[[1]]
  substring(raw, 2, nchar(raw) - 1)
}

is_placeholder_artifact <- function(artifact) {
  basename_artifact <- basename(artifact)
  grepl("^\\.[A-Za-z0-9]+$", basename_artifact) ||
    grepl("%[[:alnum:]\\.]*[A-Za-z]", artifact)
}

rows <- list()

for (script in script_files) {
  lines <- readLines(script, warn = FALSE)

  for (i in seq_along(lines)) {
    strings <- extract_quoted_strings(lines[[i]])
    artifacts <- strings[grepl(output_ext, strings, ignore.case = TRUE)]
    artifacts <- artifacts[!vapply(artifacts, is_placeholder_artifact, logical(1))]

    if (!length(artifacts)) next

    for (artifact in artifacts) {
      window_start <- max(1, i - 2)
      window_end <- min(length(lines), i + 2)
      context_window <- paste(lines[window_start:window_end], collapse = " ")
      context_lower <- tolower(context_window)
      lhs <- trimws(sub("<-.*$", "", lines[[i]]))

      is_generated <- grepl(
        "ggsave\\(|write_csv\\(|write\\.csv\\(|save\\(|png\\(|pdf\\(|jpeg\\(|jpg\\(|tiff\\(",
        context_lower
      ) ||
        grepl(
          "fig_path|paper_figures_path|graphs_path|out_dir|out_base",
          context_lower
        ) ||
        grepl(
          "^(out|output|fig|figure|summary|params|high|sum|rank|vpd_quantile).*",
          lhs,
          ignore.case = TRUE
        )

      is_input <- grepl(
        "read_csv\\(|read\\.csv\\(|fread\\(|load\\(|list\\.files\\(",
        context_lower
      )

      role <- if (is_generated) {
        "generated"
      } else if (is_input) {
        "input/reference"
      } else {
        "mentioned"
      }

      rows[[length(rows) + 1]] <- data.frame(
        artifact = artifact,
        role = role,
        script = file.path("code", basename(script)),
        line = i,
        context = trimws(lines[[i]]),
        stringsAsFactors = FALSE
      )
    }
  }
}

index <- if (length(rows)) do.call(rbind, rows) else data.frame(
  artifact = character(),
  script = character(),
  line = integer(),
  context = character()
)

index <- index[order(tolower(index$artifact), index$script, index$line), ]

generated_keys <- paste(
  index$artifact[index$role == "generated"],
  index$script[index$role == "generated"],
  sep = "\r"
)
is_redundant_mention <- index$role == "mentioned" &
  paste(index$artifact, index$script, sep = "\r") %in% generated_keys
index <- index[!is_redundant_mention, ]

if (nzchar(query)) {
  query_lower <- tolower(query)
  keep <- grepl(query_lower, tolower(index$artifact), fixed = TRUE) |
    grepl(query_lower, tolower(index$script), fixed = TRUE) |
    grepl(query_lower, tolower(index$context), fixed = TRUE)
  index <- index[keep, ]
} else {
  index <- index[index$role == "generated", ]
}

if (!nrow(index)) {
  cat("No generated outputs matched")
  if (nzchar(query)) cat(" query: ", query, sep = "")
  cat("\n")
  quit(status = 0)
}

cat(sprintf("%-72s  %-15s  %-86s  %s\n", "ARTIFACT", "ROLE", "SCRIPT", "LINE"))
cat(sprintf("%-72s  %-15s  %-86s  %s\n", paste(rep("-", 72), collapse = ""), paste(rep("-", 15), collapse = ""), paste(rep("-", 86), collapse = ""), "----"))

for (i in seq_len(nrow(index))) {
  cat(sprintf(
    "%-72s  %-15s  %-86s  %s\n",
    index$artifact[[i]],
    index$role[[i]],
    index$script[[i]],
    index$line[[i]]
  ))
}
