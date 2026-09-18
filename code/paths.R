find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "code")) &&
        dir.exists(file.path(current, "data"))
    ) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find repo root containing both code/ and data/ from: ", start)
    }
    current <- parent
  }
}

repo_root <- find_repo_root()
code_path <- file.path(repo_root, "code")

input_path <- file.path(repo_root, "data")
input_folder <- input_path
data_dir <- input_path
output_path <- input_path

fig_path <- file.path(repo_root, "figures")
paper_figures_path <- fig_path
graphs_path <- file.path(repo_root, "graphs")
out_base <- file.path(repo_root, "output")

dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
dir.create(graphs_path, showWarnings = FALSE, recursive = TRUE)
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

# u* filtering helper (thresholds + apply_ustar_filter()); shared by all scripts
source(file.path(code_path, "ustar_filter.R"))

# shared figure palette: season/ENSO colours, shapes, linetypes (see palette.R)
source(file.path(code_path, "palette.R"))

# suppress the stray Rplots.pdf that non-interactive Rscript runs create from
# diagnostic print()/hist() calls; ggsave()/png() devices and interactive
# sessions are unaffected.
if (!interactive()) grDevices::pdf(NULL)
