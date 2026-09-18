# combine the diel NEE cycle, the cumulative diel NEE curves and the
# environmental stack into manuscript Figure 2.
#
# Layout: two columns at the GCB maximum width of 180 mm (7.1 in).
#   left  column: (a) diel NEE over (b) cumulative diel NEE
#   right column: (c) environmental stack, spanning the full height
#
# Each source panel is authored by its producing script at exactly the size of
# the slot it occupies here, so every raster is placed at 1:1 and the type is the
# same size in all three. That is why the panels come from dedicated
# "Fig2<x>_..._panel.png" exports rather than from the standalone figures: the
# standalone diel panel is 12 in wide, and squeezing it into a 3.7 in slot would
# leave its labels at about 4.6 pt.
#
#   (a) code/04_diel_cycles_carbon_balance.R  -> Fig2a_diel_NEE_panel.png
#   (b) code/24_cumulative_diel_nee.R         -> Fig2b_cumulative_diel_NEE_panel.png
#   (c) code/14_env_enso_seasonal_stack.R     -> Fig2c_env_stack_panel.png

suppressPackageStartupMessages({
  library(grid)
  library(png)
})

# paths
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

panel_a_png <- file.path(fig_path, "Fig2a_diel_NEE_panel.png")
panel_b_png <- file.path(fig_path, "Fig2b_cumulative_diel_NEE_panel.png")
panel_c_png <- file.path(fig_path, "Fig2c_env_stack_panel.png")

out_png <- file.path(fig_path, "Fig2_combined_env_stack_plus_NEE_diel.png")

inputs <- c(panel_a_png, panel_b_png, panel_c_png)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("missing input panel(s):\n  - ", paste(missing_inputs, collapse = "\n  - "),
       "\nRun code/04, code/14 and code/24 first.")
}

# ─────────────────────────────────────────────────────────────────────────────
# geometry (inches). Must match the ggsave sizes in the producing scripts.

left_w   <- 3.70
right_w  <- 3.25
col_gap  <- 0.15
row_gap  <- 0.12

panel_a_h <- 3.10
panel_b_h <- 2.70
panel_c_h <- panel_a_h + row_gap + panel_b_h   # 5.92, the right column height

total_w <- left_w + col_gap + right_w          # 7.10 in = 180 mm
total_h <- panel_c_h

tag_size <- 10

img_a <- png::readPNG(panel_a_png)
img_b <- png::readPNG(panel_b_png)
img_c <- png::readPNG(panel_c_png)

# warn if a panel was not authored at the size of its slot: it would still draw,
# but stretched, and its type would no longer match the other panels
check_size <- function(img, w_in, h_in, label) {
  got_w <- dim(img)[2] / 300
  got_h <- dim(img)[1] / 300
  if (abs(got_w - w_in) > 0.02 || abs(got_h - h_in) > 0.02) {
    warning(sprintf(
      "panel %s is %.2f x %.2f in but its slot is %.2f x %.2f in; re-render it at the slot size",
      label, got_w, got_h, w_in, h_in), call. = FALSE)
  }
}
check_size(img_a, left_w,  panel_a_h, "a")
check_size(img_b, left_w,  panel_b_h, "b")
check_size(img_c, right_w, panel_c_h, "c")

draw_labeled_panel <- function(img, label) {
  grid::grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
  grid::grid.text(
    label,
    x = unit(0.02, "in"),
    y = unit(1, "npc") - unit(0.02, "in"),
    just = c("left", "top"),
    gp = grid::gpar(fontface = "bold", fontsize = tag_size, col = "black")
  )
}

draw_combined_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      nrow = 3, ncol = 3,
      widths  = unit(c(left_w, col_gap, right_w), "in"),
      heights = unit(c(panel_a_h, row_gap, panel_b_h), "in")
    )
  ))

  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw_labeled_panel(img_a, "a")
  grid::popViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  draw_labeled_panel(img_b, "b")
  grid::popViewport()

  # the environmental stack spans all three rows of the right-hand column
  grid::pushViewport(grid::viewport(layout.pos.row = 1:3, layout.pos.col = 3))
  draw_labeled_panel(img_c, "c")
  grid::popViewport()

  grid::popViewport()
}

grDevices::png(
  filename = out_png,
  width = total_w,
  height = total_h,
  units = "in",
  res = 300
)
draw_combined_figure()
grDevices::dev.off()

message("saved combined Figure 2:")
message("  - ", out_png)
message("source panels:")
message("  - a: ", panel_a_png)
message("  - b: ", panel_b_png)
message("  - c: ", panel_c_png)
message("canvas: ", total_w, " x ", round(total_h, 2), " in (", round(total_w * 25.4), " mm wide)")
