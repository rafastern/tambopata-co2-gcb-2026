# =============================================================================
# palette.R  -  Shared color / shape / linetype conventions for ALL figures.
# Sourced by paths.R so every script uses the same encoding.
#
# Design principle: SEASON and ENSO use DISJOINT color sets, so a given colour
# always means the same thing across every figure and on any channel:
#   - Season  -> BLUE (wet) / RED (dry)         [+ point SHAPE as redundant cue]
#   - ENSO    -> GOLD (El Nino) / PURPLE (La Nina) / GREY (Neutral)   [lines]
#                use `enso_fill` for filled background bands, not `enso_cols`
#                                                [+ LINETYPE as redundant cue]
# ENSO deliberately avoids blue and red so it never collides with season.
# Colours are chosen to be colour-vision-deficiency friendly (Okabe-Ito-based);
# there is no red/green pairing anywhere.
#
# Keys include every label spelling used across scripts (ASCII + accented,
# lower/upper case) so the same vector is a drop-in everywhere.
# =============================================================================

## ---- SEASON -----------------------------------------------------------------
season_cols <- c(
  "wet" = "blue3", "Wet" = "blue3",
  "dry" = "red3",  "Dry" = "red3"
)
# Redundant, colour-independent cue: circle (wet) vs triangle (dry)
season_shapes <- c(
  "wet" = 16, "Wet" = 16,
  "dry" = 17, "Dry" = 17
)
# Filled-marker shapes (for scale_shape_manual with fill aesthetics)
season_shapes_filled <- c(
  "wet" = 21, "Wet" = 21,
  "dry" = 24, "Dry" = 24
)
# Light season tints for filled areas (boxplots, background bands). Same hues as
# season_cols but pale, so black median/whisker lines stay readable inside a box.
season_box_fills <- c(
  "wet" = "#9DC3E6", "Wet" = "#9DC3E6",
  "dry" = "#F4B6B0", "Dry" = "#F4B6B0"
)

## ---- ENSO -------------------------------------------------------------------
# LINE / POINT colours. Use `enso_fill` below for filled background rectangles;
# these are deliberately dark and would swamp a background band.
#
# Do not lighten these. The previous values (goldenrod2 #EEB422 and #CC79A7)
# were unreadable as lines: goldenrod2 had a WCAG contrast ratio of 1.88 against
# plain white, before any background, and 1.32-1.36 over the season bands, where
# 3.0 is the minimum for non-text graphical objects. A co-author flagged it on
# the environmental stack figure.
#
# Contrast ratios of the values below (season bands drawn at alpha 0.25):
#                 white   wet band   dry band
#   El Nino        3.54     3.07       3.11
#   La Nina        5.87     5.09       5.15
#   Neutral        5.74     4.98       5.04
#
# La Nina had to move bluer, not merely darker: a darkened version of the old
# pink (#A8497F) reached contrast 5.35 but collapsed onto grey40 under red-blind
# simulation (separation 23, against 77 for the old pink), because magenta loses
# its red channel. Minimum pairwise separation under Vienot-Brettel-Mollon
# deuter-, prot- and tritanopia simulation is 72 for the triple below, against
# 66 for the old one.
enso_cols <- c(
  "El Nino" = "#C77400", "El Niño" = "#C77400",
  "La Nina" = "#8E44AD", "La Niña" = "#8E44AD",
  "Neutral" = "grey40",  "neutral" = "grey40"
)
# Lighter tints for background shading (geom_rect fill); pair with alpha ~0.9.
# Same hues as enso_cols but pale, so lines and points drawn on top stay legible.
#
# Do not lighten the La Nina tint or darken the Neutral one. The previous pair
# (#E4D9EE and #E6E6E6) were separated by a colour distance of only 13.9 in
# normal vision and 14.3 under deuteranopia, i.e. the same colour to everyone;
# a co-author flagged Fig. 3 as unreadable. The values below give a minimum
# pairwise separation of 48.5 across normal vision and deuter-, prot- and
# tritanopia, while keeping contrast >= 3.7 against the season point colours.
#
# Colour alone is still not enough for these pale fills: the three differ by only
# ~9 percentage points of greyscale luminance, so a figure whose ONLY ENSO cue is
# the band colour is unreadable in black and white. Fig. 3 therefore also carries
# a labelled ENSO strip (code/04). Add the same to any new ENSO-shaded figure.
enso_fill <- c(
  "El Nino" = "#F6E2B3", "El Niño" = "#F6E2B3",
  "La Nina" = "#D5C2E8", "La Niña" = "#D5C2E8",
  "Neutral" = "#F0F0F0", "neutral" = "#F0F0F0"
)
# Redundant, colour-independent cue for ENSO lines
enso_linetypes <- c(
  "El Nino" = "solid",  "El Niño" = "solid",
  "La Nina" = "dashed", "La Niña" = "dashed",
  "Neutral" = "dotted", "neutral" = "dotted"
)
