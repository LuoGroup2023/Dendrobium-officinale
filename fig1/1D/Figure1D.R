# Redraw the four haplotype tracks from MCScanX block and gene-coordinate files.
# Usage: Rscript Figure1D.R GENE_COORDINATES COLLINEARITY [OUTPUT_DIR]
library(grid)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: Rscript Figure1D.R GENE_COORDINATES COLLINEARITY [OUTPUT_DIR]")
}
output_dir <- if (length(args) == 3L) args[3] else "."
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

genes <- read.delim(
  args[1], header = FALSE, quote = "", comment.char = "",
  col.names = c("chr", "gene", "start", "end"),
  colClasses = c("character", "character", "integer", "integer")
)
stopifnot(!anyNA(genes), !anyDuplicated(genes$gene),
          all(genes$start >= 1L), all(genes$start <= genes$end))

haps <- c("A", "B", "C", "D")
tracks <- expand.grid(number = 1:19, hap = haps, stringsAsFactors = FALSE)
tracks$chr <- paste0("chr", tracks$number, tracks$hap)
stopifnot(setequal(genes$chr, tracks$chr))
# The original tracks use the last annotated gene end as their length.
length_bp <- tapply(genes$end, genes$chr, max)
tracks$length_bp <- unname(length_bp[tracks$chr])

mcscan <- readLines(args[2], warn = FALSE)
headers <- grep("^## Alignment ", mcscan)
header_fields <- strcapture(
  "^## Alignment ([0-9]+):.* N=([0-9]+) (chr[0-9]+[ABCD])&(chr[0-9]+[ABCD]) ([a-z]+)$",
  mcscan[headers],
  proto = list(id = integer(), pairs = integer(), a = character(),
               b = character(), orientation = character())
)
stopifnot(length(headers) > 0L, !anyNA(header_fields),
          all(header_fields$orientation %in% c("plus", "minus")),
          !anyDuplicated(header_fields$id), all(header_fields$pairs > 0L))
block_lines <- diff(c(headers, length(mcscan) + 1L)) - 1L
stopifnot(identical(block_lines, header_fields$pairs))

# Use each block's outermost gene coordinates on the two chromosomes.
blocks <- header_fields
for (name in c("a_start", "a_end", "b_start", "b_end")) {
  blocks[[name]] <- integer(nrow(blocks))
}
for (i in seq_len(nrow(blocks))) {
  lines <- mcscan[headers[i] + seq_len(block_lines[i])]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  stopifnot(all(lengths(fields) == 4L))
  pair_fields <- matrix(unlist(fields, use.names = FALSE), ncol = 4L, byrow = TRUE)
  a <- match(pair_fields[, 2], genes$gene)
  b <- match(pair_fields[, 3], genes$gene)
  stopifnot(!anyNA(a), !anyNA(b),
            all(genes$chr[a] == blocks$a[i]), all(genes$chr[b] == blocks$b[i]))
  blocks[i, c("a_start", "a_end", "b_start", "b_end")] <-
    list(min(genes$start[a]), max(genes$end[a]),
         min(genes$start[b]), max(genes$end[b]))
}

a_track <- match(blocks$a, tracks$chr)
b_track <- match(blocks$b, tracks$chr)
a_row <- match(tracks$hap[a_track], haps)
b_row <- match(tracks$hap[b_track], haps)
stopifnot(!anyNA(a_row), !anyNA(b_row))

# Match the original panel: all blocks for A-B, B-C and C-D,
# including blocks linking chromosomes with different numbers.
visible <- which(abs(a_row - b_row) == 1L)
counts <- tabulate(pmin(a_row[visible], b_row[visible]), nbins = 3L)
cat(sprintf("Parsed %d blocks / %d gene pairs. Plotted %d blocks / %d gene pairs.\n",
            nrow(blocks), sum(blocks$pairs), length(visible), sum(blocks$pairs[visible])))
cat(sprintf("Blocks by comparison: A-B %d, B-C %d, C-D %d.\n",
            counts[1], counts[2], counts[3]))

width_mm <- 186
height_mm <- 55
row_y <- seq(16.1, 50.5, length.out = 4)
bar_height <- 1.8
min_gap_mm <- 1.35
left_mm <- 16.5
right_mm <- 181.5
hap_bp <- tapply(tracks$length_bp, tracks$hap, sum)
mm_per_bp <- (right_mm - left_mm - 18 * min_gap_mm) / max(hap_bp)
gap_mm <- (right_mm - left_mm - hap_bp[tracks$hap] * mm_per_bp) / 18
tracks$x <- left_mm + (tracks$number - 1) * gap_mm +
  (ave(tracks$length_bp, tracks$hap, FUN = cumsum) - tracks$length_bp) * mm_per_bp
tracks$width <- tracks$length_bp * mm_per_bp

palettes <- list(
  purple = c(A = "#6595CA", B = "#FCC216", C = "#39C271", D = "#8E69B4")
)

ink <- "#1C2430"
font <- "Arial"
tint <- function(color, strength) {
  rgb <- grDevices::col2rgb(color) / 255
  rgb <- 1 - strength * (1 - rgb)
  grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ])
}

x_mm <- function(value) unit(value, "mm")
y_mm <- function(value) unit(height_mm - value, "mm")

label <- function(text, x, y, size = 8, face = "plain", just = "left", color = ink, angle = 0) {
  grid.text(text, x_mm(x), y_mm(y), just = c(just, "centre"), rot = angle,
            gp = gpar(fontfamily = font, fontsize = size, fontface = face, col = color))
}

# Both Bezier control points lie on the vertical midpoint between tracks.
ribbon <- function(top_x, bottom_x, top_y, bottom_y, fill) {
  t <- seq(0, 1, length.out = 40L)
  blend <- t^2 * (3 - 2 * t)
  y <- top_y + (bottom_y - top_y) * (1.5 * t - 1.5 * t^2 + t^3)
  edge <- function(i) top_x[i] + (bottom_x[i] - top_x[i]) * blend
  grid.polygon(x_mm(c(edge(1), rev(edge(2)))), y_mm(c(y, rev(y))),
               gp = gpar(fill = fill, col = NA))
}

draw_panel <- function(hap_color) {
  track_fill <- setNames(tint(hap_color, 0.2), haps)
  outline_rgb <- grDevices::col2rgb(hap_color) / 255 * 0.66
  track_outline <- setNames(
    grDevices::rgb(outline_rgb[1, ], outline_rgb[2, ], outline_rgb[3, ]),
    haps
  )
  ribbon_fill <- setNames(grDevices::adjustcolor(hap_color, alpha.f = 0.80), haps)
  grid.newpage()
  for (i in visible) {
    upper_is_a <- a_row[i] < b_row[i]
    upper <- if (upper_is_a) a_track[i] else b_track[i]
    lower <- if (upper_is_a) b_track[i] else a_track[i]
    a_bp <- c(blocks$a_start[i], blocks$a_end[i])
    b_bp <- c(blocks$b_start[i], blocks$b_end[i])
    top_bp <- if (upper_is_a) a_bp else b_bp
    bottom_bp <- if (upper_is_a) b_bp else a_bp
    if (blocks$orientation[i] == "minus") bottom_bp <- rev(bottom_bp)
    row <- min(a_row[i], b_row[i])
    ribbon(tracks$x[upper] + top_bp * mm_per_bp,
           tracks$x[lower] + bottom_bp * mm_per_bp,
           row_y[row] + bar_height / 2 - 0.12,
           row_y[row + 1] - bar_height / 2 + 0.12,
           ribbon_fill[haps[row]])
  }

  hap_rows <- match(tracks$hap, haps)
  center <- tracks$x + tracks$width / 2
  for (i in seq_len(nrow(tracks))) {
    hap <- tracks$hap[i]
    y <- row_y[hap_rows[i]]
    if (tracks$number[i] == 1) {
      label(paste0("hap", hap), 1.8, y, size = 8.5, face = "bold")
    }
    grid.roundrect(
      x_mm(center[i]), y_mm(y),
      x_mm(tracks$width[i]), x_mm(bar_height), r = x_mm(0.56),
      gp = gpar(fill = track_fill[hap], col = track_outline[hap], lwd = 0.67)
    )
    if (hap == "A") {
      label(paste0("chr", tracks$number[i]), center[i], y - 2.35,
            size = 8.1, just = "left", angle = 60)
    }
  }

  label("(D)", 1.65, 4.15, size = 10)
  label("Haplotypes", left_mm, 4.15, face = "bold")
  key_x <- c(34.5, 50.1, 65.7, 81.3)
  for (i in seq_along(key_x)) {
    grid.roundrect(x_mm(key_x[i]), y_mm(4.15), x_mm(3.8), x_mm(2.5), r = x_mm(0.38),
                   gp = gpar(fill = track_fill[i], col = track_outline[i], lwd = 0.8))
    label(paste0("hap", haps[i]), key_x[i] + 3.05, 4.15)
  }
  ribbon(c(101.5, 106.6), c(103.2, 108.3), 2.9, 5.4, ribbon_fill["A"])
  label("Collinear block (upper haplotype color)", 109.8, 4.15)
}

for (i in seq_along(palettes)) {
  palette <- palettes[[i]]
  prefix <- file.path(output_dir, "Figure1D_purple_chr_labels_wide")

  grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_mm / 25.4,
                       height = height_mm / 25.4, family = font, bg = "white")
  draw_panel(palette)
  invisible(dev.off())

  svglite::svglite(paste0(prefix, ".svg"), width = width_mm / 25.4,
                   height = height_mm / 25.4, bg = "white")
  draw_panel(palette)
  invisible(dev.off())

  grDevices::png(paste0(prefix, ".png"), width = width_mm, height = height_mm,
                 units = "mm", res = 600, type = "cairo", family = font, bg = "white")
  draw_panel(palette)
  invisible(dev.off())
}




