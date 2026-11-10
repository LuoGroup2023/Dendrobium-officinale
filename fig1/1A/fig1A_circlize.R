#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", argv, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
setwd(script_dir)
project_lib <- file.path(script_dir, ".Rlib")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
suppressPackageStartupMessages({
  library(circlize)
  library(grid)
})
# Change these parameters to adjust typography and export size.
font_family <- "Arial"
options(warn = 1)
base_font_pt <- 12
chr_label_cex <- 1.0   # 12 pt; original 0.6 (7.2 pt)
tick_label_cex <- 0.8  # 9.6 pt; original 0.5 (6 pt)
figure_width_in <- 8
figure_height_in <- 8
preview_dpi <- 160
args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args)) args[1] else "circos_plot3_Arial_large"

ref <- read.table("genome.txt", header = TRUE, stringsAsFactors = FALSE)
stopifnot(identical(names(ref), c("chr", "Length")), !anyDuplicated(ref$chr), all(ref$Length > 0))
# All observations are retained; read all density files without a header.
input_files <- c("gene_density_500k.txt", "Copia_coverage.txt", "Gypsy_coverage.txt", "500k_gc.txt")
for (f in input_files) {
  d <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  stopifnot(ncol(d) == 4, all(d[[1]] %in% ref$chr), setequal(d[[1]], ref$chr),
            all(is.finite(d[[4]])), all(d[[2]] >= 0), all(d[[3]] > d[[2]]),
            all(d[[3]] <= ref$Length[match(d[[1]], ref$chr)]))
  message(f, ": ", nrow(d), " windows retained")
}
gap_list <- rep(2, nrow(ref)); names(gap_list) <- ref$chr
gap_list["chr19"] <- 8

draw_circos <- function() {
  par(family = font_family, ps = base_font_pt)
circos.clear()
circos.par("track.height"=0.8,gap.degree= gap_list,start.degree = 86,clock.wise = T,
           cell.padding=c(0,0,0,0),
           track.margin = c(0.01, 0.02))



circos.initialize(factors=ref$chr,          
                  xlim = cbind(rep(0, nrow(ref)), ref$Length))

col_text <- "black" 
bg_color <- "#00ADFF" 


circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    Genome <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    circos.text(
      mean(xlim), 
      mean(ylim),
      Genome,
      cex = chr_label_cex,   
      col = col_text, 
      facing = "bending.inside",
      niceFacing = TRUE
    )
  },
  bg.col = bg_color,  
  bg.border = NA,     
  track.height = 0.06
)



brk <- seq(0, 100, by = 20) * 10^6 


circos.track(
  track.index = get.current.track.index(),
  panel.fun = function(x, y) {
    circos.axis(
      h = "top",
      major.at = brk,
      labels = paste0(brk / 10^6),  
      labels.cex = tick_label_cex,                   
      col = col_text,                    
      labels.col = col_text,              
      lwd = 0.8,
      labels.facing = "clockwise"
    )
  },
  bg.border = F
)


##gene_density
gene_density <- read.table("gene_density_500k.txt", header = FALSE)
colnames(gene_density) <- c("chr","start","end","value")
col_heatmap <- colorRamp2(
  breaks = c(min(gene_density$value), 
             median(gene_density$value),
             max(gene_density$value)),
  colors = c("#FFFFFF", "#eb8c8d", "#da3032")
)

circos.genomicHeatmap(
  gene_density,
  col = col_heatmap,
  side = "inside",
  border = NA,
  line_col = NA,
  heatmap_height = 0.08,
  connection_height = 0.005,
  track.margin = c(0.02, 0.02))
 


# Draw the same continuous legend on the active Cairo device. This avoids
# ComplexHeatmap's temporary PostScript device when measuring Arial text.
legend_at <- pretty(range(gene_density$value), n = 4)
legend_values <- seq(min(legend_at), max(legend_at), length.out = 201)
grid::pushViewport(grid::viewport(x = unit(0.84, "npc"), y = unit(0.947, "npc"),
                                width = unit(70, "pt"), height = unit(80, "pt"),
                                just = c("left", "top")))
grid::grid.text("Gene density", x = unit(0, "pt"), y = unit(80, "pt"),
                just = c("left", "bottom"),
                gp = grid::gpar(fontfamily = font_family, fontsize = 12, fontface = "bold"))
bar_y <- 6; bar_h <- 66; bar_w <- 12
for (i in seq_len(length(legend_values) - 1)) {
  grid::grid.rect(x = unit(0, "pt"),
                  y = unit(bar_y + (i - 1) / 200 * bar_h, "pt"),
                  width = unit(bar_w, "pt"), height = unit(bar_h / 200 + 0.02, "pt"),
                  just = c("left", "bottom"),
                  gp = grid::gpar(fill = col_heatmap(mean(legend_values[c(i, i + 1)])), col = NA))
}
for (v in legend_at) {
  yy <- bar_y + (v - min(legend_at)) / diff(range(legend_at)) * bar_h
  grid::grid.segments(x0 = unit(c(0, bar_w - 3), "pt"), x1 = unit(c(3, bar_w), "pt"),
                      y0 = unit(yy, "pt"), y1 = unit(yy, "pt"), gp = grid::gpar(col = "white"))
  grid::grid.text(as.character(v), x = unit(bar_w + 3, "pt"), y = unit(yy, "pt"),
                  just = "left", gp = grid::gpar(fontfamily = font_family, fontsize = 10))
}
grid::popViewport()


##Copia

Copia <- read.table("Copia_coverage.txt", header = FALSE)
colnames(Copia) <- c("chr","start","end","value")


circos.genomicTrackPlotRegion(Copia, 
                              track.height = 0.08, 
                              bg.border = "grey", 
                              bg.col = NA,
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "h", col = "#1f78b4", lwd = 1, ...)})

## Gypsy 
Gypsy <- read.table("Gypsy_coverage.txt", header = FALSE)
colnames(Gypsy) <- c("chr","start","end","value")


circos.genomicTrackPlotRegion(Gypsy,
                              track.height = 0.08, 
                              bg.border = "#CFCFCF", 
                              bg.col = NA,
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "h", col = "grey", lwd = 1, ...)})

### GC content 
GC <- read.table('500k_gc.txt')
colnames(GC) <- c("chr","start","end","value")

circos.genomicTrackPlotRegion(
  GC,
  track.height = 0.08,
  bg.col = "white",
  bg.border = "grey",  
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, type = "l", col = 'red', ...)
  }
)


circos.clear()


}
# Cairo embeds Arial as selectable text; it does not outline the labels.
if (!capabilities("cairo")) stop("R was built without Cairo support")
grDevices::cairo_pdf(paste0(output_prefix, ".pdf"), width = figure_width_in,
                    height = figure_height_in, family = font_family, pointsize = base_font_pt)
tryCatch(draw_circos(), finally = dev.off())
grDevices::png(paste0(output_prefix, ".png"), width = figure_width_in,
              height = figure_height_in, units = "in", res = preview_dpi,
              type = "cairo", family = font_family, pointsize = base_font_pt)
tryCatch(draw_circos(), finally = dev.off())
writeLines(capture.output(sessionInfo()), paste0(output_prefix, ".sessionInfo.txt"))
message("Saved ", output_prefix, ".pdf and .png")
