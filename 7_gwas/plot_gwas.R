###This script is used to plot the significant variants for each trait###

####The plot will only show the significant variants in the chromosomes tha contain them displaying the results of the three datasets (upper from lower): SVs, SNPs and SVs+SNPs (both)

#1. Load the necessary packages

library("data.table")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("scales")

#2. Set the plot colors

colors <- brewer.pal(8, "Dark2")[c(1, 4)]

#3. Load the the files that contain the significant variants detected with each dataset for a given trait

svs <- fread("/path/to/trait_str_var.txt")
snps <- fread("/path/to/trait_snps.txt")
both <- fread("/path/to/trait_both.txt")

#4. Indicate the chromosomes that will appear in the plot

chromosomes_with_data <- setdiff(1:12, c(a,b,c,d,...)) #The a,b,c,d,... are the chromosomes that don't contain significant variants. If all the chromosomes contain significant variants, indicate a number different from any chromosome contained in the genome

#5. Combine to get global chromosome position ranges

combined <- rbindlist(list(svs, snps, both), fill = TRUE)
chr_ranges <- combined[Chromosome %in% chromosomes_with_data, .(
  xmin = min(Position, na.rm = TRUE),
  xmax = max(Position, na.rm = TRUE)
), by = Chromosome]

#6. Modified function to control x-axis visibility

generate_plot_list <- function(df, colors, chr_ranges, show_x_axis = TRUE) {
  plot_list <- list()
  for (i in seq_along(chromosomes_with_data)) {
    chr_num <- chromosomes_with_data[i]
    chr_data <- df[df$Chromosome == chr_num]
    xmin <- chr_ranges[Chromosome == chr_num, xmin] / 1e6  # Convert to Mb
    xmax <- chr_ranges[Chromosome == chr_num, xmax] / 1e6
    
    show_y_axis <- (i == 1)  # Show y-axis only on the first plot of the row
    
    p <- ggplot(chr_data, aes(x = Position / 1e6, y = -log10(Pvalue))) +
      geom_point(
        aes(shape = Type, color = Type), 
        size = 1, 
        stroke = 1
      ) +
      scale_shape_manual(values = c("SNP" = 21, "SV" = 24)) +
      scale_color_manual(values = c("SNP" = colors[1], "SV" = colors[2])) +
      scale_x_continuous(limits = c(xmin, xmax), labels = function(x) paste0(round(x), " Mb")) +
      ylim(0, 45) +
      theme_bw() +
      labs(
        x = paste("Chr", chr_num),
        y = if (show_y_axis) expression(-log[10](italic(P))) else NULL,
        title = NULL
      ) +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.2),
        axis.text.x  = if (show_x_axis) element_text(angle = 90, vjust = 0.5, size = 8) else element_blank(),
        axis.title.x = if (show_x_axis) element_text(size = 8) else element_blank(),
        axis.text.y  = if (show_y_axis) element_text(size = 8) else element_blank(),
        axis.title.y = if (show_y_axis) element_text(size = 10) else element_blank(),
        axis.ticks.y = if (show_y_axis) element_line() else element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
      )
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  return(plot_list)
}

#7. Generate plots

svs_plots  <- generate_plot_list(svs,  colors, chr_ranges, show_x_axis = FALSE)
snps_plots <- generate_plot_list(snps, colors, chr_ranges, show_x_axis = FALSE)
both_plots <- generate_plot_list(both, colors, chr_ranges, show_x_axis = TRUE)

#8. Combine all plots in three rows

final_plot <- wrap_plots(
  wrap_plots(plotlist = svs_plots, nrow = 1),
  wrap_plots(plotlist = snps_plots, nrow = 1),
  wrap_plots(plotlist = both_plots, nrow = 1),
  ncol = 1
)

print(final_plot)

#9. Clean the environment

rm(list = ls())
gc()
