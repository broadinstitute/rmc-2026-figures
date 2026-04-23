library(ggplot2)
library(viridis)
library(ggpointdensity)
library(ggExtra)

prop_uncons_rmc_transcript = read_table(paste0(data_dir, "prop_uncons_per_rmc_transcript.tsv.bgz"))
prop_uncons_no_rmc_transcript = read_table(paste0(data_dir, "prop_uncons_per_transcript_no_rmc.tsv.bgz"))
prop_uncons_section = read_table(paste0(data_dir, "prop_uncons_per_rmc_section.tsv.bgz"))

# Panel a) All transcripts, OE/phyloP at transcript level
# Spearman rho : 0.4982611; p < 1e-50
all_transcript <- ggplot(rbind(prop_uncons_rmc_transcript, prop_uncons_no_rmc_transcript), aes(x = oe, y = prop_uncons)) +
    geom_pointdensity(show.legend = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_color_viridis() +
    xlim(c(0, 3)) +
    ylab("Proportion bases unconserved") +
    xlab("Missense OE") +
    # theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.title = element_text(size = 12), legend.text = element_text(size = 11)) +
    theme_classic() +
    guides(color = guide_colorbar(title = "Density")) +
    theme(
        text = element_text(color = "black", size = other_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        axis.text.x = element_text(size = axis_text_size),
        axis.text.y = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size)
    )
pa = ggMarginal(all_transcript, type = "density")
for (ff in c("png", "svg")) {
    ggsave(
        paste0(
            current_freeze_analysis_dir,
            "rmc/constraint/rmc_vs_conservation_all_transcript.",
            ff
        ),
        pa,
        width = consensus_plot_width,
        height = consensus_plot_height
    )
}


# Panel b) All transcripts, OE/phyloP at MCR-level
# Spearman rho : 0.4389949; p < 1e-50
all_transcript_mcr <- ggplot(rbind(prop_uncons_section[, c(5, 2:4,6)], prop_uncons_no_rmc_transcript), aes(x = oe, y = prop_uncons)) +
    geom_pointdensity(show.legend = TRUE) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_color_viridis() +
    xlim(c(0, 3)) +
    ylab("Proportion bases unconserved") +
    xlab("Missense OE") +
    # theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.title = element_text(size = 12), legend.text = element_text(size = 11)) +
    theme_classic() +
    guides(color = guide_colorbar(title = "Density")) +
    theme(
        text = element_text(color = "black", size = other_text_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        axis.text.x = element_text(size = axis_text_size),
        axis.text.y = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        #legend.text = element_text(size = legend_text_size)
        legend.text = element_blank(),
        legend.ticks = element_blank()
    )
pb = ggMarginal(all_transcript_mcr, type = "density")
for (ff in c("png", "svg")) {
    ggsave(
        paste0(
            current_freeze_analysis_dir,
            "rmc/constraint/rmc_vs_conservation_all_mcr.",
            ff
        ),
        pb,
        width = consensus_plot_width,
        height = consensus_plot_height
    )
}


# Panel c) 2+MCR transcripts, OE/phyloP at transcript level
# Spearman rho : 0.4504078 ; p < 1e-50
rmc_transcript <- ggplot(prop_uncons_rmc_transcript, aes(x = oe, y = prop_uncons)) +
  geom_pointdensity(show.legend = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_viridis() +
  xlim(c(0, 3)) +
  ylab("Proportion bases unconserved") +
  xlab("Missense OE") +
  # theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.title = element_text(size = 12), legend.text = element_text(size = 11)) +
  theme_classic() +
  guides(color = guide_colorbar(title = "Density")) +
  theme(
    text = element_text(color = "black", size = other_text_size),
    axis.title.x = element_text(size = axis_title_size),
    axis.title.y = element_text(size = axis_title_size),
    axis.text.x = element_text(size = axis_text_size),
    axis.text.y = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size)
  )
pc = ggMarginal(rmc_transcript, type = "density")
for (ff in c("png", "svg")) {
  ggsave(
    paste0(
      current_freeze_analysis_dir,
      "rmc/constraint/rmc_vs_conservation_rmc_transcript.",
      ff
    ),
    pc,
    width = consensus_plot_width,
    height = consensus_plot_height
  )
}


# Panel d) 2+MCR transcripts, OE/phyloP at MCR-level
# Spearman rho : 0.3224149; p < 1e-50
rmc_mcr <- ggplot(prop_uncons_section, aes(x = oe, y = prop_uncons)) +
  geom_pointdensity(show.legend = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_viridis() +
  xlim(c(0, 3)) +
  ylab("Proportion bases unconserved") +
  xlab("Missense OE") +
  # theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.title = element_text(size = 12), legend.text = element_text(size = 11)) +
  theme_classic() +
  guides(color = guide_colorbar(title = "Density")) +
  theme(
    text = element_text(color = "black", size = other_text_size),
    axis.title.x = element_text(size = axis_title_size),
    axis.title.y = element_text(size = axis_title_size),
    axis.text.x = element_text(size = axis_text_size),
    axis.text.y = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    #legend.text = element_text(size = legend_text_size)
    legend.text = element_blank(),
    legend.ticks = element_blank()
  )
pd = ggMarginal(rmc_mcr, type = "density")
for (ff in c("png", "svg")) {
  ggsave(
    paste0(
      current_freeze_analysis_dir,
      "rmc/constraint/rmc_vs_conservation_rmc_section.",
      ff
    ),
    pd,
    width = consensus_plot_width,
    height = consensus_plot_height
  )
}
