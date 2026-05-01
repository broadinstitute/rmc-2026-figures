# Helper functions for RMC and MPC analyses

library(RColorBrewer)
library(GenomicRanges)
library(ggtext)
library(tidyverse)
theme_set(theme_classic())

source("helpers.R")

# Top-level directory paths
rmc_dir = get_script_dir()
data_dir = paste0(rmc_dir, "data/")
analysis_dir = paste0(rmc_dir, "analysis/")
general_analysis_dir = paste0(analysis_dir, "general_results/")

# # Transcript to gene name
transcript_gene_coords_data = read_table(
  paste0(data_dir,"genes_to_trascript.txt")
)

