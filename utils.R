# Helper functions for RMC and MPC analyses

library(RColorBrewer)
library(GenomicRanges)
library(ggtext)
library(tidyverse)
theme_set(theme_classic())

source("helpers.R")

do_chisq_test = function(x1, x2, y1, y2) {
  #' Run 2x2 chi-square association test and print out some useful results.
  #'
  #' @param x1 double. [1,1] in matrix.
  #' @param x2 double. [1,2] in matrix.
  #' @param y1 double. [2,1] in matrix.
  #' @param y2 double. [2,2] in matrix.
  #' @return chisq.test. result.
  
  mat = matrix(c(x1, y1, x2, y2), ncol = 2)
  print(mat)
  print(mat / rowSums(mat))
  chisq = chisq.test(mat)
  print(paste0("Chi-square: ", round(chisq$statistic, 3), ", p=", signif(chisq$p.value, 3)))
  print(paste("Odds ratio:", (mat[1, 1] / mat[2, 1]) / (mat[1, 2] / mat[2, 2])))
  return(chisq)
}

do_fisher_exact_test = function(x1, x2, y1, y2) {
  #' Run 2x2 Fisher exact association test and print out some useful results.
  #'
  #' @param x1 double. [1,1] in matrix.
  #' @param x2 double. [1,2] in matrix.
  #' @param y1 double. [2,1] in matrix.
  #' @param y2 double. [2,2] in matrix.
  #' @return fisher.test. result.
  
  mat = matrix(c(x1, y1, x2, y2), ncol = 2)
  print(mat)
  print(mat / rowSums(mat))
  fisher = fisher.test(mat)
  print(paste("Fisher exact p:", signif(fisher$p.value, 3)))
  print(paste("Odds ratio:", (mat[1, 1] / mat[2, 1]) / (mat[1, 2] / mat[2, 2])))
  return(fisher)
}

# Top-level directory paths
rmc_dir = get_script_dir()
data_dir = paste0(rmc_dir, "/data")
resource_dir = paste0(rmc_dir,"/resources")
analysis_dir = paste0(rmc_dir, "/analysis")
general_analysis_dir = paste0(analysis_dir, "/general_results")

annot_38 <<- rtracklayer::import(paste0(resource_dir,'/gencode.v39.annotation.gtf.gz'))
mcols(annot_38) = mcols(annot_38)[, c("transcript_id", "type", "gene_id", "gene_name")]
mcols(annot_38)$transcript_id = rm_ensembl_version(mcols(annot_38)$transcript_id)
annot_38_canonical_cds = annot_38[
  (annot_38$type == "CDS") &
    (
      annot_38$transcript_id %in% transcript_n_sections[
        transcript_n_sections$transcript %in% qc_pass_transcripts,
        "transcript"
      ]
    )
]

# # Transcript to gene name
transcript_gene_coords_data = read_table(
  paste0(resource_dir,"/genes_to_trascript.txt")
)

# # Transcript-level constraint data from gnomAD v4
constraint_data = read_table(
  paste0(resource_dir,"/gnomad.v4.1.constraint_metrics_max_af_0.001.tsv")
)

# # Transcripts considered QC-pass in gnomAD v4
qc_pass_transcripts = constraint_data$transcript[constraint_data$constraint_flags == "[]"]

# # Autosomal dominant and recessively associated genes from OMIM
omim_data = read_table(paste0(resource_dir,"/omim_table_02_15_2023.tsv"))
omim_data = omim_data[, c("phenotype_description", "gene_symbols", "phenotype_inheritance")]
omim_data = as.data.frame(
  omim_data %>%
    mutate(gene = strsplit(gene_symbols, ",")) %>%
    tidyr::unnest(gene) %>%
    mutate(gene = str_trim(gene)) %>%
    select(-gene_symbols)
)
auto_dominant_transcripts = unique(
  c(
    constraint_data[
      constraint_data$gene %in% omim_data[
        omim_data$phenotype_inheritance == "Autosomal dominant",
        "gene"
      ],
      "transcript"
    ],
    transcript_gene_coords_data[
      transcript_gene_coords_data$gencode_gene %in% omim_data[
        omim_data$phenotype_inheritance == "Autosomal dominant",
        "gene"
      ],
      "transcript"
    ]
  )
)
auto_recessive_transcripts = unique(
  c(
    constraint_data[
      constraint_data$gene %in% omim_data[
        omim_data$phenotype_inheritance == "Autosomal recessive",
        "gene"
      ],
      "transcript"
    ],
    transcript_gene_coords_data[
      transcript_gene_coords_data$gencode_gene %in% omim_data[
        omim_data$phenotype_inheritance == "Autosomal recessive",
        "gene"
      ],
      "transcript"
    ]
  )
)

# Get G2P gene lists
g2p_genes = bind_rows(
  # Iterate through all panels
  lapply(
    list.files(
      path = "~/Desktop/resources/G2P/",
      pattern = "*.csv",
      full.names = TRUE,
    ),
    function(x) {
      df = read.csv(x, stringsAsFactors = FALSE, na.strings = c(NA, ""), colClasses = rep("character",23))
      # Set "No gene mim" to NA value
      df$gene.mim = as.integer(ifelse(df$gene.mim == "No gene mim", NA, df$gene.mim))
      return(df)
    }
  )
)
# LoF consequences
lof_consequences = c(
  "absent gene product",
  "decreased gene product level"
)
# Non-LoF consequences
non_lof_consequences = c(
  "altered gene product structure",
  "increased gene product level"
)
# Genes where mutational consequences include at least one LoF consequence
# and may also include non-LoF
lof_inclusive_g2p_genes = g2p_genes[
  do.call(
    "|",
    lapply(lof_consequences, function(x) grepl(x, g2p_genes$variant.consequence))
  ) &
    (g2p_genes$confidence.category %in% c("definitive", "strong")) &
    is.na(g2p_genes$confidence.value.flag),
]
# Genes where mutational consequences include at least one non-LoF consequence
# and may also include LoF
non_lof_inclusive_g2p_genes = g2p_genes[
  do.call(
    "|",
    lapply(non_lof_consequences, function(x) grepl(x, g2p_genes$variant.consequence))
  ) &
    (g2p_genes$confidence.category %in% c("definitive", "strong")) &
    is.na(g2p_genes$confidence.value.flag),
]
# Genes where mutational consequences are only LoF consequences
lof_exclusive_g2p_genes = lof_inclusive_g2p_genes[
  lof_inclusive_g2p_genes$gene.symbol %in%
    setdiff(lof_inclusive_g2p_genes$gene.symbol, non_lof_inclusive_g2p_genes$gene.symbol),
]
# Genes where mutational consequences are only non-LoF consequences
non_lof_exclusive_g2p_genes = non_lof_inclusive_g2p_genes[
  non_lof_inclusive_g2p_genes$gene.symbol %in%
    setdiff(non_lof_inclusive_g2p_genes$gene.symbol, lof_inclusive_g2p_genes$gene.symbol),
]

# Set RMC freeze values
rmc_freezes = c(2)
current_rmc_freeze = 2

# Set relevant directory paths
get_freeze_analysis_dir_path = function(freeze) {
  paste0(analysis_dir, "freeze", freeze, "_results/")
}
freeze_analysis_dirs = lapply(rmc_freezes, get_freeze_analysis_dir_path)
names(freeze_analysis_dirs) = as.character(rmc_freezes)
current_freeze_analysis_dir = freeze_analysis_dirs[[as.character(current_rmc_freeze)]]
freeze_comparison_dir = paste0(analysis_dir, "freeze_comparison/")
# A variable is assigned for ExAC analysis directory but not currently used
exac_analysis_dir = paste0(analysis_dir, "exac_results/")
paper_dir = paste0(rmc_dir, "paper/raw_figs/")

# Create labels for sources of missense constraint annotation
annot_sources_by_pop_ref = list(
  "gnomad" = paste0("freeze", rmc_freezes)
)
annot_sources = do.call(c, sapply(
  seq_along(annot_sources_by_pop_ref), function(i) {
    if (length(annot_sources_by_pop_ref[[i]]) == 0) {
      names(annot_sources_by_pop_ref)[i]
    } else {
      paste(names(annot_sources_by_pop_ref)[i], annot_sources_by_pop_ref[[i]], sep = "_")
    }
  },
  simplify = FALSE
))
current_annot_source = paste0("gnomad_freeze", current_rmc_freeze)

# Create labels for datasets analyzed
# "intersect" indicates that only transcripts analyzed with all annotation sources will be retained
# "nonintersect" indicates that all transcripts analyzed across any annotation source will be retained
# transcript_intersect_types = c("intersect", "nonintersect")
transcript_intersect_types = c("nonintersect")

# Number of transcripts analyzed for RMC/MPC
n_transcripts = 17841

##############
## Missense RMC utilities
##############
# Read in RMC results with missense O/E by region
mc_sections = lapply(annot_sources, function(x) {
  if (x == "exac") {
    tbl_path = paste0(data_dir, x, "_all_mc.tsv")
  } else {
    gnomad_mc_suffix = str_split(x, "gnomad_")[[1]][2]
    tbl_path = paste0(data_dir, "gnomad_all_mc_", gnomad_mc_suffix, ".tsv")
  }
  tbl = read_table(tbl_path)
  # Recompute section O/E as the file has it capped at 1
  tbl$section_oe = tbl$section_obs / tbl$section_exp
  # Remove any NAs - seems these appear when there is no gnomAD constraint data available
  tbl = tbl[!(is.na(tbl$chr)), ]
  # Filter to QC-pass transcripts
  # tbl = tbl[tbl$transcript %in% qc_pass_transcripts, ]
  return(tbl)
})
names(mc_sections) = annot_sources
mc_section_grs = lapply(mc_sections, function(x) {
  x$chr = x$chr
  x$start = pmin(x$section_start, x$section_end)
  x$end = pmax(x$section_start, x$section_end)
  return(makeGRangesFromDataFrame(x, keep.extra.columns = TRUE))
})

# Read in gnomAD data
mc_tbls = lapply(annot_sources, function(x) {
  tbl = mc_sections[[x]]
  return(tbl)
})
names(mc_tbls) = annot_sources

# Get number of sections per transcript including with no RMC
transcript_n_sections_lst = lapply(mc_tbls, function(x) {
  x %>%
    group_by(transcript) %>%
    summarize(n_sections = n())
})

# Join across annotation sources into one table and write out
transcript_n_sections = transcript_n_sections_lst %>% purrr::reduce(full_join, by = "transcript")
colnames(transcript_n_sections) =
  c("transcript", paste0("n_sections_", names(transcript_n_sections_lst)))

transcript_n_sections = as.data.frame(transcript_n_sections[order(transcript_n_sections$transcript), ])

##############
## Missense O/E utilities
##############
# Define missense O/E column names for annotation
mis_oe_col = "oe"

# Define increments on missense O/E bins for annotation
oe_bin_incrs = c(0.1, 0.2)

bin_incr_to_str = function(bin_incr) {
  #' Convert a numeric bin increment to a string representation.
  #'
  #' @param bin_incr double. Bin increment to convert.
  #' @return string. String representation of the bin increment.
  #' @example bin_incr_to_str(0.1) = "0.1".
  as.character(bin_incr)
}

bin_incr_to_int_str = function(bin_incr) {
  #' Convert a numeric bin increment to a string representation of its integer equivalent.
  #'
  #' @param bin_incr double. Bin increment to convert.
  #' @return string. String representation of the bin increment.
  #' @example bin_incr_to_int_str(0.1) = "1"
  #' @example bin_incr_to_int_str(0.2) = "2"
  as.character(as.integer(bin_incr * 10))
}

# Define O/E bin intervals
oe_bin_intervals = lapply(oe_bin_incrs, function(x) {
  seq(0, 1, x)
})
names(oe_bin_intervals) = sapply(oe_bin_incrs, bin_incr_to_str)

# Define O/E bin labels
oe_bin_labels = lapply(oe_bin_intervals, function(x) {
  bin_intervals = as.character(x)
  bin_intervals[1] = "0"
  bin_intervals[length(bin_intervals)] = "1.0+"
  return(
    sapply(1:(length(bin_intervals) - 1), function(i) {
      paste0(bin_intervals[i], "-", bin_intervals[i + 1])
    })
  )
})

# Define O/E bin column names to annotate
oe_bin_cols = sapply(oe_bin_incrs, function(x) {
  paste0("oe_bin_", str_split(as.character(x), "\\.")[[1]][2])
})
names(oe_bin_cols) = sapply(oe_bin_incrs, bin_incr_to_str)

add_oe_bin_annot = function(oe_tbl, oe_bin_incr, mis_oe_col, oe_bin_col) {
  #' Annotate a table containing O/E values with the bins those O/E values fall in
  #' based on a specific O/E column name and bin increment.
  #'
  #' @param oe_tbl dataframe. Table containing O/E values.
  #' @param oe_bin_incr double. Increment on O/E bins to use.
  #' @param mis_oe_col string. Name of table column containing O/E values.
  #' @param oe_bin_col string. Name to be used for table column containing output O/E bins.
  #' @return dataframe. `oe_tbl` annotated with a column named `oe_bin_col` containing
  #' bins corresponding to O/E values in `mis_oe_col`.
  
  oe_bin_incr_str = bin_incr_to_str(oe_bin_incr)
  
  # Get O/E bins to annotate with
  bins = oe_bin_intervals[[oe_bin_incr_str]]
  
  # Get indices of O/E bin intervals
  oe_tbl[, oe_bin_col] = sapply(oe_tbl[, mis_oe_col], function(x) {
    if (is.na(x)) {
      return(NA)
    } else if (x <= bins[2]) {
      return(1)
    } else if (x > bins[length(bins) - 1]) {
      return(length(bins) - 1)
    } else {
      return(max(which(x > bins[1:(length(bins) - 2)])))
    }
  })
  
  # Convert indices to O/E bin labels
  oe_tbl[, oe_bin_col] = ifelse(
    is.na(oe_tbl[, oe_bin_col]),
    NA,
    oe_bin_labels[[oe_bin_incr_str]][oe_tbl[, oe_bin_col]]
  )
  return(oe_tbl)
}

# Define plotting colors for O/E bins using Reds palette
bin_colors = lapply(oe_bin_incrs, function(x) {
  if (x == 0.1) {
    rev(colorRampPalette(brewer.pal(name = "Reds", n = 8))(10))
  } else if (x == 0.2) {
    rev(colorRampPalette(brewer.pal(name = "Reds", n = 5))(5))
  }
})
names(bin_colors) = sapply(oe_bin_incrs, bin_incr_to_str)

##############
## MPC utilities
##############
# Define MPC column names to annotate with
mpc_cols = lapply(annot_sources, function(x) paste0(x, "_mpc"))
names(mpc_cols) = annot_sources

trim_sort_mpc_thresholds_for_bins = function(thresholds) {
  thresholds = sort(unique(thresholds))
  # Assert all input values are nonnegative
  stopifnot(all(thresholds >= 0))
  # # Remove 0 if included (first bin will be < smallest nonzero value)
  # if (thresholds[1] == 0) {
  #   thresholds = thresholds[-1]
  # }
  return(thresholds)
}

mpc_thresholds_to_bins = function(thresholds) {
  thresholds = trim_sort_mpc_thresholds_for_bins(thresholds)
  # Create bins
  bins = c(
    paste0("<", thresholds[1]),
    unlist(
      Map(
        function(x, y) paste(x, y, sep = "<=MPC<"),
        thresholds[-length(thresholds)],
        thresholds[-1]
      )
    ),
    paste0(">=", thresholds[length(thresholds)])
  )
  return(bins)
}

add_mpc_bin_annot = function(mpc_tbl, mpc_thresholds = seq(1, 2), mpc_col, mpc_bin_col) {
  #' Annotate a table containing MPC values with the bins those MPC values fall in
  #' based on specified bin thresholds.
  #'
  #' @param mpc_tbl dataframe. Table containing MPC values.
  #' @param mpc_thresholds vector[double]. Thresholds separating MPC bins. Must be nonnegative.
  #' @param mpc_col string. Name of table column containing MPC values.
  #' @param mpc_bin_col string. Name to be used for table column containing output MPC bins.
  #' @return dataframe. `mpc_tbl` annotated with a column named `mpc_bin_col` containing
  #' bins corresponding to MPC values in `mpc_col`.
  
  mpc_thresholds = trim_sort_mpc_thresholds_for_bins(mpc_thresholds)
  bins = mpc_thresholds_to_bins(mpc_thresholds)
  
  # Get indices of MPC bin intervals
  mpc_tbl[, mpc_bin_col] = sapply(mpc_tbl[, mpc_col], function(x) {
    if (is.na(x)) {
      return(NA)
    } else if (x < mpc_thresholds[1]) {
      return(1)
    } else if (x >= mpc_thresholds[length(mpc_thresholds)]) {
      return(length(mpc_thresholds) + 1)
    } else {
      return(max(which(x >= mpc_thresholds)) + 1)
    }
  })
  
  # Convert indices to MPC bin labels
  mpc_tbl[, mpc_bin_col] = ifelse(
    is.na(mpc_tbl[, mpc_bin_col]),
    NA,
    bins[mpc_tbl[, mpc_bin_col]]
  )
  return(mpc_tbl)
}

##############
## Plot utilities
##############
# Plot dimensions
consensus_plot_height = 5
consensus_plot_width = 6

# Extension type to save plot as
plot_extension = "svg"

# Text characteristics
axis_title_size = 15
axis_text_size = 12
legend_title_size = 12
legend_text_size = 11
other_text_size = 11
consensus_text_color = "#000000"

# Function to automate above
ggsave_plot = function(filename,
                       p,
                       extension = NA,
                       height = NA,
                       width = NA,
                       x_title_size = NA,
                       y_title_size = NA,
                       x_text_size = NA,
                       y_text_size = NA,
                       l_title_size = NA,
                       l_text_size = NA,
                       text_color = NA,
                       ...) {
  #' Export plots to files with ggplot using consensus parameters unless specified.
  #' NOTE: Axis and legend text are in form `element_markdown` not `element_text`.
  #'
  #' @param filename string. File path without extension. E.g. "myplot".
  #' @param p ggplot. Plot to save.
  #' @param extension string. File extension.
  #' @param height double. Plot height in inches.
  #' @param width double. Plot width in inches.
  #' @param x_title_size double. Size of x-axis title font.
  #' @param y_title_size double. Size of y-axis title font.
  #' @param x_text_size double. Size of x-axis text font.
  #' @param y_text_size double. Size of y-axis text font.
  #' @param l_title_size double. Size of legend title font.
  #' @param l_text_size double. Size of legend text font.
  #' @param text_color string. Color for all text in plot.
  #' @param ... Remaining arguments to pass to ggsave.
  
  x_title = ifelse(is.na(x_title_size), axis_title_size, x_title_size)
  y_title = ifelse(is.na(y_title_size), axis_title_size, y_title_size)
  x_text = ifelse(is.na(x_text_size), axis_text_size, x_text_size)
  y_text = ifelse(is.na(y_text_size), axis_text_size, y_text_size)
  l_title = ifelse(is.na(l_title_size), legend_title_size, l_title_size)
  l_text = ifelse(is.na(l_text_size), legend_text_size, l_text_size)
  c = ifelse(is.na(text_color), consensus_text_color, text_color)
  
  p = p +
    theme(
      text = element_text(color = c, size = other_text_size),
      axis.title.x = element_markdown(size = x_title, color = c),
      axis.title.y = element_markdown(size = y_title, color = c),
      axis.text.x = element_markdown(size = x_text, color = c),
      axis.text.y = element_markdown(size = y_text, color = c),
      legend.title = element_markdown(size = l_title, color = c),
      legend.text = element_markdown(size = l_text, color = c)
    )
  
  e = ifelse(is.na(extension), plot_extension, extension)
  ggsave(
    paste0(filename, ".", e),
    p,
    height = ifelse(is.na(height), consensus_plot_height, height),
    width = ifelse(is.na(width), consensus_plot_width, width),
    units = "in",
    ...
  )
  # Save again as png if not already
  if (e != "png") {
    ggsave(
      paste0(filename, ".png"),
      p,
      height = ifelse(is.na(height), consensus_plot_height, height),
      width = ifelse(is.na(width), consensus_plot_width, width),
      units = "in",
      ...
    )
  }
}

# Color schemes
plp_blb_colors = list("P/LP" = "#9D1309", "B/LB" = "#87A4DC", "NA" = "#D1D1D1")