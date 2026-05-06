#!/bin/r env

source(file.path(dirname(normalizePath(
  c(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)),
    tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL),
    tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL))[1]
)), "utils.R"), chdir = TRUE)

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

#########################
### Add annotations to transcripts
#########################
transcript_mc_annots = transcript_n_sections[
  transcript_n_sections$transcript %in% qc_pass_transcripts,
]

# Add transcript coding length
transcript_cds_lengths = as.data.frame(annot_38_canonical_cds) %>%
  group_by(transcript_id) %>%
  summarize(sum(width))
colnames(transcript_cds_lengths) = c("transcript_id", "cds_length")
transcript_mc_annots = left_join(
  transcript_mc_annots,
  transcript_cds_lengths,
  by = c("transcript" = "transcript_id")
)
transcript_mc_annots$log10_cds_length = log10(transcript_mc_annots$cds_length)
transcript_mc_annots$cds_length_pct = ceiling(
  100 * rank(transcript_mc_annots$cds_length) / nrow(transcript_mc_annots)
)

# Add gene name
transcript_mc_annots = left_join(
  transcript_mc_annots,
  constraint_data[, c("transcript", "gene")]
)