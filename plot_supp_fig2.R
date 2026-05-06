#!/bin/r env

source(file.path(dirname(normalizePath(
  c(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)),
    tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL),
    tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL))[1]
)), "utils.R"), chdir = TRUE)

source("process_data.R")

