args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."), mustWork = TRUE)
} else normalizePath(getwd(), mustWork = TRUE)
input_file <- file.path(repo_root, "analysis", "outputs", "osa_inputs.rds")
output_dir <- file.path(repo_root, "analysis", "outputs", "osa")
if (!file.exists(input_file)) stop("Missing OSA inputs: ", input_file)
suppressPackageStartupMessages({ library(afscOSA); library(dplyr); library(ggplot2); library(cowplot) })
inputs <- readRDS(input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
runs <- lapply(names(inputs), function(fleet) {
  x <- inputs[[fleet]]
  afscOSA::run_osa(obs = x$obs, exp = x$exp, N = x$n_eff, fleet = fleet,
                   index = x$index, years = x$years, index_label = "Age")
})
names(runs) <- names(inputs)
plots <- afscOSA::plot_osa(runs, outpath = output_dir, figheight = 9, figwidth = 12)
summary <- dplyr::bind_rows(lapply(runs, `[[`, "res")) |>
  dplyr::group_by(fleet) |>
  dplyr::summarise(residuals = dplyr::n(), mean = mean(resid, na.rm = TRUE),
    sdnr = stats::sd(resid, na.rm = TRUE),
    lower_2.5_pct = stats::quantile(resid, .025, na.rm = TRUE),
    upper_97.5_pct = stats::quantile(resid, .975, na.rm = TRUE), .groups = "drop")
utils::write.csv(summary, file.path(output_dir, "osa_summary.csv"), row.names = FALSE)
saveRDS(list(runs = runs, plots = plots, summary = summary, created = Sys.time(),
  r_version = R.version.string, afscOSA_version = as.character(utils::packageVersion("afscOSA")),
  input_md5 = unname(tools::md5sum(input_file)), lineage = attr(inputs, "lineage")),
  file.path(output_dir, "sporc_ebswp_osa_residuals.rds"))
message("Wrote SPoRC OSA diagnostics to ", output_dir)
