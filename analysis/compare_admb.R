suppressPackageStartupMessages({
  library(yaml)
})

source("R/extract_results.R")

cfg <- yaml::read_yaml("config/base.yml")
model_path <- file.path(cfg$paths$outputs_dir, "base.rds")

if (!file.exists(model_path)) {
  stop("Missing SPoRC base fit: ", model_path, "\nRun Rscript analysis/run_scenarios.R first.")
}
if (!file.exists(cfg$paths$admb_results_rds)) {
  stop("Missing ADMB comparison file: ", cfg$paths$admb_results_rds)
}

data <- readRDS(cfg$paths$data_rds)
sporc <- readRDS(model_path)
admb <- readRDS(cfg$paths$admb_results_rds)

sporc_ts <- extract_sporc_timeseries(sporc, years = data$years)
comparison <- compare_timeseries(sporc_ts, admb$timeseries)

dir.create(cfg$paths$comparison_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(
    meta = list(
      sporc_model = model_path,
      admb_model = cfg$paths$admb_results_rds,
      created = as.character(Sys.time())
    ),
    sporc_timeseries = sporc_ts,
    admb_timeseries = admb$timeseries,
    comparison = comparison
  ),
  file.path(cfg$paths$comparison_dir, "sporc_vs_admb_2024.rds")
)
write.csv(
  comparison,
  file.path(cfg$paths$comparison_dir, "sporc_vs_admb_2024.csv"),
  row.names = FALSE
)

summary_table <- aggregate(
  abs(relative_difference) ~ quantity,
  data = comparison,
  FUN = function(x) max(x, na.rm = TRUE)
)
names(summary_table)[2] <- "max_abs_relative_difference"
print(summary_table)
