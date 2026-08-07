args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."), mustWork = TRUE)
} else normalizePath(getwd(), mustWork = TRUE)

suppressPackageStartupMessages({ library(SPoRC); library(yaml) })
source(file.path(repo_root, "R", "build_inputs.R"))
cfg_file <- file.path(repo_root, "config", "base.yml")
cfg <- yaml::read_yaml(cfg_file)
data_file <- file.path(repo_root, cfg$paths$data_rds)
base_file <- file.path(repo_root, cfg$paths$outputs_dir, "base.rds")
if (!file.exists(base_file)) stop("Missing saved SPoRC fit: ", base_file)

raw <- readRDS(data_file)
inputs <- build_pollock_inputs(cfg, raw)
fit <- readRDS(base_file)

slice_comp <- function(x, fleet, survey = FALSE, drop_age1 = FALSE) {
  out <- if (survey) x[1, , , 1, fleet, drop = TRUE] else x[1, , , 1, 1, drop = TRUE]
  out <- matrix(out, nrow = length(raw$years), ncol = length(raw$ages))
  if (drop_age1) out <- out[, -1, drop = FALSE]
  out
}
normalize_rows <- function(x) sweep(x, 1, rowSums(x, na.rm = TRUE), "/")
make_fleet <- function(obs, exp, use, n_eff, index) {
  keep <- as.numeric(use) > 0
  out <- list(
    obs = normalize_rows(obs[keep, , drop = FALSE]),
    exp = normalize_rows(exp[keep, , drop = FALSE]),
    n_eff = as.numeric(n_eff)[keep], index = index, years = raw$years[keep]
  )
  totals <- rowSums(out$obs, na.rm = TRUE)
  rounded <- round(sweep(out$obs, 1, out$n_eff / totals, `*`), 0)
  valid <- is.finite(out$n_eff) & out$n_eff >= 1 & totals > 0 &
    rowSums(rounded, na.rm = TRUE) >= 1 & apply(out$exp, 1, function(z) all(is.finite(z)))
  list(data = lapply(out, function(z) if (is.matrix(z)) z[valid, , drop = FALSE]
    else if (length(z) == length(valid)) z[valid] else z), dropped = out$years[!valid])
}

fish <- make_fleet(
  slice_comp(inputs$data$ObsFishAgeComps, 1), slice_comp(fit$rep$CAA, 1),
  inputs$data$UseFishAgeComps[1, , 1], inputs$data$ISS_FishAgeComps[1, , 1, 1], raw$ages
)
bts <- make_fleet(
  slice_comp(inputs$data$ObsSrvAgeComps, 1, TRUE, TRUE),
  slice_comp(fit$rep$SrvIAA, 1, TRUE, TRUE), inputs$data$UseSrvAgeComps[1, , 1],
  inputs$data$ISS_SrvAgeComps[1, , 1, 1], raw$ages[-1]
)
ats <- make_fleet(
  slice_comp(inputs$data$ObsSrvAgeComps, 2, TRUE, TRUE),
  slice_comp(fit$rep$SrvIAA, 2, TRUE, TRUE), inputs$data$UseSrvAgeComps[1, , 2],
  inputs$data$ISS_SrvAgeComps[1, , 1, 2], raw$ages[-1]
)
osa_inputs <- list(Fishery = fish$data, BTS = bts$data, ATS = ats$data)
attr(osa_inputs, "lineage") <- list(
  base_file = normalizePath(base_file), base_md5 = unname(tools::md5sum(base_file)),
  data_file = normalizePath(data_file), data_md5 = unname(tools::md5sum(data_file)),
  config_file = normalizePath(cfg_file), config_md5 = unname(tools::md5sum(cfg_file)),
  prepared = Sys.time(), r_version = R.version.string,
  SPoRC_version = as.character(utils::packageVersion("SPoRC")),
  excluded_zero_information_years = list(Fishery = fish$dropped, BTS = bts$dropped, ATS = ats$dropped)
)
output_file <- file.path(repo_root, cfg$paths$outputs_dir, "osa_inputs.rds")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(osa_inputs, output_file)
message("Wrote OSA inputs to ", output_file)
