# Build the data list expected by SPoRC and save to data/ebs_pollock.rds.
# Default uses the 2024 ADMB base-run data file.

suppressPackageStartupMessages({
  library(yaml)
})

parse_admb_dat <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  lines <- lines[!grepl("^#[[:space:]]*[-0-9.]", lines)]

  headers <- grep("^#", lines)
  if (length(headers) == 0) {
    stop("No headers found in ", path)
  }

  names <- sub("^#", "", lines[headers])
  names <- trimws(names)
  values <- vector("list", length(names))

  for (i in seq_along(headers)) {
    start <- headers[i] + 1
    end <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
    block <- lines[start:end]
    if (length(block) == 0) {
      values[[i]] <- numeric()
      next
    }
    values[[i]] <- scan(text = block, quiet = TRUE)
  }

  setNames(values, names)
}

as_matrix <- function(values, n_rows, n_cols, label) {
  if (length(values) != n_rows * n_cols) {
    stop(label, " has ", length(values), " values; expected ", n_rows * n_cols)
  }
  matrix(values, nrow = n_rows, ncol = n_cols, byrow = TRUE)
}

match_years <- function(src_years, years, label) {
  idx <- match(src_years, years)
  if (any(is.na(idx))) {
    missing <- src_years[is.na(idx)]
    stop("Missing years for ", label, ": ", paste(missing, collapse = ", "))
  }
  idx
}

cfg <- yaml::read_yaml("config/base.yml")
source_root <- Sys.getenv("EBS_POLLOCK_ADMB", unset = cfg$paths$admb_source_root)
run_dir <- file.path(source_root,  "runs", "2024")
pm_dat_path <- file.path(run_dir, "pm.dat")

if (!file.exists(pm_dat_path)) {
  stop("pm.dat not found at ", pm_dat_path)
}

pm_lines <- readLines(pm_dat_path, warn = FALSE)
pm_lines <- pm_lines[pm_lines != ""]
if (length(pm_lines) < 2) {
  stop("pm.dat does not contain a data file reference")
}

data_path <- normalizePath(file.path(run_dir, pm_lines[2]), mustWork = TRUE)
message("Using ADMB data file: ", data_path)

admb <- parse_admb_dat(data_path)

styr <- as.integer(admb$styr)
endyr <- as.integer(admb$endyr)
years <- seq(styr, endyr)
n_years <- length(years)

n_ages <- as.integer(admb$nages)
ages <- seq_len(n_ages)

waa_mat <- as_matrix(admb$wt_ssb, n_years, n_ages, "wt_ssb")
WAA <- array(waa_mat, dim = c(1, n_years, n_ages, 1))

mat_mat <- matrix(rep(admb$p_mature, n_years), nrow = n_years, byrow = TRUE)
MatAA <- array(mat_mat, dim = c(1, n_years, n_ages, 1))

age_err <- as_matrix(admb$age_err, n_ages, n_ages, "age_err")

waa_fish_mat <- as_matrix(admb$wt_fsh, n_years, n_ages, "wt_fsh")
WAA_fish <- array(waa_fish_mat, dim = c(1, n_years, n_ages, 1, 1))

WAA_srv <- array(0, dim = c(1, n_years, n_ages, 1, 3))
for (f in 1:3) {
  WAA_srv[,,,, f] <- WAA
}

bts_years <- as.integer(admb$yrs_bts_data)
ats_years <- as.integer(admb$yrs_ats_data)
avo_years <- as.integer(admb$yrs_avo)

bts_idx <- match_years(bts_years, years, "BTS WAA")
ats_idx <- match_years(ats_years, years, "ATS WAA")
avo_idx <- match_years(avo_years, years, "AVO WAA")

WAA_srv[1, bts_idx, , 1, 1] <- as_matrix(admb$wt_bts, length(bts_years), n_ages, "wt_bts")
WAA_srv[1, ats_idx, , 1, 2] <- as_matrix(admb$wt_ats, length(ats_years), n_ages, "wt_ats")
WAA_srv[1, avo_idx, , 1, 3] <- as_matrix(admb$wt_avo, length(avo_years), n_ages, "wt_avo")

obs_catch <- rep(NA_real_, n_years)
if (length(admb$obs_catch) == n_years) {
  obs_catch <- admb$obs_catch
} else if (length(admb$obs_catch) == length(admb$yrs_fsh_data)) {
  fsh_idx <- match_years(as.integer(admb$yrs_fsh_data), years, "fishery catch")
  obs_catch[fsh_idx] <- admb$obs_catch
} else {
  stop("obs_catch length does not match years or yrs_fsh_data")
}

ObsCatch <- array(obs_catch, dim = c(1, n_years, 1))
Catch_Type <- array(1, dim = c(1, n_years, 1))
UseCatch <- array(ifelse(is.na(obs_catch), 0, 1), dim = c(1, n_years, 1))

ObsFishIdx <- array(NA_real_, dim = c(1, n_years, 1))
ObsFishIdx_SE <- array(NA_real_, dim = c(1, n_years, 1))
UseFishIdx <- array(0, dim = c(1, n_years, 1))
if (length(admb$yrs_cpue) > 0) {
  cpue_idx <- match_years(as.integer(admb$yrs_cpue), years, "CPUE")
  ObsFishIdx[1, cpue_idx, 1] <- admb$obs_cpue
  ObsFishIdx_SE[1, cpue_idx, 1] <- admb$obs_cpue_std
  UseFishIdx[1, cpue_idx, 1] <- 1
}

fsh_years <- as.integer(admb$yrs_fsh_data)
fsh_idx <- match_years(fsh_years, years, "fishery age comps")
obs_fsh_age <- as_matrix(admb$oac_fsh_data, length(fsh_years), n_ages, "oac_fsh_data")

ObsFishAgeComps <- array(NA_real_, dim = c(1, n_years, n_ages, 1, 1))
ObsFishAgeComps[1, fsh_idx, , 1, 1] <- obs_fsh_age

UseFishAgeComps <- array(0, dim = c(1, n_years, 1))
UseFishAgeComps[1, fsh_idx, 1] <- 1

ISS_FishAgeComps <- array(NA_real_, dim = c(1, n_years, 1, 1))
ISS_FishAgeComps[1, fsh_idx, 1, 1] <- admb$sam_fsh

n_lens <- 1L
ObsFishLenComps <- array(NA_real_, dim = c(1, n_years, n_lens, 1, 1))
UseFishLenComps <- array(0, dim = c(1, n_years, 1))
ISS_FishLenComps <- NULL

ObsSrvIdx <- array(NA_real_, dim = c(1, n_years, 3))
ObsSrvIdx_SE <- array(NA_real_, dim = c(1, n_years, 3))
UseSrvIdx <- array(0, dim = c(1, n_years, 3))

ObsSrvIdx[1, bts_idx, 1] <- admb$ob_bts
ObsSrvIdx_SE[1, bts_idx, 1] <- admb$ob_bts_std
UseSrvIdx[1, bts_idx, 1] <- 1

ObsSrvIdx[1, ats_idx, 2] <- admb$ob_ats
ObsSrvIdx_SE[1, ats_idx, 2] <- admb$ob_ats_std
UseSrvIdx[1, ats_idx, 2] <- 1

ObsSrvIdx[1, avo_idx, 3] <- admb$ob_avo
ObsSrvIdx_SE[1, avo_idx, 3] <- admb$ob_avo_std
UseSrvIdx[1, avo_idx, 3] <- 1

ObsSrvAgeComps <- array(NA_real_, dim = c(1, n_years, n_ages, 1, 3))
UseSrvAgeComps <- array(0, dim = c(1, n_years, 3))
ISS_SrvAgeComps <- array(NA_real_, dim = c(1, n_years, 1, 3))

obs_bts_age <- as_matrix(admb$oac_bts, length(bts_years), n_ages, "oac_bts")
ObsSrvAgeComps[1, bts_idx, , 1, 1] <- obs_bts_age
UseSrvAgeComps[1, bts_idx, 1] <- 1
ISS_SrvAgeComps[1, bts_idx, 1, 1] <- admb$sam_bts

obs_ats_age <- as_matrix(admb$oac_ats, length(ats_years), n_ages, "oac_ats")
ObsSrvAgeComps[1, ats_idx, , 1, 2] <- obs_ats_age
UseSrvAgeComps[1, ats_idx, 2] <- 1
ISS_SrvAgeComps[1, ats_idx, 1, 2] <- admb$sam_ats

ObsSrvLenComps <- array(NA_real_, dim = c(1, n_years, n_lens, 1, 3))
UseSrvLenComps <- array(0, dim = c(1, n_years, 3))
ISS_SrvLenComps <- NULL

ebs_data <- list(
  years = years,
  ages = ages,
  WAA = WAA,
  WAA_fish = WAA_fish,
  WAA_srv = WAA_srv,
  MatAA = MatAA,
  AgeingError_fish = age_err,
  AgeingError_srv = age_err,
  AgeingError = age_err,
  ObsCatch = ObsCatch,
  Catch_Type = Catch_Type,
  UseCatch = UseCatch,
  ObsFishIdx = ObsFishIdx,
  ObsFishIdx_SE = ObsFishIdx_SE,
  UseFishIdx = UseFishIdx,
  ObsFishAgeComps = ObsFishAgeComps,
  UseFishAgeComps = UseFishAgeComps,
  ISS_FishAgeComps = ISS_FishAgeComps,
  ObsFishLenComps = ObsFishLenComps,
  UseFishLenComps = UseFishLenComps,
  ISS_FishLenComps = ISS_FishLenComps,
  ObsSrvIdx = ObsSrvIdx,
  ObsSrvIdx_SE = ObsSrvIdx_SE,
  UseSrvIdx = UseSrvIdx,
  ObsSrvAgeComps = ObsSrvAgeComps,
  ISS_SrvAgeComps = ISS_SrvAgeComps,
  UseSrvAgeComps = UseSrvAgeComps,
  ObsSrvLenComps = ObsSrvLenComps,
  UseSrvLenComps = UseSrvLenComps,
  ISS_SrvLenComps = ISS_SrvLenComps
)

dir.create(dirname(cfg$paths$data_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(ebs_data, cfg$paths$data_rds)
message("Saved ", cfg$paths$data_rds)
