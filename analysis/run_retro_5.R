suppressPackageStartupMessages({
  library(SPoRC)
  library(yaml)
})

source("R/build_inputs.R")

cfg <- yaml::read_yaml("config/base.yml")
raw_data <- readRDS(cfg$paths$data_rds)
input_list <- build_pollock_inputs(cfg, raw_data)

retro <- SPoRC::do_retrospective(
  n_retro = 5,
  data = input_list$data,
  parameters = input_list$par,
  mapping = input_list$map,
  random = cfg$selectivity_random_effects$random_effects,
  do_par = FALSE,
  n_cores = 1,
  newton_loops = cfg$fit$newton_loops,
  do_sdrep = FALSE
)

out_dir <- cfg$paths$outputs_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(retro, file.path(out_dir, "retro_5_peel.rds"))
