library(SPoRC)
library(yaml)

source("R/config.R")
source("R/build_inputs.R")

configs <- load_scenario_configs("config/base.yml", "config/scenarios.yml")
base_cfg <- read_config("config/base.yml")
if (isTRUE(base_cfg$fit$base_only)) {
  configs <- configs[vapply(configs, function(cfg) cfg$scenario_name == "base", logical(1))]
  if (length(configs) == 0) {
    stop("Base scenario not found in config/scenarios.yml")
  }
}
data <- readRDS(base_cfg$paths$data_rds)

models <- vector("list", length(configs))
names(models) <- vapply(configs, function(cfg) cfg$scenario_name, "")

out_dir <- base_cfg$paths$outputs_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(configs)) {
  cfg <- configs[[i]]

  input_list <- build_pollock_inputs(cfg, data)
  data_list <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  model <- fit_model(
    data_list,
    parameters,
    mapping,
    random = cfg$selectivity_random_effects$random_effects,
    newton_loops = cfg$fit$newton_loops,
    silent = cfg$fit$silent
  )

  model$sdrep <- RTMB::sdreport(model)
  models[[i]] <- model

  saveRDS(model, file.path(out_dir, paste0(cfg$scenario_name, ".rds")))
}


saveRDS(models, file.path(out_dir, "all_models.rds"))
