read_config <- function(path) {
  yaml::read_yaml(path)
}

apply_scenario <- function(base, scenario) {
  name <- scenario$name
  scenario$name <- NULL
  cfg <- modifyList(base, scenario, keep.null = TRUE)
  cfg$scenario_name <- name
  cfg
}

load_scenario_configs <- function(base_path, scenarios_path) {
  base <- read_config(base_path)
  scenario_data <- read_config(scenarios_path)
  scenarios <- scenario_data$scenarios

  if (is.null(scenarios) || length(scenarios) == 0) {
    stop("No scenarios found in ", scenarios_path)
  }

  lapply(scenarios, function(scenario) apply_scenario(base, scenario))
}
