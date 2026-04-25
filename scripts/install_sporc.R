library(yaml)

cfg <- yaml::read_yaml("config/base.yml")
ref <- cfg$sporc$ref

remotes::install_github("afsc-assessments/SPoRC", ref = ref)
