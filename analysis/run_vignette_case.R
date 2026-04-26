library(SPoRC)
library(RTMB)

data("sgl_rg_ebswp_data")

input_list <- Setup_Mod_Dim(
  years = sgl_rg_ebswp_data$years,
  ages = sgl_rg_ebswp_data$ages,
  lens = NA,
  n_regions = 1,
  n_sexes = 1,
  n_fish_fleets = 1,
  n_srv_fleets = 3
)

inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

input_list <- Setup_Mod_Rec(
  input_list = input_list,
  do_rec_bias_ramp = 0,
  sigmaR_switch = 1,
  ln_sigmaR = log(c(1, 1)),
  rec_model = "bh_rec",
  steepness_h = inv_steepness(0.623013),
  h_spec = "fix",
  sigmaR_spec = "fix",
  init_age_strc = 1,
  ln_global_R0 = 10,
  t_spawn = 0.25,
  equil_init_age_strc = 2
)

fix_natmort <- array(
  0,
  dim = c(
    input_list$data$n_regions,
    length(input_list$data$years),
    length(input_list$data$ages),
    1
  )
)
fix_natmort[, , 1, ] <- 0.9
fix_natmort[, , 2, ] <- 0.45
fix_natmort[, , -(1:2), ] <- 0.3

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,
  WAA = sgl_rg_ebswp_data$WAA,
  MatAA = sgl_rg_ebswp_data$MatAA,
  AgeingError = sgl_rg_ebswp_data$AgeingError,
  fit_lengths = 0,
  M_spec = "fix",
  Fixed_natmort = fix_natmort
)

input_list <- Setup_Mod_Movement(
  input_list = input_list,
  use_fixed_movement = 1,
  Fixed_Movement = NA,
  do_recruits_move = 0
)

input_list <- Setup_Mod_Tagging(
  input_list = input_list,
  UseTagging = 0
)

input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,
  ObsCatch = sgl_rg_ebswp_data$ObsCatch,
  Catch_Type = sgl_rg_ebswp_data$Catch_Type,
  UseCatch = sgl_rg_ebswp_data$UseCatch,
  Use_F_pen = 1,
  sigmaC_spec = "fix",
  ln_sigmaC = array(log(0.05), dim = c(1, length(input_list$data$years), 1))
)

input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,
  ObsFishIdx = sgl_rg_ebswp_data$ObsFishIdx,
  ObsFishIdx_SE = sgl_rg_ebswp_data$ObsFishIdx_SE,
  UseFishIdx = sgl_rg_ebswp_data$UseFishIdx,
  ObsFishAgeComps = sgl_rg_ebswp_data$ObsFishAgeComps,
  UseFishAgeComps = sgl_rg_ebswp_data$UseFishAgeComps,
  ISS_FishAgeComps = sgl_rg_ebswp_data$ISS_FishAgeComps,
  ISS_FishLenComps = NULL,
  ObsFishLenComps = array(
    NA_real_,
    dim = c(1, length(input_list$data$years), length(input_list$data$lens), 1, 1)
  ),
  UseFishLenComps = array(0, dim = c(1, length(input_list$data$years), 1)),
  fish_idx_type = c("biom"),
  FishAgeComps_LikeType = c("Multinomial"),
  FishLenComps_LikeType = c("none"),
  FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
  FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
)

input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,
  ObsSrvIdx = sgl_rg_ebswp_data$ObsSrvIdx,
  ObsSrvIdx_SE = sgl_rg_ebswp_data$ObsSrvIdx_SE,
  UseSrvIdx = sgl_rg_ebswp_data$UseSrvIdx,
  ObsSrvAgeComps = sgl_rg_ebswp_data$ObsSrvAgeComps,
  ISS_SrvAgeComps = sgl_rg_ebswp_data$ISS_SrvAgeComps,
  UseSrvAgeComps = sgl_rg_ebswp_data$UseSrvAgeComps,
  ObsSrvLenComps = array(
    NA_real_,
    dim = c(1, length(input_list$data$years), length(input_list$data$lens), 1, 3)
  ),
  UseSrvLenComps = array(0, dim = c(1, length(input_list$data$years), 3)),
  ISS_SrvLenComps = NULL,
  srv_idx_type = c("biom", "biom", "biom"),
  SrvAgeComps_LikeType = c("Multinomial", "Multinomial", "none"),
  SrvLenComps_LikeType = c("none", "none", "none"),
  SrvAgeComps_Type = c(
    "agg_Year_1-terminal_Fleet_1",
    "agg_Year_1-terminal_Fleet_2",
    "none_Year_1-terminal_Fleet_3"
  ),
  SrvLenComps_Type = c(
    "none_Year_1-terminal_Fleet_1",
    "none_Year_1-terminal_Fleet_2",
    "none_Year_1-terminal_Fleet_3"
  )
)

input_list <- Setup_Mod_Fishsel_and_Q(
  input_list = input_list,
  cont_tv_fish_sel = c("2dar1_Fleet_1"),
  fishsel_pe_pars_spec = "fix",
  fish_sel_devs_spec = "est_all",
  corr_opt_semipar = "corr_zero_y_b",
  fish_sel_blocks = c("none_Fleet_1"),
  fish_sel_model = c("logist1_Fleet_1"),
  fish_q_blocks = c("none_Fleet_1"),
  fish_fixed_sel_pars_spec = c("est_all"),
  fish_q_spec = c("est_all")
)

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,
  cont_tv_srv_sel = c("iid_Fleet_1", "2dar1_Fleet_2", "2dar1_Fleet_3"),
  srvsel_pe_pars_spec = c("fix", "fix", "fix"),
  srv_sel_devs_spec = c("est_all", "est_all", "est_shared_f_2"),
  corr_opt_semipar = c(NA, "corr_zero_y_b", "corr_zero_y_b"),
  srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_sel_model = c("logist1_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3"),
  srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_shared_f_2"),
  srv_q_spec = c("est_all", "est_all", "est_all")
)

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = 1,
  Wt_FishIdx = 1,
  Wt_SrvIdx = 1,
  Wt_Rec = 1,
  Wt_F = 1,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(
    1,
    dim = c(input_list$data$n_regions, length(input_list$data$years),
            input_list$data$n_sexes, input_list$data$n_fish_fleets)
  ),
  Wt_FishLenComps = array(
    1,
    dim = c(input_list$data$n_regions, length(input_list$data$years),
            input_list$data$n_sexes, input_list$data$n_fish_fleets)
  ),
  Wt_SrvAgeComps = array(
    1,
    dim = c(input_list$data$n_regions, length(input_list$data$years),
            input_list$data$n_sexes, input_list$data$n_srv_fleets)
  ),
  Wt_SrvLenComps = array(
    1,
    dim = c(input_list$data$n_regions, length(input_list$data$years),
            input_list$data$n_sexes, input_list$data$n_srv_fleets)
  )
)

data_list <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

parameters$fishsel_pe_pars[, 4, , ] <- log(0.075)
parameters$srvsel_pe_pars[, 1:2, , 1] <- log(0.075)
parameters$srvsel_pe_pars[, 4, , 2] <- log(0.15)

model <- fit_model(
  data_list,
  parameters,
  mapping,
  random = NULL,
  newton_loops = 3,
  silent = TRUE
)

model$sdrep <- RTMB::sdreport(model)
model$case <- list(
  name = "vign",
  source = "https://chengmatt.github.io/SPoRC/articles/f_single_region_ebs_pollock_case_study.html"
)

dir.create("analysis/outputs", recursive = TRUE, showWarnings = FALSE)
saveRDS(model, "analysis/outputs/vign.rds")
saveRDS(input_list, "analysis/outputs/vign_input_list.rds")
