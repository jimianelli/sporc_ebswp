build_pollock_inputs <- function(cfg, data) {
  input_list <- Setup_Mod_Dim(
    years = data$years,
    ages = data$ages,
    lens = NA,
    n_regions = cfg$dimensions$n_regions,
    n_sexes = cfg$dimensions$n_sexes,
    n_fish_fleets = cfg$dimensions$n_fish_fleets,
    n_srv_fleets = cfg$dimensions$n_srv_fleets,
    verbose = FALSE
  )

  steepness_par <- cfg$recruitment$steepness_h
  if (isTRUE(cfg$recruitment$steepness_is_h)) {
    steepness_par <- qlogis((steepness_par - 0.2) / 0.8)
  }

  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = cfg$recruitment$do_rec_bias_ramp,
    sigmaR_switch = cfg$recruitment$sigmaR_switch,
    ln_sigmaR = cfg$recruitment$ln_sigmaR,
    rec_model = cfg$recruitment$rec_model,
    steepness_h = steepness_par,
    h_spec = cfg$recruitment$h_spec,
    sigmaR_spec = cfg$recruitment$sigmaR_spec,
    init_age_strc = cfg$recruitment$init_age_strc,
    ln_global_R0 = cfg$recruitment$ln_global_R0,
    t_spawn = cfg$recruitment$t_spawn,
    equil_init_age_strc = cfg$recruitment$equil_init_age_strc
  )

  if (isTRUE(cfg$recruitment$free_init_age_composition)) {
    init_dev_map <- input_list$par$ln_InitDevs
    init_dev_map[] <- seq_along(init_dev_map)
    input_list$map$ln_InitDevs <- factor(init_dev_map)
    input_list$data$equil_init_age_strc <- 0
  }

  fix_natmort <- array(
    0,
    dim = c(
      input_list$data$n_regions,
      length(input_list$data$years),
      length(input_list$data$ages),
      1
    )
  )
  fix_natmort[, , 1, ] <- cfg$biology$fixed_M$age1
  fix_natmort[, , 2, ] <- cfg$biology$fixed_M$age2
  if (length(input_list$data$ages) > 2) {
    fix_natmort[, , -(1:2), ] <- cfg$biology$fixed_M$age3plus
  }

  waa_fish <- if (!is.null(data$WAA_fish)) data$WAA_fish else NULL
  waa_srv <- if (!is.null(data$WAA_srv)) data$WAA_srv else NULL
  ageing_error <- if (!is.null(data$AgeingError)) data$AgeingError else NULL
  ageing_error_fish <- if (!is.null(data$AgeingError_fish)) data$AgeingError_fish else ageing_error
  ageing_error_srv <- if (!is.null(data$AgeingError_srv)) data$AgeingError_srv else ageing_error

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = data$WAA,
    WAA_fish = waa_fish,
    WAA_srv = waa_srv,
    MatAA = data$MatAA,
    AgeingError_fish = ageing_error_fish,
    AgeingError_srv = ageing_error_srv,
    AgeingError = ageing_error,
    fit_lengths = cfg$biology$fit_lengths,
    M_spec = cfg$biology$M_spec,
    Fixed_natmort = fix_natmort
  )

  use_fixed_movement <- isTRUE(as.logical(cfg$movement$use_fixed_movement))
  fixed_movement <- cfg$movement$Fixed_Movement
  expected_dim <- c(
    input_list$data$n_regions,
    input_list$data$n_regions,
    length(input_list$data$years),
    length(input_list$data$ages),
    input_list$data$n_sexes
  )

  has_valid_dim <- function(x, expected) {
    if (is.null(x) || is.null(dim(x))) {
      return(FALSE)
    }
    dims <- dim(x)
    length(dims) == length(expected) && all(dims == expected)
  }

  if (use_fixed_movement && !has_valid_dim(fixed_movement, expected_dim)) {
    fixed_movement <- array(0, dim = expected_dim)
    for (r in seq_len(expected_dim[1])) {
      fixed_movement[r, r, , , ] <- 1
    }
  }

  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = use_fixed_movement,
    Fixed_Movement = fixed_movement,
    do_recruits_move = cfg$movement$do_recruits_move
  )

  input_list <- Setup_Mod_Tagging(
    input_list = input_list,
    UseTagging = cfg$tagging$UseTagging
  )

  sigmaC <- cfg$catch$ln_sigmaC
  if (is.null(sigmaC)) {
    sigmaC <- log(cfg$catch$sigmaC)
  }
  ln_sigmaC <- array(
    sigmaC,
    dim = c(
      input_list$data$n_regions,
      length(input_list$data$years),
      input_list$data$n_fish_fleets
    )
  )

  ensure_region_year_fleet <- function(x, n_regions, n_years, n_fleets, label) {
    if (is.null(x)) {
      return(NULL)
    }
    dims <- dim(x)
    if (!is.null(dims) && length(dims) == 3 &&
      all(dims == c(n_regions, n_years, n_fleets))) {
      return(x)
    }
    if (!is.null(dims) && length(dims) == 2 && all(dims == c(n_years, n_fleets))) {
      return(array(x, dim = c(n_regions, n_years, n_fleets)))
    }
    if (!is.null(dims) && length(dims) == 3) {
      x <- drop(x)
    }
    vec <- as.numeric(x)
    if (length(vec) == n_regions * n_years * n_fleets) {
      return(array(vec, dim = c(n_regions, n_years, n_fleets)))
    }
    if (length(vec) == n_years * n_fleets) {
      return(array(vec, dim = c(n_regions, n_years, n_fleets)))
    }
    warning(
      label,
      " has length ",
      length(vec),
      " but expected ",
      n_regions * n_years * n_fleets
    )
    array(NA_real_, dim = c(n_regions, n_years, n_fleets))
  }

  ensure_year_fleet <- function(x, n_years, n_fleets, label) {
    if (is.null(x)) {
      return(NULL)
    }
    dims <- dim(x)
    if (!is.null(dims) && length(dims) == 2 && all(dims == c(n_years, n_fleets))) {
      return(x)
    }
    if (!is.null(dims) && length(dims) == 3) {
      x <- drop(x)
    }
    vec <- as.numeric(x)
    if (length(vec) == n_years * n_fleets) {
      return(matrix(vec, nrow = n_years, ncol = n_fleets))
    }
    warning(label, " has length ", length(vec), " but expected ", n_years * n_fleets)
    matrix(NA_real_, nrow = n_years, ncol = n_fleets)
  }

  n_years <- length(input_list$data$years)
  n_regions <- input_list$data$n_regions
  n_fleets <- input_list$data$n_fish_fleets
  obs_catch <- ensure_region_year_fleet(
    data$ObsCatch,
    n_regions,
    n_years,
    n_fleets,
    "ObsCatch"
  )
  catch_type <- ensure_year_fleet(data$Catch_Type, n_years, n_fleets, "Catch_Type")
  use_catch <- ensure_region_year_fleet(
    data$UseCatch,
    n_regions,
    n_years,
    n_fleets,
    "UseCatch"
  )

  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = obs_catch,
    Catch_Type = catch_type,
    UseCatch = use_catch,
    Use_F_pen = cfg$catch$Use_F_pen,
    sigmaC_spec = cfg$catch$sigmaC_spec,
    ln_sigmaC = ln_sigmaC
  )

  fishage_comp_agg_type <- cfg$fishery$FishAge_comp_agg_type
  if (is.null(fishage_comp_agg_type)) {
    fishage_comp_agg_type <- rep(1, input_list$data$n_fish_fleets)
  }
  fishlen_comp_agg_type <- cfg$fishery$FishLen_comp_agg_type
  if (is.null(fishlen_comp_agg_type)) {
    fishlen_comp_agg_type <- rep(0, input_list$data$n_fish_fleets)
  }

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = data$ObsFishIdx,
    ObsFishIdx_SE = data$ObsFishIdx_SE,
    UseFishIdx = data$UseFishIdx,
    ObsFishAgeComps = data$ObsFishAgeComps,
    UseFishAgeComps = data$UseFishAgeComps,
    ISS_FishAgeComps = data$ISS_FishAgeComps,
    ObsFishLenComps = data$ObsFishLenComps,
    UseFishLenComps = data$UseFishLenComps,
    ISS_FishLenComps = data$ISS_FishLenComps,
    fish_idx_type = cfg$fishery$fish_idx_type,
    FishAgeComps_LikeType = cfg$fishery$FishAgeComps_LikeType,
    FishLenComps_LikeType = cfg$fishery$FishLenComps_LikeType,
    FishAgeComps_Type = cfg$fishery$FishAgeComps_Type,
    FishLenComps_Type = cfg$fishery$FishLenComps_Type,
    FishAge_comp_agg_type = fishage_comp_agg_type,
    FishLen_comp_agg_type = fishlen_comp_agg_type
  )

  srvage_comp_agg_type <- cfg$survey$SrvAge_comp_agg_type
  if (is.null(srvage_comp_agg_type)) {
    srvage_comp_agg_type <- rep(1, input_list$data$n_srv_fleets)
  }
  srvlen_comp_agg_type <- cfg$survey$SrvLen_comp_agg_type
  if (is.null(srvlen_comp_agg_type)) {
    srvlen_comp_agg_type <- rep(0, input_list$data$n_srv_fleets)
  }

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = data$ObsSrvIdx,
    ObsSrvIdx_SE = data$ObsSrvIdx_SE,
    UseSrvIdx = data$UseSrvIdx,
    ObsSrvAgeComps = data$ObsSrvAgeComps,
    ISS_SrvAgeComps = data$ISS_SrvAgeComps,
    UseSrvAgeComps = data$UseSrvAgeComps,
    ObsSrvLenComps = data$ObsSrvLenComps,
    UseSrvLenComps = data$UseSrvLenComps,
    ISS_SrvLenComps = data$ISS_SrvLenComps,
    srv_idx_type = cfg$survey$srv_idx_type,
    SrvAgeComps_LikeType = cfg$survey$SrvAgeComps_LikeType,
    SrvLenComps_LikeType = cfg$survey$SrvLenComps_LikeType,
    SrvAgeComps_Type = cfg$survey$SrvAgeComps_Type,
    SrvLenComps_Type = cfg$survey$SrvLenComps_Type,
    SrvAge_comp_agg_type = srvage_comp_agg_type,
    SrvLen_comp_agg_type = srvlen_comp_agg_type
  )

  corr_opt_fish_semipar <- cfg$selectivity_random_effects$corr_opt_fish_semipar
  if (is.null(corr_opt_fish_semipar)) {
    corr_opt_fish_semipar <- cfg$selectivity_random_effects$corr_opt_semipar
  }

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    cont_tv_fish_sel = cfg$selectivity_random_effects$cont_tv_fish_sel,
    fishsel_pe_pars_spec = cfg$selectivity_random_effects$fishsel_pe_pars_spec,
    fish_sel_devs_spec = cfg$selectivity_random_effects$fish_sel_devs_spec,
    corr_opt_semipar = corr_opt_fish_semipar,
    fish_sel_blocks = cfg$fishery_selectivity$fish_sel_blocks,
    fish_sel_model = cfg$fishery_selectivity$fish_sel_model,
    fish_q_blocks = cfg$fishery_selectivity$fish_q_blocks,
    fish_fixed_sel_pars_spec = cfg$fishery_selectivity$fish_fixed_sel_pars_spec,
    fish_q_spec = cfg$fishery_selectivity$fish_q_spec
  )

  if (!is.null(cfg$fishery_selectivity$fish_fixed_sel_pars_start) &&
      !is.null(input_list$par$ln_fish_fixed_sel_pars)) {
    values <- as.numeric(unlist(cfg$fishery_selectivity$fish_fixed_sel_pars_start))
    n_values <- min(length(values), dim(input_list$par$ln_fish_fixed_sel_pars)[2])
    input_list$par$ln_fish_fixed_sel_pars[1, seq_len(n_values), 1, 1, 1] <-
      values[seq_len(n_values)]
  }

  extract_fleet_id <- function(x) {
    match <- regexec("Fleet_(\\d+)$", x)
    reg <- regmatches(x, match)
    if (length(reg) == 0 || length(reg[[1]]) < 2) {
      return(NA_integer_)
    }
    as.integer(reg[[1]][2])
  }

  cont_tv_srv_sel <- cfg$selectivity_random_effects$cont_tv_srv_sel
  if (is.null(cont_tv_srv_sel)) {
    cont_tv_srv_sel <- character(0)
  }
  cont_tv_srv_sel <- as.character(cont_tv_srv_sel)
  fleet_ids <- vapply(cont_tv_srv_sel, extract_fleet_id, integer(1))
  missing_fleets <- setdiff(
    seq_len(input_list$data$n_srv_fleets),
    fleet_ids[!is.na(fleet_ids)]
  )
  if (length(missing_fleets) > 0) {
    cont_tv_srv_sel <- c(
      cont_tv_srv_sel,
      paste0("none_Fleet_", missing_fleets)
    )
  }

  srvsel_pe_pars_spec <- cfg$survey$srvsel_pe_pars_spec
  if (is.null(srvsel_pe_pars_spec)) {
    srvsel_pe_pars_spec <- cfg$selectivity_random_effects$srvsel_pe_pars_spec
  }
  if (is.null(srvsel_pe_pars_spec)) {
    srvsel_pe_pars_spec <- rep("none", input_list$data$n_srv_fleets)
  }
  corr_opt_srv_semipar <- cfg$selectivity_random_effects$corr_opt_srv_semipar
  if (is.null(corr_opt_srv_semipar)) {
    corr_opt_srv_semipar <- cfg$selectivity_random_effects$corr_opt_semipar
  }

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    cont_tv_srv_sel = cont_tv_srv_sel,
    srv_sel_blocks = cfg$survey$srv_sel_blocks,
    srv_sel_model = cfg$survey$srv_sel_model,
    srv_q_blocks = cfg$survey$srv_q_blocks,
    srvsel_pe_pars_spec = srvsel_pe_pars_spec,
    srv_fixed_sel_pars_spec = cfg$survey$srv_fixed_sel_pars_spec,
    srv_q_spec = cfg$survey$srv_q_spec,
    srv_sel_devs_spec = cfg$selectivity_random_effects$srv_sel_devs_spec,
    corr_opt_semipar = corr_opt_srv_semipar
  )

  if (!is.null(cfg$survey$srv_fixed_sel_pars_start) &&
      !is.null(input_list$par$ln_srv_fixed_sel_pars)) {
    srv_fleet_names <- c("BTS", "ATS", "AVO")
    for (fleet_name in names(cfg$survey$srv_fixed_sel_pars_start)) {
      fleet <- match(fleet_name, srv_fleet_names)
      if (is.na(fleet)) {
        stop("Unknown survey fleet in srv_fixed_sel_pars_start: ", fleet_name)
      }
      values <- as.numeric(unlist(cfg$survey$srv_fixed_sel_pars_start[[fleet_name]]))
      n_values <- min(length(values), dim(input_list$par$ln_srv_fixed_sel_pars)[2])
      input_list$par$ln_srv_fixed_sel_pars[1, seq_len(n_values), 1, 1, fleet] <-
        log(values[seq_len(n_values)])
    }
  }

  if (identical(cfg$selectivity_random_effects$fishsel_pe_pars_spec, "fix") &&
      !is.null(cfg$selectivity_random_effects$fishsel_pe_pars_log_sigma) &&
      !is.null(input_list$par$fishsel_pe_pars)) {
    input_list$par$fishsel_pe_pars[, 4, , ] <-
      log(cfg$selectivity_random_effects$fishsel_pe_pars_log_sigma)
  }
  if (any(srvsel_pe_pars_spec == "fix") &&
      !is.null(cfg$selectivity_random_effects$srvsel_bts_logistic_log_sigma) &&
      !is.null(input_list$par$srvsel_pe_pars)) {
    input_list$par$srvsel_pe_pars[, 1:2, , 1] <-
      log(cfg$selectivity_random_effects$srvsel_bts_logistic_log_sigma)
  }
  if (any(srvsel_pe_pars_spec == "fix") &&
      !is.null(cfg$selectivity_random_effects$srvsel_ats_avo_log_sigma) &&
      !is.null(input_list$par$srvsel_pe_pars)) {
    input_list$par$srvsel_pe_pars[, 4, , 2] <-
      log(cfg$selectivity_random_effects$srvsel_ats_avo_log_sigma)
  }

  make_weight_array <- function(weight, n_regions, n_years, n_sexes, n_fleets) {
    array(weight, dim = c(n_regions, n_years, n_sexes, n_fleets))
  }

  n_regions <- input_list$data$n_regions
  n_years <- length(input_list$data$years)
  n_sexes <- input_list$data$n_sexes
  n_fish_fleets <- input_list$data$n_fish_fleets
  n_srv_fleets <- input_list$data$n_srv_fleets

  wt_fish_age <- make_weight_array(
    cfg$weights$Wt_FishAgeComps,
    n_regions,
    n_years,
    n_sexes,
    n_fish_fleets
  )
  wt_fish_len <- make_weight_array(
    cfg$weights$Wt_FishLenComps,
    n_regions,
    n_years,
    n_sexes,
    n_fish_fleets
  )
  wt_srv_age <- make_weight_array(
    cfg$weights$Wt_SrvAgeComps,
    n_regions,
    n_years,
    n_sexes,
    n_srv_fleets
  )
  wt_srv_len <- make_weight_array(
    cfg$weights$Wt_SrvLenComps,
    n_regions,
    n_years,
    n_sexes,
    n_srv_fleets
  )

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = cfg$weights$Wt_Catch,
    Wt_FishIdx = cfg$weights$Wt_FishIdx,
    Wt_SrvIdx = cfg$weights$Wt_SrvIdx,
    Wt_Rec = cfg$weights$Wt_Rec,
    Wt_F = cfg$weights$Wt_F,
    Wt_Tagging = cfg$weights$Wt_Tagging,
    Wt_FishAgeComps = wt_fish_age,
    Wt_FishLenComps = wt_fish_len,
    Wt_SrvAgeComps = wt_srv_age,
    Wt_SrvLenComps = wt_srv_len
  )

  input_list
}
