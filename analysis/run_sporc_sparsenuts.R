#!/usr/bin/env Rscript

parse_arg <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[[1]])
}

parse_int_arg <- function(args, name, default) {
  value <- parse_arg(args, name, as.character(default))
  as.integer(value)
}

args <- commandArgs(trailingOnly = TRUE)

model_name <- parse_arg(args, "model-name", "sporc_base")
fit_stem <- parse_arg(args, "fit-stem", "base")
objective_mode <- parse_arg(args, "objective-mode", "saved")
cores <- parse_int_arg(args, "cores", 1L)
chains <- parse_int_arg(args, "chains", 4L)
num_samples <- parse_arg(args, "num-samples", NULL)
num_warmup <- parse_arg(args, "num-warmup", NULL)
default_output_stem <- if (identical(fit_stem, "base")) "sporc_sparsenuts_default" else paste0("sporc_sparsenuts_", fit_stem)
output_stem <- parse_arg(args, "output-stem", default_output_stem)
default_fig_prefix <- if (identical(fit_stem, "base")) "sporc_sparsenuts" else output_stem
fig_prefix <- parse_arg(args, "fig-prefix", default_fig_prefix)

suppressPackageStartupMessages({
  library(RTMB)
  library(SparseNUTS)
})

if (getRversion() < "4.6.0") {
  stop("SparseNUTS diagnostic MCMC must be run with R >= 4.6.0; found R ", getRversion(), ".")
}

fit_path <- file.path("analysis", "outputs", paste0(fit_stem, ".rds"))
if (!file.exists(fit_path)) {
  stop("Missing SPoRC fit: ", fit_path, ". Run analysis/run_scenarios.R or analysis/run_vignette_case.R first.")
}

out_dir <- file.path("analysis", "outputs", "sparsenuts")
fig_dir <- file.path("analysis", "outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, paste0(output_stem, ".rds"))

fit <- readRDS(fit_path)
saved_fit <- fit

sample_args <- list(
  obj = fit,
  cores = cores,
  chains = chains,
  model_name = model_name
)
if (!is.null(num_samples)) {
  sample_args$num_samples <- as.integer(num_samples)
}
if (!is.null(num_warmup)) {
  sample_args$num_warmup <- as.integer(num_warmup)
}

if (identical(objective_mode, "global-data")) {
  suppressPackageStartupMessages({
    library(SPoRC)
    library(yaml)
  })
  input_list_path <- file.path("analysis", "outputs", paste0(fit_stem, "_input_list.rds"))
  if (file.exists(input_list_path)) {
    input_list <- readRDS(input_list_path)
    random_effects <- NULL
  } else {
    source("R/config.R")
    source("R/build_inputs.R")
    cfg <- read_config("config/base.yml")
    raw_data <- readRDS(cfg$paths$data_rds)
    input_list <- build_pollock_inputs(cfg, raw_data)
    random_effects <- cfg$selectivity_random_effects$random_effects
  }

  SPoRC_rtmb <<- SPoRC:::SPoRC_rtmb
  .sporc_sparsenuts_data <<- input_list$data
  SPoRC_rtmb_global_data <<- function(pars) {
    SPoRC_rtmb(pars, .sporc_sparsenuts_data)
  }

  fit <- RTMB::MakeADFun(
    SPoRC_rtmb_global_data,
    parameters = input_list$par,
    map = input_list$map,
    random = random_effects,
    silent = TRUE
  )
  if (!identical(names(fit$par), names(saved_fit$par))) {
    stop("Rebuilt objective has different fixed parameter names than the saved fit.")
  }
  if (!identical(names(fit$env$last.par.best), names(saved_fit$env$last.par.best))) {
    stop("Rebuilt objective has different full parameter names than the saved fit.")
  }
  fit$par <- saved_fit$par
  fit$env$last.par.best <- saved_fit$env$last.par.best
  fit$optim <- saved_fit$optim

  sample_args$obj <- fit
  sample_args$globals <- list(
    SPoRC_rtmb = SPoRC_rtmb,
    .sporc_sparsenuts_data = .sporc_sparsenuts_data,
    SPoRC_rtmb_global_data = SPoRC_rtmb_global_data
  )
} else if (!identical(objective_mode, "saved")) {
  stop("Unknown --objective-mode: ", objective_mode)
}

snuts_fit <- do.call(SparseNUTS::sample_snuts, sample_args)

diagnostics <- SparseNUTS::check_snuts_diagnostics(snuts_fit, print = FALSE)
out <- list(
  snuts_fit = snuts_fit,
  diagnostics = diagnostics,
  sampler_settings = list(
    model_name = model_name,
    fit_stem = fit_stem,
    sampler_call = paste0(
      "SparseNUTS::sample_snuts(obj = fit, cores = ", cores,
      ", chains = ", chains, ", model_name = model_name)"
    ),
    objective_mode = objective_mode,
    defaults = "SparseNUTS package defaults for samples, warmup, metric, init, and control unless overridden by command-line arguments.",
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    sparsenuts_version = as.character(utils::packageVersion("SparseNUTS")),
    pairs_call = "graphics::pairs(snuts_fit, pars = match(pair_pars, available_pars), inc_warmup = FALSE, order = 'orig') dispatches to SparseNUTS::pairs.tmbfit.",
    plotting = "Pairs, native marginal, sampler-parameter, uncertainty, and Q figures use SparseNUTS plotting methods; interval and trace figures use base R."
  ),
  selected_parameters = NULL,
  figures = NULL
)
saveRDS(out, out_path)

monitor <- as.data.frame(snuts_fit$monitor)
monitor$variable <- snuts_fit$monitor$variable
monitor$rank_rhat <- ifelse(is.finite(monitor$rhat), monitor$rhat, Inf)
monitor$rank_ess <- ifelse(is.finite(monitor$ess_bulk), monitor$ess_bulk, Inf)
monitor <- monitor[order(-monitor$rank_rhat, monitor$rank_ess), ]
available_pars <- dimnames(snuts_fit$samples)[[3]]
monitor <- monitor[monitor$variable %in% available_pars, ]

pick_pars <- function(n) {
  out <- monitor$variable[seq_len(min(n, nrow(monitor)))]
  if (length(out) == 0) {
    out <- available_pars[seq_len(min(n, length(available_pars)))]
  }
  out
}

pair_pars <- pick_pars(6)
marginal_pars <- pick_pars(12)
trace_pars <- pair_pars

postwarmup <- snuts_fit$samples[-seq_len(snuts_fit$warmup), , , drop = FALSE]

figures <- list()

plot_interval_summaries <- function(samples, pars, prob = 0.8) {
  probs <- c((1 - prob) / 2, 0.5, 1 - (1 - prob) / 2)
  summary <- t(vapply(pars, function(par) {
    stats::quantile(as.vector(samples[, , par, drop = TRUE]), probs = probs, na.rm = TRUE)
  }, numeric(length(probs))))
  y <- seq_along(pars)
  graphics::plot(
    range(summary),
    range(y),
    type = "n",
    yaxt = "n",
    xlab = "Parameter value",
    ylab = "",
    main = paste0(round(100 * prob), "% posterior intervals")
  )
  graphics::axis(2, at = y, labels = pars, las = 2, cex.axis = 0.7)
  graphics::segments(summary[, 1], y, summary[, 3], y, lwd = 4, col = "grey70")
  graphics::points(summary[, 2], y, pch = 19)
}

plot_traces <- function(samples, pars) {
  old_par <- graphics::par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(graphics::par(old_par), add = TRUE)
  for (par in pars) {
    chain_samples <- samples[, , par, drop = TRUE]
    graphics::matplot(
      chain_samples,
      type = "l",
      lty = 1,
      lwd = 0.6,
      xlab = "Post-warmup iteration",
      ylab = "Value",
      main = par
    )
  }
  graphics::mtext("Trace plots", outer = TRUE, cex = 1.1)
}

figures$pairs <- file.path(fig_dir, paste0(fig_prefix, "_pairs.png"))
png(figures$pairs, width = 1800, height = 1800, res = 180)
graphics::pairs(
  snuts_fit,
  pars = match(pair_pars, available_pars),
  inc_warmup = FALSE,
  order = "orig"
)
dev.off()

figures$marginals <- file.path(fig_dir, paste0(fig_prefix, "_marginals.png"))
png(figures$marginals, width = 1800, height = 1400, res = 180)
plot_interval_summaries(postwarmup, marginal_pars, prob = 0.8)
dev.off()

figures$trace <- file.path(fig_dir, paste0(fig_prefix, "_trace.png"))
png(figures$trace, width = 1800, height = 1400, res = 180)
plot_traces(postwarmup, trace_pars)
dev.off()

figures$sampler_params <- file.path(fig_dir, paste0(fig_prefix, "_sampler_params.png"))
png(figures$sampler_params, width = 1800, height = 1400, res = 180)
SparseNUTS::plot_sampler_params(snuts_fit, plot = TRUE)
dev.off()

figures$native_marginals <- file.path(fig_dir, paste0(fig_prefix, "_native_marginals.png"))
png(figures$native_marginals, width = 1800, height = 1400, res = 180)
SparseNUTS::plot_marginals(
  snuts_fit,
  pars = match(marginal_pars, available_pars),
  mfrow = c(3, 4)
)
dev.off()

figures$uncertainties <- file.path(fig_dir, paste0(fig_prefix, "_uncertainties.png"))
png(figures$uncertainties, width = 1600, height = 1200, res = 180)
SparseNUTS::plot_uncertainties(snuts_fit, plot = TRUE)
dev.off()

figures$Q <- file.path(fig_dir, paste0(fig_prefix, "_Q.png"))
png(figures$Q, width = 1600, height = 1200, res = 180)
try(SparseNUTS::plot_Q(snuts_fit, Q = solve(snuts_fit$mle$Qinv)), silent = TRUE)
dev.off()

out$selected_parameters <- list(
  pairs = pair_pars,
  marginals = marginal_pars,
  trace = trace_pars
)
out$figures <- figures
saveRDS(out, out_path)

message("Wrote: ", out_path)
message("Wrote figures to: ", fig_dir)
