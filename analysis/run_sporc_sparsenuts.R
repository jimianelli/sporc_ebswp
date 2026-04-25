#!/usr/bin/env Rscript

parse_arg <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[[1]])
}

args <- commandArgs(trailingOnly = TRUE)

model_name <- parse_arg(args, "model-name", "sporc_base")

suppressPackageStartupMessages({
  library(SparseNUTS)
  library(bayesplot)
  library(posterior)
})

fit_path <- file.path("analysis", "outputs", "base.rds")
if (!file.exists(fit_path)) {
  stop("Missing SPoRC fit: ", fit_path, ". Run analysis/run_scenarios.R first.")
}

out_dir <- file.path("analysis", "outputs", "sparsenuts")
fig_dir <- file.path("analysis", "outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(fit_path)

snuts_fit <- SparseNUTS::sample_snuts(
  obj = fit,
  cores = 1,
  model_name = model_name
)

diagnostics <- SparseNUTS::check_snuts_diagnostics(snuts_fit, print = FALSE)

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

figures$pairs <- file.path(fig_dir, "sporc_sparsenuts_pairs.png")
png(figures$pairs, width = 1800, height = 1800, res = 180)
print(bayesplot::mcmc_pairs(postwarmup[, , pair_pars, drop = FALSE]))
dev.off()

figures$marginals <- file.path(fig_dir, "sporc_sparsenuts_marginals.png")
png(figures$marginals, width = 1800, height = 1400, res = 180)
print(bayesplot::mcmc_areas(postwarmup[, , marginal_pars, drop = FALSE], prob = 0.8))
dev.off()

figures$trace <- file.path(fig_dir, "sporc_sparsenuts_trace.png")
png(figures$trace, width = 1800, height = 1400, res = 180)
print(bayesplot::mcmc_trace(postwarmup[, , trace_pars, drop = FALSE]))
dev.off()

sampler_params <- do.call(
  rbind,
  lapply(seq_along(snuts_fit$sampler_params), function(chain) {
    x <- as.data.frame(snuts_fit$sampler_params[[chain]])
    x$iteration <- seq_len(nrow(x))
    x$chain <- paste0("Chain ", chain)
    x
  })
)
sampler_long <- reshape(
  sampler_params,
  varying = setdiff(names(sampler_params), c("iteration", "chain")),
  v.names = "value",
  timevar = "parameter",
  times = setdiff(names(sampler_params), c("iteration", "chain")),
  direction = "long"
)
rownames(sampler_long) <- NULL

figures$sampler_params <- file.path(fig_dir, "sporc_sparsenuts_sampler_params.png")
png(figures$sampler_params, width = 1800, height = 1400, res = 180)
print(
  ggplot2::ggplot(sampler_long, ggplot2::aes(x = iteration, y = value, color = chain)) +
    ggplot2::geom_line(linewidth = 0.35, alpha = 0.8) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = "Iteration", y = "Value", color = "Chain") +
    ggthemes::theme_few()
)
dev.off()

figures$native_marginals <- file.path(fig_dir, "sporc_sparsenuts_native_marginals.png")
png(figures$native_marginals, width = 1800, height = 1400, res = 180)
SparseNUTS::plot_marginals(
  snuts_fit,
  pars = match(marginal_pars, available_pars),
  mfrow = c(3, 4)
)
dev.off()

figures$uncertainties <- file.path(fig_dir, "sporc_sparsenuts_uncertainties.png")
png(figures$uncertainties, width = 1600, height = 1200, res = 180)
SparseNUTS::plot_uncertainties(snuts_fit, plot = TRUE)
dev.off()

figures$Q <- file.path(fig_dir, "sporc_sparsenuts_Q.png")
png(figures$Q, width = 1600, height = 1200, res = 180)
try(SparseNUTS::plot_Q(snuts_fit, Q = solve(snuts_fit$mle$Qinv)), silent = TRUE)
dev.off()

out <- list(
  snuts_fit = snuts_fit,
  diagnostics = diagnostics,
  sampler_settings = list(
    model_name = model_name,
    sampler_call = "SparseNUTS::sample_snuts(obj = fit, cores = 1, model_name = model_name)",
    defaults = "SparseNUTS package defaults for chains, samples, warmup, metric, init, and control; cores set to 1 for RTMB serial execution."
  ),
  selected_parameters = list(
    pairs = pair_pars,
    marginals = marginal_pars,
    trace = trace_pars
  ),
  figures = figures
)

out_path <- file.path(out_dir, "sporc_sparsenuts_default.rds")
saveRDS(out, out_path)

message("Wrote: ", out_path)
message("Wrote figures to: ", fig_dir)
