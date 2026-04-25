extract_sporc_timeseries <- function(model, years = NULL) {
  rep <- model$rep
  if (is.null(rep)) {
    stop("SPoRC model does not include a rep component")
  }

  n_years <- length(drop(rep$Aggregated_SSB))
  if (is.null(years)) {
    years <- seq_len(n_years)
  }
  if (length(years) != n_years) {
    stop("years has length ", length(years), "; expected ", n_years)
  }

  fmort <- drop(rep$Fmort)
  if (length(fmort) == n_years) {
    fmort <- as.numeric(fmort)
  } else {
    fmort <- rep(NA_real_, n_years)
  }

  data.frame(
    year = rep(years, 4),
    quantity = rep(c("F", "SSB", "Total_Biom", "Recruit"), each = n_years),
    value = c(
      fmort,
      as.numeric(drop(rep$Aggregated_SSB)),
      as.numeric(drop(rep$Total_Biom)),
      as.numeric(drop(rep$Rec))
    ),
    model = "sporc",
    row.names = NULL
  )
}

compare_timeseries <- function(sporc_ts, admb_ts) {
  admb_ts$model <- "admb"

  merged <- merge(
    sporc_ts[, c("year", "quantity", "value")],
    admb_ts[, c("year", "quantity", "value")],
    by = c("year", "quantity"),
    suffixes = c("_sporc", "_admb")
  )

  merged$difference <- merged$value_sporc - merged$value_admb
  merged$relative_difference <- merged$difference / merged$value_admb
  merged
}
