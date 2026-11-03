trap_weights <- function(grid) {
  stopifnot(length(grid) >= 3L, all(diff(grid) > 0))
  h <- diff(grid)
  if (max(abs(h - h[1L])) > 1e-12) {

    w <- numeric(length(grid))
    w[1L] <- h[1L] / 2
    w[length(grid)] <- h[length(h)] / 2
    if (length(grid) > 2L) {
      w[2L:(length(grid) - 1L)] <- (h[-length(h)] + h[-1L]) / 2
    }
    return(w)
  }
  c(h[1L] / 2, rep(h[1L], length(grid) - 2L), h[1L] / 2)
}

normalize_density <- function(x, qweights) {
  if (is.null(dim(x))) {
    den <- sum(qweights * x)
    if (!is.finite(den) || den <= 0) stop("Density has non-positive integral.")
    return(x / den)
  }
  den <- colSums(sweep(x, 1L, qweights, "*"))
  if (any(!is.finite(den)) || any(den <= 0)) stop("A density column has non-positive integral.")
  sweep(x, 2L, den, "/")
}

clr_ref <- function(density, lambda, pweights, tiny = 1e-300) {
  if (is.null(dim(density))) density <- matrix(density, ncol = 1L)
  log_ratio <- sweep(log(pmax(density, tiny)), 1L, log(lambda), "-")
  spatial_mean <- colSums(sweep(log_ratio, 1L, pweights, "*"))
  out <- sweep(log_ratio, 2L, spatial_mean, "-")
  out
}

inverse_clr_ref <- function(Y, lambda, qweights) {
  if (is.null(dim(Y))) Y <- matrix(Y, ncol = 1L)
  col_max <- apply(Y, 2L, max)
  expo <- exp(sweep(Y, 2L, col_max, "-"))
  numerator <- sweep(expo, 1L, lambda, "*")
  normalize_density(numerator, qweights)
}

weighted_inner <- function(f, g, pweights) sum(pweights * f * g)

weighted_norm <- function(f, pweights) sqrt(weighted_inner(f, f, pweights))

silverman_bw <- function(x, grid) {
  sx <- sd(x)
  iq <- IQR(x) / 1.34
  candidates <- c(sx, iq)
  candidates <- candidates[is.finite(candidates) & candidates > 0]
  s <- if (length(candidates)) min(candidates) else diff(range(grid)) / 6
  h <- 0.9 * s * length(x)^(-1 / 5)
  max(h, 0.75 * min(diff(grid)))
}

tail_score <- function(x, bounds = c(0.20, 0.80), transition = 0.05) {
  a <- bounds[1L]
  b <- bounds[2L]
  if (a >= b) stop("tail bounds must satisfy lower < upper.")

  if (transition <= 0) return(as.numeric(x <= a | x >= b))

  left <- ifelse(
    x <= a, 1,
    ifelse(x >= a + transition, 0, (a + transition - x) / transition)
  )
  right <- ifelse(
    x >= b, 1,
    ifelse(x <= b - transition, 0, (x - (b - transition)) / transition)
  )
  pmax(left, right)
}

default_sim_config <- function(setting = c("current", "new", "nonstationary")) {
  setting <- match.arg(setting)
  nonstationary <- setting == "nonstationary"
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(cores)) cores <- 1L
  list(
    setting = setting,
    T_values = c(200L, 500L, 1000L),
    N_values = c(100L, 200L, 400L),
    d1_values = if (nonstationary) c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4) else seq(0, 0.45, by = 0.05),
    nrep = 200L,
    grid_size = 101L,
    J = 5L,

    tail_strength = if (setting == "current") 20 else 0,
    phi_vec = c(0.25, 0.15, 0.10, 0, 0),
    theta_vec = c(-0.20, 0.10, 0, 0, 0),
    sigma_vec = c(0.34, 0.24, 0.17, 0.11, 0.08),
    burn = if (nonstationary) 0L else 1200L,
    score_scale_method = if (nonstationary) "common_endpoint_sd" else "none",
    score_scale_horizon = 1000L,
    score_scale_target_sd = 1,
    observation_method = if (nonstationary) "fixed_bins" else "iid",
    bin_subdivisions = 10L,
    grid_bw_multiplier = 2,
    kernel = "gaussian",
    kde_t_df = 3,
    bw_adjust = 1,

    density_floor = 1e-300,

    r1_method = "tail_bandwidth_blend",
    tail_bw_multiplier = 2,
    tail_bounds = c(0.20, 0.80),
    tail_transition = 0.05,
    eps_const = 0.10,
    eps_power = 0.30,
    m_power = 0.65,
    lw_interval = if (nonstationary) c(-0.49, 2.49) else c(-0.45, 0.49),
    diff_tuning_const = 0.4,
    max_diff_order = 3L,

    direction_cov = "lrc",
    lrc_bw_power = 0.33,
    nominal_coverage = 0.95,
    seed_base = 137220731L,
    seed_block = 10000L,
    reps_per_task = 20L,
    n_cores = max(4L, max(1L, as.integer(cores) - 1L)),
    output_dir = file.path("simulation_outputs", setting),
    save_rds = FALSE
  )
}

reflected_kernel_matrix <- function(x, grid, h, kernel, df = 3) {
  if (!is.finite(h) || h <= 0) stop("Invalid KDE bandwidth.")
  u <- outer(grid, x, "-")
  left <- outer(grid, -x, "-")
  right <- outer(grid, 2 - x, "-")
  if (kernel == "epanechnikov") {
    a <- sqrt(5) * h
    K <- function(z) 0.75 * pmax(1 - (z / a)^2, 0) / a
  } else if (kernel == "gaussian") {
    cst <- 1 / (h * sqrt(2 * pi))
    K <- function(z) cst * exp(-0.5 * (z / h)^2)
  } else if (kernel == "student") {
    a <- h * sqrt((df - 2) / df)
    K <- function(z) dt(z / a, df = df) / a
  } else stop("Unknown kernel.")
  K(u) + K(left) + K(right)
}

clr_for_fit <- function(f, config, setup) {

  f <- normalize_density(pmax(f, config$density_floor * setup$lambda), setup$qweights)
  clr_ref(f, setup$lambda, setup$pweights)
}

lag_covariance_matrix <- function(Z, lag) {
  T <- ncol(Z)
  lag <- as.integer(lag)
  if (lag < 0L || lag >= T) stop("lag must satisfy 0 <= lag < T.")

  if (lag == 0L) return(tcrossprod(Z) / T)
  Z[, (lag + 1L):T, drop = FALSE] %*%
    t(Z[, 1L:(T - lag), drop = FALSE]) / T
}

bartlett_long_run_covariance <- function(Z, bw_power = 0.33, bandwidth = NULL) {

  Z <- as.matrix(Z)
  T <- ncol(Z)
  if (T < 2L || any(!is.finite(Z))) {
    stop("Z must contain finite observations at two or more time points.")
  }
  if (length(bw_power) != 1L || !is.finite(bw_power) ||
      bw_power <= 0 || bw_power >= 1) {
    stop("bw_power must lie strictly between 0 and 1.")
  }
  Z <- sweep(Z, 1L, rowMeans(Z), "-")
  b <- if (is.null(bandwidth)) max(1L, floor(T^bw_power)) else bandwidth
  if (length(b) != 1L || !is.finite(b) || b != as.integer(b) || b < 1L || b >= T) {
    stop("The Bartlett bandwidth must be an integer in [1, T-1].")
  }
  Cmat <- lag_covariance_matrix(Z, 0L)
  max_lag <- min(T - 1L, b - 1L)
  if (max_lag >= 1L) {
    for (h in seq_len(max_lag)) {
      Gh <- lag_covariance_matrix(Z, h)
      Cmat <- Cmat + (1 - h / b) * (Gh + t(Gh))
    }
  }
  list(Cmat = (Cmat + t(Cmat)) / 2, bw_used = b)
}

local_whittle_fit <- function(score, m, interval) {
  x <- score - mean(score)
  n <- length(x)
  m <- min(as.integer(m), floor((n - 1L) / 2L))
  if (m < 2L || any(!is.finite(x)) || sum(x^2) <= .Machine$double.xmin) {
    stop("Degenerate score in local Whittle fit.")
  }
  I <- pmax(Mod(fft(x)[2L:(m + 1L)])^2 / (2 * pi * n), .Machine$double.xmin)
  log_omega <- log(2 * pi * seq_len(m) / n)
  objective <- function(d) {
    z <- log(I) + 2 * d * log_omega
    zmax <- max(z)
    zmax + log(mean(exp(z - zmax))) - 2 * d * mean(log_omega)
  }
  d_hat <- optimize(objective, interval = interval)$minimum
  nu <- log(seq_len(m))
  se <- 1 / (2 * sqrt(sum((nu - mean(nu))^2)))
  if (!is.finite(d_hat)) stop("Non-finite local Whittle estimate.")
  list(d_hat = d_hat, se = se)
}

memory_score_fit <- function(score, config) {
  n <- length(score)
  nonstationary <- config$setting == "nonstationary"
  m <- if (nonstationary) floor(1 + n^config$m_power) else floor(n^config$m_power)
  m <- min(m, floor((n - 1L) / 2L))
  fit <- local_whittle_fit(score, m, config$lw_interval)
  initial_d_hat <- fit$d_hat
  diff_order <- 0L
  if (nonstationary) {
    threshold <- 1 - config$diff_tuning_const / log(m)
    x <- as.numeric(score)
    while (fit$d_hat >= threshold) {
      if (diff_order >= config$max_diff_order) {
        stop("Adaptive differencing did not terminate by max_diff_order.")
      }
      x <- diff(x)
      diff_order <- diff_order + 1L
      fit <- local_whittle_fit(x, m, config$lw_interval)
    }
    m_used <- min(m, floor((length(x) - 1L) / 2L))
    fit$d_hat <- fit$d_hat + diff_order
    fit$se <- 1 / (2 * sqrt(m_used))
  } else m_used <- m
  fit$diff_order <- diff_order
  fit$initial_d_hat <- initial_d_hat
  fit$m_used <- m_used
  fit
}

memory_confidence_interval <- function(fit, confidence_level = 0.95) {
  if (length(confidence_level) != 1L || !is.finite(confidence_level) ||
      confidence_level <= 0 || confidence_level >= 1) stop("Invalid confidence level.")
  crit <- qnorm((1 + confidence_level) / 2)
  c(lower = fit$d_hat - crit * fit$se, upper = fit$d_hat + crit * fit$se)
}

default_real_data_config <- function() {
  list(
    input_file = file.path("daily_returns.csv"),
    output_dir = "real_data_results",
    support = NULL,
    support_padding = 0.05,
    outside_support = "error",
    grid_size = 401L,
    kernel = "gaussian",
    bandwidth_mode = "daily",
    bandwidth_return = NULL,
    bw_adjust = 1,
    tail_bw_multiplier = 2,
    tail_bounds = c(0.20, 0.80),
    tail_transition = 0.05,
    eps_const = 0.10,
    eps_power = 0.30,
    density_floor = 1e-300,
    kde_t_df = 3,
    memory_mode = "stationary",
    curve_difference = 0L,
    m_powers = c(0.55, 0.65, 0.75),
    lrc_bw_power = 0.33,
    nominal_coverage = 0.95,
    min_observations = 10L,
    make_plots = TRUE
  )
}

validate_real_config <- function(config) {
  choices <- list(kernel = c("gaussian", "student", "epanechnikov"),
    bandwidth_mode = c("daily", "common"), memory_mode = c("stationary", "adaptive"),
    outside_support = c("error", "truncate"))
  for (key in names(choices)) {
    if (length(config[[key]]) != 1L || !(config[[key]] %in% choices[[key]])) stop("Invalid ", key, ".")
  }
  for (key in c("grid_size", "curve_difference", "min_observations")) {
    value <- config[[key]]
    if (length(value) != 1L || !is.finite(value) || value != as.integer(value)) stop("Invalid ", key, ".")
  }
  if (config$grid_size < 21L || config$curve_difference < 0L || config$curve_difference > 3L ||
      config$min_observations < 2L) stop("Invalid grid size, difference order, or sample size threshold.")
  if (config$memory_mode == "adaptive" && config$curve_difference != 0L) {
    stop("Use either curve differencing or the original adaptive scalar-score procedure.")
  }
  for (key in c("bw_adjust", "tail_bw_multiplier", "eps_power", "density_floor", "kde_t_df")) {
    value <- config[[key]]
    if (length(value) != 1L || !is.finite(value) || value <= 0) stop("Invalid ", key, ".")
  }
  if (config$tail_bw_multiplier < 1 || config$kde_t_df <= 2 || config$density_floor >= 0.01 ||
      length(config$eps_const) != 1L || !is.finite(config$eps_const) || config$eps_const < 0 ||
      length(config$support_padding) != 1L || !is.finite(config$support_padding) || config$support_padding < 0) {
    stop("Invalid smoothing, floor, or support parameters.")
  }
  if (!length(config$m_powers) || any(!is.finite(config$m_powers)) ||
      any(config$m_powers <= 0 | config$m_powers >= 1) || anyDuplicated(config$m_powers)) {
    stop("m_powers must contain distinct numbers in (0, 1).")
  }
  if (length(config$lrc_bw_power) != 1L || !is.finite(config$lrc_bw_power) ||
      config$lrc_bw_power <= 0 || config$lrc_bw_power >= 1 ||
      config$nominal_coverage <= 0 || config$nominal_coverage >= 1) stop("Invalid inference parameters.")
  bounds <- config$tail_bounds
  if (length(bounds) != 2L || any(!is.finite(bounds)) || bounds[1L] < 0 || bounds[2L] > 1 ||
      bounds[1L] >= bounds[2L] || !is.finite(config$tail_transition) || config$tail_transition < 0 ||
      bounds[1L] + config$tail_transition > bounds[2L] - config$tail_transition) stop("Invalid tail profile.")
  if (!is.null(config$support) && (length(config$support) != 2L || any(!is.finite(config$support)) ||
      config$support[1L] >= config$support[2L])) stop("support must be two increasing finite endpoints.")
  if (!is.null(config$bandwidth_return) && (length(config$bandwidth_return) != 1L ||
      !is.finite(config$bandwidth_return) || config$bandwidth_return <= 0)) stop("Invalid fixed bandwidth.")
  invisible(config)
}

prepare_return_samples <- function(path, config) {
  data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("day", "asset_name", "log_return")
  if (!all(required %in% names(data))) stop("Input must contain: ", paste(required, collapse = ", "))
  if (!nrow(data) || !is.numeric(data$log_return) || any(!is.finite(data$log_return))) {
    stop("Input returns must be nonempty, finite and numeric.")
  }
  if (anyNA(data$day) || anyNA(data$asset_name) || length(unique(data$asset_name)) != 1L) {
    stop("Analyze one asset at a time, with non-missing day and asset labels.")
  }
  if (anyNA(as.Date(data$day, format = "%Y-%m-%d"))) stop("Invalid day labels.")
  if ("timestamp_utc" %in% names(data)) {
    if (anyNA(data$timestamp_utc) || anyDuplicated(data$timestamp_utc)) stop("Invalid or duplicated return timestamps.")
    data <- data[order(data$day, data$timestamp_utc), , drop = FALSE]
  } else data <- data[order(data$day), , drop = FALSE]
  days <- sort(unique(data$day))
  samples <- split(data$log_return, factor(data$day, levels = days))
  if (length(samples) - config$curve_difference < 30L) stop("At least 30 curves are required after differencing.")
  before <- lengths(samples)
  support <- config$support
  automatic_support <- is.null(support)
  if (automatic_support) {
    support <- range(data$log_return)
    width <- diff(support)
    if (!is.finite(width) || width <= 0) stop("Returns have zero pooled range.")
    support <- support + c(-1, 1) * config$support_padding * width
    message("Using a common working support estimated from the pooled return range: [",
      format(support[1L], digits = 7), ", ", format(support[2L], digits = 7), "].")
  }
  outside <- vapply(samples, function(x) sum(x < support[1L] | x > support[2L]), integer(1))
  if (any(outside > 0L) && config$outside_support == "error") {
    stop(sum(outside), " returns lie outside support. Widen it or explicitly choose --outside-support=truncate.")
  }
  if (any(outside > 0L)) {
    samples <- lapply(samples, function(x) x[x >= support[1L] & x <= support[2L]])
    message("Explicitly truncating ", sum(outside), " returns; the target is the conditional distribution on support.")
  }
  if (any(lengths(samples) < config$min_observations)) {
    stop("Too few within-support returns in day(s): ", paste(days[lengths(samples) < config$min_observations], collapse = ", "))
  }
  width <- diff(support)
  unit_samples <- lapply(samples, function(x) (x - support[1L]) / width)
  list(samples = unit_samples, days = days, support = support, support_width = width,
    automatic_support = automatic_support, asset_name = unique(data$asset_name),
    diagnostics = data.frame(day = days, n_before_support = before,
      n_outside_support = outside, n_used = lengths(samples), row.names = NULL))
}

estimate_observed_densities <- function(samples, estimator_config, setup, config, functions) {
  n <- length(samples)
  bandwidths <- vapply(samples, functions$silverman_bw, numeric(1), grid = setup$grid)
  if (config$bandwidth_mode == "common") bandwidths[] <- median(bandwidths)
  if (!is.null(config$bandwidth_return)) bandwidths[] <- config$bandwidth_return / setup$support_width
  bandwidths <- config$bw_adjust * bandwidths
  tail_blend <- functions$tail_score(setup$grid, config$tail_bounds, config$tail_transition)
  R0 <- R1 <- matrix(0, length(setup$grid), n, dimnames = list(NULL, names(samples)))
  for (tt in seq_len(n)) {
    x <- samples[[tt]]
    h <- bandwidths[tt]
    f_base <- rowMeans(functions$reflected_kernel_matrix(x, setup$grid, h,
      config$kernel, config$kde_t_df))
    f_tail <- rowMeans(functions$reflected_kernel_matrix(x, setup$grid,
      config$tail_bw_multiplier * h, config$kernel, config$kde_t_df))
    R0[, tt] <- functions$normalize_density(f_base, setup$qweights)
    R1[, tt] <- functions$normalize_density((1 - tail_blend) * f_base + tail_blend * f_tail,
      setup$qweights)
    if (tt %% 100L == 0L || tt == n) message("KDE: ", tt, "/", n, " days.")
  }
  epsilon <- min(0.25, estimator_config$eps_const * n^(-estimator_config$eps_power))
  mixture <- function(f) functions$normalize_density(
    sweep((1 - epsilon) * f, 1L, epsilon * setup$lambda, "+"), setup$qweights)
  list(densities = list(R0 = R0, R1 = R1, R2 = mixture(R0), R3 = mixture(R1)),
    bandwidths = bandwidths, epsilon = epsilon, tail_blend = tail_blend)
}

default_forecast_config <- function() {
  list(input_file = file.path("daily_returns.csv"),
    output_dir = "forecast_results", holdout = 15L, cv_holdout = 15L,
    p_values = 0:2, q_values = 0:2, components = 1L, grid_size = 401L,
    memory_mode = "stationary", m_power = 0.65, lrc_bw_power = 0.33,
    confidence_level = 0.95, kernel = "gaussian", bw_adjust = 1,
    tail_bw_multiplier = 2, tail_bounds = c(0.20, 0.80), tail_transition = 0.05,
    eps_const = 0.10, eps_power = 0.30, density_floor = 1e-300, kde_t_df = 3,
    support = NULL, support_padding = 0.05, outside_support = "error",
    d0_mode = "each", make_plots = TRUE, save_forecast_densities = TRUE)
}

validate_forecast_config <- function(config) {
  for (key in c("holdout", "cv_holdout", "components", "grid_size")) {
    x <- config[[key]]
    if (length(x) != 1L || !is.finite(x) || x < 1L || x != as.integer(x)) stop("Invalid ", key, ".")
  }
  for (key in c("p_values", "q_values")) {
    x <- config[[key]]
    if (!length(x) || any(!is.finite(x)) || any(x < 0 | x != as.integer(x)) || anyDuplicated(x)) stop("Invalid ", key, ".")
  }
  if (!(config$d0_mode %in% c("each", "R0"))) stop("d0_mode must be each or R0.")
  if (!(config$memory_mode %in% c("stationary", "adaptive"))) stop("Invalid memory_mode.")
  if (!is.finite(config$confidence_level) || config$confidence_level <= 0 || config$confidence_level >= 1) stop("Invalid confidence level.")
  if (!file.exists(config$input_file)) stop("Cannot find the daily-return data: ", config$input_file)
  for (key in c("m_power", "lrc_bw_power")) {
    x <- config[[key]]
    if (length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) stop("Invalid ", key, ".")
  }
  invisible(config)
}

load_forecast_helpers <- function(config) {
  helpers <- environment(load_forecast_helpers)
  list(helpers = helpers, functions = helpers)
}

prepare_forecast_data <- function(config, helpers, functions) {
  raw <- read.csv(config$input_file, stringsAsFactors = FALSE)
  if (!all(c("day", "asset_name", "log_return") %in% names(raw)) ||
      !is.numeric(raw$log_return) || any(!is.finite(raw$log_return))) stop("Invalid daily-return input.")
  days <- sort(unique(raw$day))
  n <- length(days)
  calibration_end <- n - config$holdout - config$cv_holdout
  if (calibration_end < 30L) stop("Too few observations for the requested outer and inner rolling windows.")
  support <- config$support
  if (is.null(support)) {
    values <- raw$log_return[raw$day %in% days[seq_len(calibration_end)]]
    support <- range(values)
    if (diff(support) <= 0) stop("The calibration returns have zero range.")
    support <- support + c(-1, 1) * config$support_padding * diff(support)
  }
  basic <- helpers$default_real_data_config()
  for (key in intersect(names(config), names(basic))) basic[[key]] <- config[[key]]
  basic$support <- support
  basic$curve_difference <- 0L
  basic$bandwidth_mode <- "daily"
  helpers$validate_real_config(basic)
  prepared <- helpers$prepare_return_samples(config$input_file, basic)
  raw_samples <- split(raw$log_return, factor(raw$day, levels = prepared$days))
  grid <- seq(0, 1, length.out = config$grid_size)
  qweights <- functions$trap_weights(grid)
  setup <- list(grid = grid, qweights = qweights, lambda = rep(1, length(grid)),
    pweights = qweights, support_width = prepared$support_width)
  estimator <- functions$default_sim_config(if (config$memory_mode == "adaptive") "nonstationary" else "current")
  for (key in intersect(names(config), names(estimator))) estimator[[key]] <- config[[key]]
  estimator$nominal_coverage <- config$confidence_level
  kde <- helpers$estimate_observed_densities(prepared$samples, estimator, setup, basic, functions)
  list(prepared = prepared, raw_samples = raw_samples,
    unit_targets = lapply(raw_samples, function(x) (x - support[1L]) / diff(support)),
    base_densities = kde$densities[c("R0", "R1")], bandwidths = kde$bandwidths,
    setup = setup, estimator = estimator, calibration_end = calibration_end,
    support = support, n = n)
}

make_window_representation <- function(indices, data, config, functions) {
  setup <- data$setup
  n <- length(indices)
  epsilon <- min(0.25, config$eps_const * n^(-config$eps_power))
  base <- lapply(data$base_densities, function(x) x[, indices, drop = FALSE])
  mixture <- function(x) functions$normalize_density(
    sweep((1 - epsilon) * x, 1L, epsilon * setup$lambda, "+"), setup$qweights)
  densities <- c(base, list(R2 = mixture(base$R0), R3 = mixture(base$R1)))
  result <- list()
  for (method in names(densities)) {
    Y <- functions$clr_for_fit(densities[[method]], data$estimator, setup)
    mu <- rowMeans(Y)
    centered <- sweep(Y, 1L, mu, "-")
    Z <- sweep(centered, 1L, sqrt(setup$pweights), "*")
    lrc <- functions$bartlett_long_run_covariance(Z, bw_power = config$lrc_bw_power)
    eigenfit <- eigen(lrc$Cmat, symmetric = TRUE)
    if (!is.finite(eigenfit$values[1L]) || eigenfit$values[1L] <= 0) stop("Degenerate leading direction in ", method, ".")
    K <- config$components
    if (K >= min(n, length(setup$grid))) stop("Too many components.")
    V <- sweep(eigenfit$vectors[, seq_len(K), drop = FALSE], 1L, sqrt(setup$pweights), "/")
    for (j in seq_len(K)) V[, j] <- V[, j] / functions$weighted_norm(V[, j], setup$pweights)
    scores <- crossprod(sweep(V, 1L, setup$pweights, "*"), centered)
    fits <- lapply(seq_len(K), function(j) functions$memory_score_fit(as.numeric(scores[j, ]), data$estimator))
    result[[method]] <- list(mean = mu, basis = V, scores = scores, memory_fits = fits,
      epsilon = epsilon, lrc_bandwidth = lrc$bw_used)
  }
  result
}

prepare_ecdf_quadrature <- function(observations, grid) {
  observations <- sort(as.numeric(observations))
  if (!length(observations) || any(!is.finite(observations))) stop("Invalid future observations.")
  knots <- sort(unique(c(grid, observations)))
  left <- head(knots, -1L)
  width <- diff(knots)
  nodes <- as.vector(outer(left, rep(1, 3)) + outer(width, c(0.5 - sqrt(3 / 5) / 2, 0.5, 0.5 + sqrt(3 / 5) / 2)))
  weights <- as.vector(outer(width, c(5 / 18, 8 / 18, 5 / 18)))
  location <- findInterval(nodes, grid)
  inside <- location >= 1L & location < length(grid)
  list(nodes = nodes, weights = weights, empirical = findInterval(nodes, observations) / length(observations),
    location = location, inside = inside, grid = grid, n_observations = length(observations))
}

empirical_cdf_prediction_error <- function(density, quadrature, support_width = 1) {
  grid <- quadrature$grid
  widths <- diff(grid)
  mass <- 0.5 * (head(density, -1L) + tail(density, -1L)) * widths
  density <- density / sum(mass)
  cumulative <- c(0, cumsum(mass / sum(mass)))
  predicted <- numeric(length(quadrature$nodes))
  predicted[quadrature$location >= length(grid)] <- 1
  use <- quadrature$inside
  j <- quadrature$location[use]
  dx <- quadrature$nodes[use] - grid[j]
  predicted[use] <- cumulative[j] + density[j] * dx +
    0.5 * (density[j + 1L] - density[j]) / widths[j] * dx^2
  support_width * sum(quadrature$weights * (predicted - quadrature$empirical)^2)
}

fixed_d_arfima_forecast <- function(score, d, p, q, horizon) {
  if (!is.finite(d)) stop("Non-finite memory parameter.")
  integer_d <- max(0L, floor(d + 0.5))
  fractional_d <- d - integer_d
  if (fractional_d <= -0.5 || fractional_d >= 0.5) stop("Memory parameter is on or outside the stationary fractional boundary.")
  score_scale <- sd(score)
  if (!is.finite(score_scale) || score_scale <= 0) stop("Degenerate score process.")
  warnings <- character()
  fit <- withCallingHandlers(
    arfima::arfima(as.numeric(score) / score_scale, order = c(p, integer_d, q),
      fixed = list(frac = fractional_d), lmodel = "d", numeach = c(1, 1),
      dmean = FALSE, cpus = 1, quiet = TRUE),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  fit <- arfima::bestModes(fit, bestn = 1L)
  if (abs(fit$modes[[1L]]$dfrac - fractional_d) > 1e-10) stop("The ARFIMA engine changed the fixed d parameter.")
  prediction <- withCallingHandlers(predict(fit, n.ahead = horizon, n.use = length(score), exact = TRUE),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  value <- as.numeric(prediction[[1L]]$Forecast) * score_scale
  if (length(value) != horizon || any(!is.finite(value))) stop("Invalid ARFIMA forecasts.")
  list(mean = value, integer_d = integer_d, fractional_d = fractional_d,
    warnings = paste(unique(warnings), collapse = " | "), coefficients = coef(fit))
}

forecast_window_density <- function(representation, memory, p, q, horizon, setup, functions) {
  K <- nrow(representation$scores)
  d <- if (memory == "estimated") vapply(representation$memory_fits, `[[`, numeric(1), "d_hat") else rep(0, K)
  forecasts <- lapply(seq_len(K), function(j) fixed_d_arfima_forecast(
    representation$scores[j, ], d[j], p, q, horizon))
  score_forecast <- do.call(rbind, lapply(forecasts, `[[`, "mean"))
  Y <- sweep(representation$basis %*% score_forecast, 1L, representation$mean, "+")
  list(density = functions$inverse_clr_ref(Y, setup$lambda, setup$qweights),
    d = d, scalar_fits = forecasts)
}

write_memory_table <- function(representation, data, config, output_dir) {
  rows <- lapply(names(representation), function(method) {
    fit <- representation[[method]]$memory_fits[[1L]]
    interval <- memory_confidence_interval(fit, config$confidence_level)
    fractional_fit <- fit$d_hat - fit$diff_order
    data.frame(method = method, n = data$n, d_hat = fit$d_hat, se = fit$se,
      confidence_level = config$confidence_level, ci_lower = unname(interval["lower"]),
      ci_upper = unname(interval["upper"]), m_used = fit$m_used,
      score_difference_order = fit$diff_order, initial_d_hat = fit$initial_d_hat,
      near_optimization_boundary = min(abs(fractional_fit - data$estimator$lw_interval)) < 0.015,
      stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, rows)
  write.csv(table, file.path(output_dir, "full_sample_memory.csv"), row.names = FALSE)
  lines <- c("\\begin{table}[ht]", "\\centering",
    paste0("\\caption{Full-sample memory estimates and ", 100 * config$confidence_level, "\\% confidence intervals.}"),
    "\\begin{tabular}{lrrr}", "\\hline", "Method & $\\widehat d$ & Standard error & Confidence interval \\\\", "\\hline")
  for (i in seq_len(nrow(table))) {
    z <- table[i, ]
    lines <- c(lines, sprintf("%s & %.3f & %.3f & [%.3f, %.3f] \\\\", z$method, z$d_hat, z$se, z$ci_lower, z$ci_upper))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file.path(output_dir, "full_sample_memory.tex"))
  markdown <- c("| Method | Memory estimate | Standard error | Confidence interval |",
    "|---|---:|---:|---|", vapply(seq_len(nrow(table)), function(i) {
      z <- table[i, ]; sprintf("| %s | %.3f | %.3f | [%.3f, %.3f] |", z$method, z$d_hat, z$se, z$ci_lower, z$ci_upper)
    }, character(1)))
  writeLines(markdown, file.path(output_dir, "full_sample_memory.md"))
  table
}

forecast_cases <- function(config) {
  cases <- expand.grid(method = c("R0", "R1", "R2", "R3"), memory = c("estimated", "zero"), stringsAsFactors = FALSE)
  if (config$d0_mode == "R0") cases <- cases[cases$memory == "estimated" | cases$method == "R0", , drop = FALSE]
  cases$case <- paste(cases$method, cases$memory, sep = "_")
  cases
}

summarise_forecasts <- function(errors) {
  summarize <- function(z, horizon = FALSE) {
    result <- data.frame(method = z$method[1L], memory = z$memory[1L],
      n_forecasts = nrow(z), mean_error = mean(z$prediction_error), median_error = median(z$prediction_error),
      min_error = min(z$prediction_error), max_error = max(z$prediction_error), stringsAsFactors = FALSE)
    if (horizon) result$horizon <- z$horizon[1L]
    result
  }
  overall <- do.call(rbind, lapply(split(errors, interaction(errors$method, errors$memory, drop = TRUE)), summarize))
  by_horizon <- do.call(rbind, lapply(split(errors, interaction(errors$method, errors$memory, errors$horizon, drop = TRUE)), summarize, horizon = TRUE))
  rownames(overall) <- rownames(by_horizon) <- NULL
  by_horizon <- by_horizon[order(by_horizon$memory, by_horizon$method, by_horizon$horizon), ]
  overall <- overall[order(overall$memory, overall$mean_error), ]
  paired <- merge(errors[errors$memory == "estimated", ], errors[errors$memory == "zero", ],
    by = c("method", "split", "target_index", "horizon"), suffixes = c("_estimated", "_zero"))
  if (nrow(paired)) {
    paired$improvement <- paired$prediction_error_zero - paired$prediction_error_estimated
    paired_summary <- do.call(rbind, lapply(split(paired, paired$method), function(z) data.frame(
      method = z$method[1L], n_pairs = nrow(z), mean_error_estimated = mean(z$prediction_error_estimated),
      mean_error_zero = mean(z$prediction_error_zero), mean_improvement = mean(z$improvement),
      relative_improvement = 1 - mean(z$prediction_error_estimated) / mean(z$prediction_error_zero),
      fraction_improved = mean(z$improvement > 0), stringsAsFactors = FALSE)))
    rownames(paired_summary) <- NULL
  } else paired_summary <- data.frame()
  list(overall = overall, by_horizon = by_horizon, paired = paired_summary)
}

plot_forecast_performance <- function(summary, errors, config, n, output_dir) {
  save_plot <- function(name, plot, width, height) {
    extension <- tools::file_ext(name)
    temporary <- tempfile(pattern = "figure_", tmpdir = output_dir, fileext = paste0(".", extension))
    on.exit(unlink(temporary), add = TRUE)
    device <- "png"
    ggplot2::ggsave(temporary, plot, device = device, width = width, height = height, dpi = 320, bg = "white")
    target <- file.path(output_dir, name)
    if (file.exists(target)) unlink(target)
    if (!file.rename(temporary, target)) stop("Could not save figure: ", target)
  }
  theme <- ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.subtitle = ggplot2::element_text(color = "#586170", margin = ggplot2::margin(b = 12)),
      legend.position = "top", legend.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(15, 18, 13, 12))
  colors <- c(R0 = "#485563", R1 = "#0072B2", R2 = "#D88900", R3 = "#AD4381")
  shapes <- c(R0 = 16, R1 = 17, R2 = 15, R3 = 18)
  line_types <- c(R0 = "solid", R1 = "solid", R2 = "solid", R3 = "longdash")
  z <- summary$by_horizon
  z$memory_label <- factor(z$memory, levels = c("estimated", "zero"), labels = c("Estimated memory", "Memory fixed at zero"))
  p <- ggplot2::ggplot(z, ggplot2::aes(x = horizon, y = mean_error, color = method, shape = method, linetype = method, group = method)) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2.1) +
    ggplot2::facet_wrap(~memory_label, nrow = 1) + ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_shape_manual(values = shapes) + ggplot2::scale_linetype_manual(values = line_types) +
    ggplot2::scale_x_continuous(breaks = sort(unique(c(1, 5, 10, config$holdout)))) +
    ggplot2::scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
    ggplot2::labs(title = "Forecasting daily return distributions",
      subtitle = sprintf("%d rolling windows | %d training days | AR and MA orders chosen by nested rolling validation", config$holdout, n - config$holdout),
      x = "Forecast horizon (trading days)", y = "Mean empirical-CDF prediction error") + theme
  save_plot("prediction_performance.png", p, width = 11.4, height = 5.6)
  if (nrow(summary$paired)) {
    pair <- merge(z[z$memory == "estimated", c("method", "horizon", "mean_error")],
      z[z$memory == "zero", c("method", "horizon", "mean_error")], by = c("method", "horizon"), suffixes = c("_estimated", "_zero"))
    pair$improvement <- pair$mean_error_zero - pair$mean_error_estimated
    g <- ggplot2::ggplot(pair, ggplot2::aes(horizon, improvement, color = method, shape = method, linetype = method)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.45, color = "#8A929C", linetype = "dashed") +
      ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_shape_manual(values = shapes) + ggplot2::scale_linetype_manual(values = line_types) +
      ggplot2::scale_x_continuous(breaks = sort(unique(c(1, 5, 10, config$holdout)))) +
      ggplot2::scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
      ggplot2::labs(title = "Does estimating memory improve forecasting?",
        subtitle = "Paired comparison within each density-estimation method",
        x = "Forecast horizon (trading days)", y = "Error reduction from estimating memory") + theme
    save_plot("memory_forecast_improvement.png", g, width = 10.2, height = 5.6)
  }
  p$plot_env <- globalenv()
  for (name in names(p$mapping)) {
    if (rlang::is_quosure(p$mapping[[name]])) p$mapping[[name]] <- rlang::quo_set_env(p$mapping[[name]], globalenv())
  }
  saveRDS(p, file.path(output_dir, "prediction_performance_ggplot.rds"))
  invisible(p)
}

run_forecast_analysis <- function(config = default_forecast_config()) {
  validate_forecast_config(config)
  needed <- c("arfima", "ggplot2", "scales")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
  loaded <- load_forecast_helpers(config)
  helpers <- loaded$helpers
  functions <- loaded$functions
  data <- prepare_forecast_data(config, helpers, functions)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Full-sample estimation: ", data$n, " daily distributions.")
  full <- make_window_representation(seq_len(data$n), data, config, functions)
  memory_table <- write_memory_table(full, data, config, config$output_dir)
  print(memory_table, row.names = FALSE)
  cases <- forecast_cases(config)
  orders <- expand.grid(p = config$p_values, q = config$q_values)
  orders <- orders[order(orders$p + orders$q, orders$p, orders$q), , drop = FALSE]
  rownames(orders) <- NULL
  quadrature <- lapply(data$unit_targets, prepare_ecdf_quadrature, grid = data$setup$grid)
  outer_length <- data$n - config$holdout
  inner_length <- outer_length - config$cv_holdout
  cache <- new.env(parent = emptyenv())
  cv_rows <- error_rows <- order_rows <- density_records <- list()
  unique_cv_windows <- 0L
  compute_inner_window <- function(start) {
    key <- as.character(start)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache, inherits = FALSE))
    end <- start + inner_length - 1L
    horizon <- min(config$cv_holdout, data$n - 1L - end)
    rep <- make_window_representation(start:end, data, config, functions)
    losses <- array(Inf, dim = c(nrow(cases), nrow(orders), horizon))
    issues <- character()
    for (a in seq_len(nrow(cases))) {
      for (b in seq_len(nrow(orders))) {
        forecast <- tryCatch(forecast_window_density(rep[[cases$method[a]]], cases$memory[a],
          orders$p[b], orders$q[b], horizon, data$setup, functions), error = function(e) e)
        if (inherits(forecast, "error")) {
          issues <- c(issues, paste(cases$case[a], orders$p[b], orders$q[b], conditionMessage(forecast), sep = ": "))
        } else {
          for (h in seq_len(horizon)) losses[a, b, h] <- empirical_cdf_prediction_error(
            forecast$density[, h], quadrature[[end + h]], data$setup$support_width)
        }
      }
    }
    value <- list(losses = losses, start = start, end = end, issues = issues)
    assign(key, value, cache)
    unique_cv_windows <<- unique_cv_windows + 1L
    message("Inner CV window ", unique_cv_windows, ": days ", start, "-", end, ".")
    value
  }
  for (split in seq_len(config$holdout)) {
    start <- split
    end <- start + outer_length - 1L
    horizon <- data$n - end
    candidate_sum <- matrix(0, nrow(cases), nrow(orders))
    candidate_count <- 0L
    for (fold in seq_len(config$cv_holdout)) {
      inner <- compute_inner_window(start + fold - 1L)
      h <- config$cv_holdout - fold + 1L
      stopifnot(inner$end + h == end, inner$end < end)
      candidate_sum <- candidate_sum + apply(inner$losses[, , seq_len(h), drop = FALSE], c(1, 2), sum)
      candidate_count <- candidate_count + h
    }
    cv_loss <- candidate_sum / candidate_count
    outer_rep <- make_window_representation(start:end, data, config, functions)
    for (a in seq_len(nrow(cases))) {
      for (b in seq_len(nrow(orders))) cv_rows[[length(cv_rows) + 1L]] <- data.frame(
        split = split, train_start = start, train_end = end, method = cases$method[a], memory = cases$memory[a],
        p = orders$p[b], q = orders$q[b], cv_error = cv_loss[a, b], n_cv_forecasts = candidate_count,
        max_cv_target = end, stringsAsFactors = FALSE)
      ranking <- order(cv_loss[a, ], orders$p + orders$q, orders$p, orders$q)
      chosen <- NULL
      chosen_rank <- NA_integer_
      outer_failures <- character()
      for (rank in seq_along(ranking)) {
        b <- ranking[rank]
        if (!is.finite(cv_loss[a, b])) next
        fit <- tryCatch(forecast_window_density(outer_rep[[cases$method[a]]], cases$memory[a],
          orders$p[b], orders$q[b], horizon, data$setup, functions), error = function(e) e)
        if (!inherits(fit, "error")) { chosen <- fit; chosen_rank <- rank; break }
        outer_failures <- c(outer_failures, conditionMessage(fit))
      }
      if (is.null(chosen)) stop("No valid order for split ", split, ", ", cases$case[a], ".")
      b <- ranking[chosen_rank]
      method_rep <- outer_rep[[cases$method[a]]]
      order_rows[[length(order_rows) + 1L]] <- data.frame(split = split, train_start = start, train_end = end,
        first_forecast_index = end + 1L, last_forecast_index = data$n,
        method = cases$method[a], memory = cases$memory[a], p = orders$p[b], q = orders$q[b],
        d_used = paste(format(chosen$d, digits = 12), collapse = ";"),
        epsilon = method_rep$epsilon, cv_error = cv_loss[a, b], chosen_cv_rank = chosen_rank,
        fallback_errors = paste(unique(outer_failures), collapse = " | "), stringsAsFactors = FALSE)
      for (h in seq_len(horizon)) {
        loss <- empirical_cdf_prediction_error(chosen$density[, h], quadrature[[end + h]], data$setup$support_width)
        error_rows[[length(error_rows) + 1L]] <- data.frame(split = split, train_start = start, train_end = end,
          origin_day = data$prepared$days[end], target_index = end + h, target_day = data$prepared$days[end + h],
          horizon = h, method = cases$method[a], memory = cases$memory[a], p = orders$p[b], q = orders$q[b],
          leading_d_used = chosen$d[1L], prediction_error = loss,
          n_future_points = quadrature[[end + h]]$n_observations, stringsAsFactors = FALSE)
      }
      if (config$save_forecast_densities) density_records[[paste(split, cases$case[a], sep = "_")]] <- list(
        split = split, method = cases$method[a], memory = cases$memory[a], origin_index = end,
        target_indices = seq.int(end + 1L, data$n), density_unit = chosen$density)
    }
    message("Outer window ", split, "/", config$holdout, ": days ", start, "-", end,
      "; forecast days ", end + 1L, "-", data$n, ".")
    write.csv(do.call(rbind, error_rows), file.path(config$output_dir, "prediction_errors.csv"), row.names = FALSE)
    write.csv(do.call(rbind, order_rows), file.path(config$output_dir, "selected_orders.csv"), row.names = FALSE)
    write.csv(do.call(rbind, cv_rows), file.path(config$output_dir, "order_cross_validation.csv"), row.names = FALSE)
  }
  errors <- do.call(rbind, error_rows)
  expected <- nrow(cases) * config$holdout * (config$holdout + 1L) / 2L
  stopifnot(nrow(errors) == expected, all(errors$train_end < errors$target_index), all(is.finite(errors$prediction_error)))
  summary <- summarise_forecasts(errors)
  write.csv(summary$overall, file.path(config$output_dir, "prediction_summary.csv"), row.names = FALSE)
  write.csv(summary$by_horizon, file.path(config$output_dir, "prediction_by_horizon.csv"), row.names = FALSE)
  write.csv(summary$paired, file.path(config$output_dir, "estimated_vs_zero.csv"), row.names = FALSE)
  diagnostic <- data$prepared$diagnostics
  diagnostic$bandwidth_unit <- data$bandwidths
  diagnostic$bandwidth_return <- data$bandwidths * data$setup$support_width
  write.csv(diagnostic, file.path(config$output_dir, "density_diagnostics.csv"), row.names = FALSE)
  issues <- unlist(lapply(as.list(cache), `[[`, "issues"), use.names = FALSE)
  writeLines(unique(issues), file.path(config$output_dir, "cv_fit_failures.txt"))
  result <- list(config = config, full_sample_memory = memory_table, errors = errors, summary = summary,
    selected_orders = do.call(rbind, order_rows), cross_validation = do.call(rbind, cv_rows),
    forecasts = density_records, grid_unit = data$setup$grid,
    grid_return = data$support[1L] + diff(data$support) * data$setup$grid,
    support = data$support, support_calibration_end = data$calibration_end,
    dates = data$prepared$days, n = data$n, unique_cv_windows = unique_cv_windows,
    provenance = list(data_md5 = unname(tools::md5sum(config$input_file)),
      memory_functions_source = "Unchanged functions embedded from the supplied sim_run.R",
      R_session = sessionInfo()))
  saveRDS(result, file.path(config$output_dir, "forecast_analysis.rds"))
  if (config$make_plots) plot_forecast_performance(summary, errors, config, data$n, config$output_dir)
  capture.output(sessionInfo(), file = file.path(config$output_dir, "R_session_info.txt"))
  print(summary$overall, row.names = FALSE)
  message("Analysis completed: ", normalizePath(config$output_dir))
  invisible(result)
}

forecast_command_line <- function(args = commandArgs(trailingOnly = TRUE)) {
  if ("--help" %in% args) {
    cat("Rscript real_data_analysis_complete.R --analysis=forecast [--input=preprocessed/daily_returns.csv]\n",
      "  [--output-dir=forecast_results]\n",
      "  [--holdout=15] [--cv-holdout=15] [--p-values=0,1,2] [--q-values=0,1,2]\n",
      "  [--components=1] [--grid-size=401] [--memory-mode=adaptive|stationary]\n",
      "  [--m-power=0.65] [--lrc-bw-power=0.33] [--confidence-level=0.95]\n",
      "  [--support=LOWER,UPPER] [--support-padding=0.05] [--outside-support=error|truncate]\n",
      "  [--d0-mode=each|R0] [--bw-adjust=1] [--eps-const=0.10] [--eps-power=0.30]\n",
      "  [--no-plots] [--no-save-densities]\n", sep = "")
    return(invisible(NULL))
  }
  config <- default_forecast_config()
  fields <- c("input" = "input_file",
    "output-dir" = "output_dir", "holdout" = "holdout", "cv-holdout" = "cv_holdout",
    "p-values" = "p_values", "q-values" = "q_values", "components" = "components", "grid-size" = "grid_size",
    "memory-mode" = "memory_mode", "m-power" = "m_power", "lrc-bw-power" = "lrc_bw_power",
    "confidence-level" = "confidence_level", "support" = "support", "support-padding" = "support_padding",
    "outside-support" = "outside_support", "d0-mode" = "d0_mode", "bw-adjust" = "bw_adjust",
    "eps-const" = "eps_const", "eps-power" = "eps_power")
  strings <- c("input_file", "output_dir", "memory_mode", "outside_support", "d0_mode")
  for (arg in args) {
    if (arg == "--no-plots") { config$make_plots <- FALSE; next }
    if (arg == "--no-save-densities") { config$save_forecast_densities <- FALSE; next }
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) < 2L || !(pair[1L] %in% names(fields))) stop("Unknown argument: ", arg)
    field <- unname(fields[pair[1L]])
    value <- paste(pair[-1L], collapse = "=")
    config[[field]] <- if (field %in% strings) value else as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
  }
  run_forecast_analysis(config)
}

default_additional_analysis_config <- function() {
  list(
    input_file = file.path("daily_returns.csv"),
    output_dir = "functional_diagnostics",
    memory_overrides = list(),
    nlags = 100L,
    acf_confidence = 0.95,
    acf_estimation = "MC",
    acf_seed = 20261019L,
    fractional_center = "mean",
    trim = 0L,
    test_alpha = 0.05,
    test_basis_dimension = 26L,
    test_bandwidth_power = 1 / 5,
    test_bandwidth = NULL,
    test_replications = 200000L,
    test_grid_size = 1000L,
    test_batch_size = 1000L,
    test_seed = 137220731L,
    make_plots = TRUE
  )
}

fd_with_seed <- function(seed, expr) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  force(expr)
}

fd_validate_config <- function(config) {
  for (key in c("input_file", "output_dir")) {
    if (!is.character(config[[key]]) || length(config[[key]]) != 1L ||
        is.na(config[[key]]) || !nzchar(config[[key]])) stop("Invalid ", key, ".")
  }
  if (!file.exists(config$input_file)) stop("Cannot find preprocessed returns: ", config$input_file)
  integer_minimum <- c(nlags = 1L, trim = 0L, test_basis_dimension = 2L,
    test_replications = 1000L, test_grid_size = 50L, test_batch_size = 1L,
    acf_seed = 1L, test_seed = 1L)
  for (key in names(integer_minimum)) {
    x <- config[[key]]
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
        x < integer_minimum[[key]] || x > .Machine$integer.max || x != floor(x)) {
      stop("Invalid ", key, ".")
    }
  }
  if (config$test_basis_dimension %% 2L != 0L) stop("test_basis_dimension must be even.")
  for (key in c("acf_confidence", "test_alpha")) {
    x <- config[[key]]
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0 || x >= 1) {
      stop("Invalid ", key, ".")
    }
  }
  if (length(config$test_bandwidth_power) != 1L ||
      !is.finite(config$test_bandwidth_power) || config$test_bandwidth_power <= 0 ||
      config$test_bandwidth_power >= 0.5) stop("test_bandwidth_power must be in (0, 0.5).")
  if (!is.null(config$test_bandwidth) && (length(config$test_bandwidth) != 1L ||
      !is.numeric(config$test_bandwidth) || !is.finite(config$test_bandwidth) ||
      config$test_bandwidth < 1L || config$test_bandwidth != floor(config$test_bandwidth))) {
    stop("Invalid test_bandwidth.")
  }
  if (length(config$fractional_center) != 1L ||
      !(config$fractional_center %in% c("none", "mean"))) stop("fractional_center must be none or mean.")
  if (length(config$acf_estimation) != 1L ||
      !(config$acf_estimation %in% c("MC", "Imhof"))) stop("acf_estimation must be MC or Imhof.")
  if (!is.list(config$memory_overrides) || (length(config$memory_overrides) &&
      (is.null(names(config$memory_overrides)) || any(!nzchar(names(config$memory_overrides))) ||
       anyDuplicated(names(config$memory_overrides))))) stop("memory_overrides must be a named list.")
  if (!is.logical(config$make_plots) || length(config$make_plots) != 1L ||
      is.na(config$make_plots)) stop("Invalid make_plots.")
  invisible(config)
}

fd_core_environment <- function() {
  env <- environment(default_forecast_config)
  required <- c("default_forecast_config", "validate_forecast_config", "prepare_forecast_data",
    "make_window_representation", "write_memory_table", "normalize_density", "clr_for_fit", "trap_weights")
  missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(missing)) stop("Embedded analysis functions are missing: ", paste(missing, collapse = ", "))
  env
}

fd_trapezoidal_weights <- function(grid) {
  h <- diff(grid)
  if (length(grid) < 3L || any(!is.finite(grid)) || any(h <= 0)) stop("Invalid functional grid.")
  c(h[1L] / 2, (h[-length(h)] + h[-1L]) / 2, h[length(h)] / 2)
}

fd_fourier_basis <- function(grid, dimension) {
  if (dimension >= length(grid) - 1L) stop("Increase grid_size or reduce test_basis_dimension.")
  u <- (grid - min(grid)) / diff(range(grid))
  basis <- matrix(0, length(grid), dimension)
  for (j in seq_len(dimension / 2L)) {
    basis[, 2L * j - 1L] <- sqrt(2) * sin(2 * pi * j * u)
    basis[, 2L * j] <- sqrt(2) * cos(2 * pi * j * u)
  }
  w <- fd_trapezoidal_weights(grid)
  norms <- sqrt(colSums(sweep(basis^2, 1L, w, "*")))
  basis <- sweep(basis, 2L, norms, "/")
  gram <- crossprod(basis, sweep(basis, 1L, w, "*"))
  if (max(abs(gram - diag(dimension))) > 1e-7) {
    stop("The Fourier test requires an adequately resolved, equally spaced grid.")
  }
  basis
}

fd_d0_statistic <- function(Y, grid, basis_dimension = 26L,
                            bandwidth_power = 1 / 5, bandwidth = NULL) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  if (n < 30L || ncol(Y) != length(grid) || any(!is.finite(Y))) stop("Invalid test matrix.")
  w <- fd_trapezoidal_weights(grid)
  basis <- fd_fourier_basis(grid, basis_dimension)
  coefficients <- Y %*% sweep(basis, 1L, w, "*")
  centered <- sweep(coefficients, 2L, colMeans(coefficients), "-")
  partial_sums <- apply(centered, 2L, cumsum)
  K <- crossprod(partial_sums) / n^2
  eig <- eigen(K, symmetric = TRUE)
  if (!is.finite(eig$values[1L]) || eig$values[1L] <= .Machine$double.xmin) {
    stop("Degenerate partial-sum covariance in the d = 0 test.")
  }
  direction <- eig$vectors[, 1L]
  score <- as.numeric(centered %*% direction)
  b <- if (is.null(bandwidth)) max(1L, floor(n^bandwidth_power)) else bandwidth
  if (b >= n) stop("test_bandwidth must be smaller than the retained sample size.")
  lrv <- sum(score^2) / n
  if (b > 1L) {
    for (h in seq_len(b - 1L)) {
      lrv <- lrv + 2 * (1 - h / b) *
        sum(score[seq_len(n - h)] * score[seq.int(h + 1L, n)]) / n
    }
  }
  if (!is.finite(lrv) || lrv <= .Machine$double.xmin) stop("Nonpositive projected long-run variance.")
  statistic <- sum(cumsum(score)^2) / (n^2 * lrv)
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  energy <- sum(sweep(Yc^2, 2L, w, "*"))
  list(statistic = statistic, n = n, bandwidth = b, positive_lags = b - 1L,
    long_run_variance = lrv, basis_dimension = basis_dimension,
    projection_variance_fraction = sum(centered^2) / energy,
    direction = direction, projected_score = score)
}

fd_bridge_reference <- function(config) {
  parameters <- list(replications = as.integer(config$test_replications),
    grid_size = as.integer(config$test_grid_size), seed = as.integer(config$test_seed),
    algorithm = "Brownian_bridge_mean_square_grid_v1")
  filename <- sprintf("d0_bridge_reference_B%d_G%d_seed%d.rds",
    parameters$replications, parameters$grid_size, parameters$seed)
  path <- file.path(config$output_dir, filename)
  if (file.exists(path)) {
    saved <- readRDS(path)
    if (identical(saved$parameters, parameters) &&
        length(saved$values) == parameters$replications && all(is.finite(saved$values))) {
      return(saved)
    }
  }
  message("Simulating the demeaned-test null reference once: ", parameters$replications,
    " Brownian bridges on ", parameters$grid_size, " points.")
  values <- fd_with_seed(parameters$seed, {
    result <- numeric(parameters$replications)
    m <- parameters$grid_size - 1L
    r <- seq_len(m) / m
    starts <- seq.int(1L, parameters$replications, by = config$test_batch_size)
    for (j in seq_along(starts)) {
      index <- seq.int(starts[j], min(parameters$replications,
        starts[j] + config$test_batch_size - 1L))
      increments <- matrix(rnorm(m * length(index), sd = 1 / sqrt(m)), nrow = m)
      W <- matrix(apply(increments, 2L, cumsum), nrow = m)
      bridge <- W - tcrossprod(r, W[m, ])
      result[index] <- colSums(bridge^2) / parameters$grid_size
      if (j %% 20L == 0L || j == length(starts)) {
        message("Brownian-bridge reference: ", max(index), "/", parameters$replications)
      }
    }
    sort(result)
  })
  saved <- list(parameters = parameters, values = values)
  saveRDS(saved, path)
  saved
}

fd_test_row <- function(fit, reference, series, method, d_used, config) {
  limits <- unname(quantile(reference$values,
    c(config$test_alpha / 2, 1 - config$test_alpha / 2), type = 1))
  B <- length(reference$values)
  lower_p <- (1 + sum(reference$values <= fit$statistic)) / (B + 1)
  upper_p <- (1 + sum(reference$values >= fit$statistic)) / (B + 1)
  p_value <- min(1, 2 * min(lower_p, upper_p))
  decision <- if (fit$statistic < limits[1L]) "Reject: d < 0" else if (
    fit$statistic > limits[2L]) "Reject: d > 0" else "Do not reject d = 0"
  data.frame(series = series, method = method, d_used = d_used, n = fit$n,
    statistic = fit$statistic, p_value = p_value, critical_lower = limits[1L],
    critical_upper = limits[2L], alpha = config$test_alpha, decision = decision,
    bartlett_bandwidth = fit$bandwidth, positive_lags = fit$positive_lags,
    fourier_dimension = fit$basis_dimension,
    projection_variance_fraction = fit$projection_variance_fraction,
    calibration = if (series == "Original") "Paper: unknown-mean V0" else "Nominal after estimated d",
    fractional_center = if (series == "Original") NA_character_ else config$fractional_center,
    trim = config$trim, null_replications = B, stringsAsFactors = FALSE)
}

fd_write_test_table <- function(table, output_dir, stem) {
  write.csv(table, file.path(output_dir, paste0(stem, ".csv")), row.names = FALSE)
  display <- lapply(seq_len(nrow(table)), function(i) {
    x <- table[i, ]
    c(x$series, x$method, sprintf("%.3f", x$d_used), sprintf("%.4f", x$statistic),
      format.pval(x$p_value, digits = 3, eps = 2 / (x$null_replications + 1)), x$decision)
  })
  md <- c("| Series | Method | Difference order | V0 statistic | p-value | Decision |",
    "|---|---|---:|---:|---:|---|",
    vapply(display, function(x) paste0("| ", paste(x, collapse = " | "), " |"), character(1)))
  writeLines(md, file.path(output_dir, paste0(stem, ".md")))
  tex_display <- lapply(display, function(x) {
    x[5L] <- if (startsWith(trimws(x[5L]), "<")) paste0("$", trimws(x[5L]), "$") else x[5L]
    x[6L] <- switch(x[6L], "Reject: d < 0" = "Reject: $d<0$",
      "Reject: d > 0" = "Reject: $d>0$", "Do not reject d = 0" = "Do not reject $d=0$")
    x
  })
  tex <- c("\\begin{table}[ht]", "\\centering",
    "\\caption{Two-sided tests of $H_0:d=0$ using demeaned clr curves.}",
    "\\begin{tabular}{llrrrl}", "\\hline",
    "Series & Method & $d_{\\mathrm{used}}$ & $\\widetilde V_0$ & $p$-value & Decision \\\\",
    "\\hline", vapply(tex_display, function(x) paste0(paste(x, collapse = " & "), " \\\\"), character(1)),
    "\\hline", "\\end{tabular}",
    "\\end{table}")
  writeLines(tex, file.path(output_dir, paste0(stem, ".tex")))
}

fd_plot_one <- function(acf_result, confidence, main, small = FALSE) {
  fdaACF::plot_FACF(rho = acf_result$rho, Blueline = acf_result$Blueline,
    ci = confidence, main = main, ylim = c(0, max(1, acf_result$rho,
      acf_result$Blueline) * 1.02), xlim = c(0, length(acf_result$rho) + 1),
    lwd = if (small) 1.5 else 2, cex.main = if (small) 1 else 1.1,
    xlab = "Lag", ylab = "Autocorrelation", xaxs = "i")
}

fd_save_acf_plots <- function(acfs, d_hats, config) {
  save_plot <- function(stem, width, height, draw) {
    path <- file.path(config$output_dir, paste0(stem, ".png"))
    grDevices::png(path, width = width, height = height, units = "in", res = 220)
    tryCatch(draw(), finally = grDevices::dev.off())
  }
  save_plot("acf_original", 7.2, 6, function() {
    par(mar = c(4.3, 4.5, 2.6, 1.2))
    fd_plot_one(acfs$original$R0, config$acf_confidence, "Original clr series (R0)")
  })
  panels <- function(which) {
    par(mfrow = c(2, 2), mar = c(4, 4.1, 2.7, 1), mgp = c(2.5, 0.8, 0), cex = 0.87)
    for (method in names(acfs[[which]])) {
      title <- if (which == "original") paste0("Original clr: ", method) else {
        sprintf("%s: fractional difference, d = %.3f", method, d_hats[[method]])
      }
      fd_plot_one(acfs[[which]][[method]], config$acf_confidence, title, small = TRUE)
    }
  }
  save_plot("acf_original_by_method", 10.5, 8, function() panels("original"))
  save_plot("acf_fractionally_differenced", 10.5, 8, function() panels("differenced"))
}

run_additional_analysis <- function(config = default_additional_analysis_config()) {
  config <- modifyList(default_additional_analysis_config(), config, keep.null = TRUE)
  fd_validate_config(config)
  needed <- c("fdaACF", "FCVAR")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages first: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))")
  original <- fd_core_environment()
  memory_config <- original$default_forecast_config()
  unknown <- setdiff(names(config$memory_overrides), names(memory_config))
  if (length(unknown)) stop("Unknown memory_overrides: ", paste(unknown, collapse = ", "))
  for (key in names(config$memory_overrides)) memory_config[key] <- config$memory_overrides[key]
  memory_config$input_file <- config$input_file
  memory_config$output_dir <- config$output_dir
  memory_config$make_plots <- FALSE
  original$validate_forecast_config(memory_config)
  data <- original$prepare_forecast_data(memory_config, original, original)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Full-sample memory estimation: ", data$n, " daily distributions.")
  representation <- original$make_window_representation(seq_len(data$n), data,
    memory_config, original)
  memory_table <- original$write_memory_table(representation, data, memory_config, config$output_dir)
  print(memory_table, row.names = FALSE)
  methods <- c("R0", "R1", "R2", "R3")
  d_hats <- setNames(memory_table$d_hat[match(methods, memory_table$method)], methods)
  if (any(!is.finite(d_hats))) stop("Memory estimates must be finite for all four methods.")
  if (any(memory_table$near_optimization_boundary)) {
    warning("A memory estimate is near its search boundary. Check full_sample_memory.csv ",
      "and the memory_mode setting before interpreting its differenced diagnostics.", call. = FALSE)
  }
  epsilon <- representation$R0$epsilon
  mixture <- function(x) original$normalize_density(
    sweep((1 - epsilon) * x, 1L, epsilon * data$setup$lambda, "+"), data$setup$qweights)
  densities <- c(data$base_densities, list(R2 = mixture(data$base_densities$R0),
    R3 = mixture(data$base_densities$R1)))
  X <- lapply(densities[methods], function(f) t(original$clr_for_fit(f, data$estimator, data$setup)))
  Z <- lapply(methods, function(method) {
    x <- X[[method]]
    if (config$fractional_center == "mean") x <- sweep(x, 2L, colMeans(x), "-")
    z <- FCVAR::FracDiff(x = x, d = d_hats[[method]])
    if (!identical(dim(z), dim(x)) || any(!is.finite(z))) stop("Invalid FracDiff output for ", method)
    dimnames(z) <- dimnames(x)
    z
  })
  names(Z) <- methods
  n_kept <- data$n - config$trim
  if (n_kept < 30L) stop("At least 30 curves must remain after trim.")
  keep <- seq.int(config$trim + 1L, data$n)
  X_kept <- lapply(X, function(x) x[keep, , drop = FALSE])
  Z_kept <- lapply(Z, function(x) x[keep, , drop = FALSE])
  nlags <- min(as.integer(config$nlags), n_kept - 1L)
  if (nlags < config$nlags) warning("nlags reduced to ", nlags, " for the retained sample.", call. = FALSE)
  acfs <- list(original = list(), differenced = list())
  acf_rows <- list()
  series_number <- 0L
  for (series in names(acfs)) {
    curves <- if (series == "original") X_kept else Z_kept
    for (method in methods) {
      series_number <- series_number + 1L
      message("Functional ACF: ", series, " / ", method)
      acf_result <- fd_with_seed(config$acf_seed + series_number,
        fdaACF::obtain_FACF(Y = curves[[method]], v = data$setup$grid,
          nlags = nlags, ci = config$acf_confidence,
          estimation = config$acf_estimation, figure = FALSE))
      if (length(acf_result$rho) != nlags || any(!is.finite(acf_result$rho)) ||
          length(acf_result$Blueline) != 1L || !is.finite(acf_result$Blueline)) {
        stop("Non-finite or unexpected fdaACF output: ", series, " / ", method)
      }
      acfs[[series]][[method]] <- acf_result
      acf_rows[[length(acf_rows) + 1L]] <- data.frame(series = series, method = method,
        d_used = if (series == "original") 0 else d_hats[[method]],
        n = n_kept, lag = seq_len(nlags), rho = as.numeric(acf_result$rho),
        iid_bound = as.numeric(acf_result$Blueline), confidence = config$acf_confidence)
    }
  }
  acf_table <- do.call(rbind, acf_rows)
  write.csv(acf_table, file.path(config$output_dir, "functional_acf.csv"), row.names = FALSE)
  if (config$make_plots) fd_save_acf_plots(acfs, d_hats, config)
  reference <- fd_bridge_reference(config)
  test_details <- list(original = list(), differenced = list())
  rows <- list()
  for (series in names(test_details)) {
    curves <- if (series == "original") X_kept else Z_kept
    for (method in methods) {
      fit <- fd_d0_statistic(curves[[method]], data$setup$grid,
        config$test_basis_dimension, config$test_bandwidth_power, config$test_bandwidth)
      test_details[[series]][[method]] <- fit
      rows[[length(rows) + 1L]] <- fd_test_row(fit, reference,
        if (series == "original") "Original" else "Differenced", method,
        if (series == "original") 0 else d_hats[[method]], config)
    }
  }
  all_tests <- do.call(rbind, rows)
  main_tests <- all_tests
  rownames(main_tests) <- NULL
  fd_write_test_table(main_tests, config$output_dir, "d0_tests")
  fd_write_test_table(all_tests, config$output_dir, "d0_tests_all_methods")
  result <- list(config = config, memory_config = memory_config,
    memory_table = memory_table, acf = acfs, acf_table = acf_table,
    tests = main_tests, tests_all_methods = all_tests, test_details = test_details,
    clr_original = X, clr_fractionally_differenced = Z,
    dates = data$prepared$days, retained_indices = keep, retained_dates = data$prepared$days[keep],
    grid = data$setup$grid, support = data$support, epsilon = epsilon,
    null_reference_parameters = reference$parameters,
    provenance = list(embedded_original_code_md5 = "6dacac6a968fe825978d9333a95e9ae9",
      input_md5 = unname(tools::md5sum(config$input_file)),
      test_repository_revision = "13b6f063700a23a9314732801437e6060052dc60",
      package_versions = vapply(needed, function(x) as.character(utils::packageVersion(x)), character(1)),
      R_session = sessionInfo()))
  saveRDS(result, file.path(config$output_dir, "functional_diagnostics.rds"))
  capture.output(sessionInfo(), file = file.path(config$output_dir, "R_session_info.txt"))
  print(main_tests[, c("series", "method", "d_used", "statistic", "p_value", "decision")], row.names = FALSE)
  message("Completed: ", normalizePath(config$output_dir))
  invisible(result)
}

additional_analysis_command_line <- function(args = commandArgs(trailingOnly = TRUE)) {
  if ("--help" %in% args) {
    cat("Rscript real_data_analysis_complete.R [--analysis=diagnostics]\n",
      "  --input=preprocessed/daily_returns.csv\n",
      "  [--output-dir=functional_diagnostics] [--nlags=100]\n",
      "  [--memory-mode=stationary|adaptive] [--m-power=0.65] [--grid-size=401]\n",
      "  [--support=LOWER,UPPER] [--outside-support=error|truncate]\n",
      "  [--holdout=15] [--cv-holdout=15] [--lrc-bw-power=0.33]\n",
      "  [--fractional-center=none|mean] [--trim=0] [--test-alpha=0.05]\n",
      "  [--test-basis-dimension=26] [--test-replications=200000] [--test-grid-size=1000]\n",
      "  [--no-plots]\n", sep = "")
    return(invisible(NULL))
  }
  config <- default_additional_analysis_config()
  fields <- c("input" = "input_file",
    "output-dir" = "output_dir", "nlags" = "nlags", "fractional-center" = "fractional_center",
    "trim" = "trim", "test-alpha" = "test_alpha", "test-basis-dimension" = "test_basis_dimension",
    "test-replications" = "test_replications", "test-grid-size" = "test_grid_size")
  memory_fields <- c("memory-mode" = "memory_mode", "m-power" = "m_power", "grid-size" = "grid_size",
    "support" = "support", "outside-support" = "outside_support", "holdout" = "holdout",
    "cv-holdout" = "cv_holdout", "lrc-bw-power" = "lrc_bw_power")
  string_fields <- c("input_file", "output_dir", "fractional_center",
    "memory_mode", "outside_support")
  for (arg in args) {
    if (arg == "--no-plots") { config$make_plots <- FALSE; next }
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) < 2L || !(pair[1L] %in% c(names(fields), names(memory_fields)))) stop("Unknown argument: ", arg)
    key <- pair[1L]
    field <- if (key %in% names(fields)) unname(fields[key]) else unname(memory_fields[key])
    text <- paste(pair[-1L], collapse = "=")
    value <- if (field %in% string_fields) text else as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]])
    if (key %in% names(fields)) config[[field]] <- value else config$memory_overrides[[field]] <- value
  }
  run_additional_analysis(config)
}

real_data_command_line <- function(args = commandArgs(trailingOnly = TRUE)) {
  mode_args <- args[startsWith(args, "--analysis=")]
  if (length(mode_args) > 1L) stop("Specify --analysis only once.")
  mode <- if (length(mode_args)) sub("^--analysis=", "", mode_args) else "diagnostics"
  if (!(mode %in% c("diagnostics", "forecast"))) {
    stop("--analysis must be diagnostics or forecast.")
  }
  forwarded <- args[!startsWith(args, "--analysis=")]
  if ("--help" %in% forwarded) {
    cat("Self-contained analysis: --analysis=diagnostics (default) or --analysis=forecast.\n\n")
  }
  if (mode == "forecast") forecast_command_line(forwarded)
  else additional_analysis_command_line(forwarded)
}

if (sys.nframe() == 0L) real_data_command_line()
