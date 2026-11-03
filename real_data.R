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
    input_file = "daily_returns.csv",
    output_dir = "real_data_results",
    support = NULL,
    support_padding = 0.05,
    support_calibration_trim = 30L,
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
    m_power = 0.65,
    lrc_bw_power = 0.33,
    confidence_level = 0.95,
    min_observations = 10L,
    fractional_center = "mean",
    trim = 0L,
    test_alpha = 0.05,
    test_basis_dimension = 26L,
    test_bandwidth_power = 1 / 5,
    test_bandwidth = NULL,
    test_replications = 200000L,
    test_grid_size = 1000L,
    test_batch_size = 1000L,
    test_seed = 137220731L
  )
}

validate_real_config <- function(config) {
  for (key in c("input_file", "output_dir")) {
    x <- config[[key]]
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) stop("Invalid ", key, ".")
  }
  if (!file.exists(config$input_file)) stop("Cannot find daily-return data: ", config$input_file)
  choices <- list(kernel = c("gaussian", "student", "epanechnikov"),
    bandwidth_mode = c("daily", "common"), memory_mode = c("stationary", "adaptive"),
    outside_support = c("error", "truncate"), fractional_center = c("none", "mean"))
  for (key in names(choices)) {
    x <- config[[key]]
    if (length(x) != 1L || is.na(x) || !(x %in% choices[[key]])) stop("Invalid ", key, ".")
  }
  integer_minimum <- c(grid_size = 21L, support_calibration_trim = 0L,
    min_observations = 2L, trim = 0L, test_basis_dimension = 2L,
    test_replications = 1000L, test_grid_size = 50L, test_batch_size = 1L, test_seed = 1L)
  for (key in names(integer_minimum)) {
    x <- config[[key]]
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
        x < integer_minimum[[key]] || x > .Machine$integer.max || x != floor(x)) stop("Invalid ", key, ".")
  }
  if (config$test_basis_dimension %% 2L != 0L) stop("test_basis_dimension must be even.")
  for (key in c("bw_adjust", "tail_bw_multiplier", "eps_power", "density_floor", "kde_t_df")) {
    x <- config[[key]]
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0) stop("Invalid ", key, ".")
  }
  if (config$tail_bw_multiplier < 1 || config$kde_t_df <= 2 || config$density_floor >= 0.01) {
    stop("Invalid smoothing or density-floor parameters.")
  }
  for (key in c("eps_const", "support_padding")) {
    x <- config[[key]]
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 0) stop("Invalid ", key, ".")
  }
  for (key in c("m_power", "lrc_bw_power", "confidence_level", "test_alpha")) {
    x <- config[[key]]
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0 || x >= 1) stop("Invalid ", key, ".")
  }
  if (length(config$test_bandwidth_power) != 1L || !is.finite(config$test_bandwidth_power) ||
      config$test_bandwidth_power <= 0 || config$test_bandwidth_power >= 0.5) {
    stop("test_bandwidth_power must be in (0, 0.5).")
  }
  if (!is.null(config$test_bandwidth) && (length(config$test_bandwidth) != 1L ||
      !is.numeric(config$test_bandwidth) || !is.finite(config$test_bandwidth) ||
      config$test_bandwidth < 1L || config$test_bandwidth != floor(config$test_bandwidth))) {
    stop("Invalid test_bandwidth.")
  }
  bounds <- config$tail_bounds
  if (length(bounds) != 2L || any(!is.finite(bounds)) || bounds[1L] < 0 || bounds[2L] > 1 ||
      bounds[1L] >= bounds[2L] || length(config$tail_transition) != 1L ||
      !is.finite(config$tail_transition) || config$tail_transition < 0 ||
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
  if (length(samples) < 30L) stop("At least 30 daily distributions are required.")
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
  list(samples = unit_samples, days = days, support = support, support_width = width)
}

estimate_observed_densities <- function(samples, setup, config) {
  n <- length(samples)
  bandwidths <- vapply(samples, silverman_bw, numeric(1), grid = setup$grid)
  if (config$bandwidth_mode == "common") bandwidths[] <- median(bandwidths)
  if (!is.null(config$bandwidth_return)) bandwidths[] <- config$bandwidth_return / setup$support_width
  bandwidths <- config$bw_adjust * bandwidths
  tail_blend <- tail_score(setup$grid, config$tail_bounds, config$tail_transition)
  R0 <- R1 <- matrix(0, length(setup$grid), n, dimnames = list(NULL, names(samples)))
  for (tt in seq_len(n)) {
    x <- samples[[tt]]
    h <- bandwidths[tt]
    f_base <- rowMeans(reflected_kernel_matrix(x, setup$grid, h,
      config$kernel, config$kde_t_df))
    f_tail <- rowMeans(reflected_kernel_matrix(x, setup$grid,
      config$tail_bw_multiplier * h, config$kernel, config$kde_t_df))
    R0[, tt] <- normalize_density(f_base, setup$qweights)
    R1[, tt] <- normalize_density((1 - tail_blend) * f_base + tail_blend * f_tail,
      setup$qweights)
    if (tt %% 100L == 0L || tt == n) message("KDE: ", tt, "/", n, " days.")
  }
  epsilon <- min(0.25, config$eps_const * n^(-config$eps_power))
  mixture <- function(f) normalize_density(
    sweep((1 - epsilon) * f, 1L, epsilon * setup$lambda, "+"), setup$qweights)
  list(R0 = R0, R1 = R1, R2 = mixture(R0), R3 = mixture(R1))
}

prepare_memory_data <- function(config) {
  raw <- read.csv(config$input_file, stringsAsFactors = FALSE)
  if (!all(c("day", "asset_name", "log_return") %in% names(raw)) ||
      !is.numeric(raw$log_return) || any(!is.finite(raw$log_return))) stop("Invalid daily-return input.")
  days <- sort(unique(raw$day))
  n <- length(days)
  support <- config$support
  if (is.null(support)) {
    calibration_end <- n - config$support_calibration_trim
    if (calibration_end < 30L) stop("At least 30 days are required for support calibration.")
    values <- raw$log_return[raw$day %in% days[seq_len(calibration_end)]]
    support <- range(values)
    if (diff(support) <= 0) stop("The calibration returns have zero range.")
    support <- support + c(-1, 1) * config$support_padding * diff(support)
  }
  config$support <- support
  prepared <- prepare_return_samples(config$input_file, config)
  grid <- seq(0, 1, length.out = config$grid_size)
  qweights <- trap_weights(grid)
  setup <- list(grid = grid, qweights = qweights, lambda = rep(1, length(grid)),
    pweights = qweights, support_width = prepared$support_width)
  estimator <- list(setting = if (config$memory_mode == "adaptive") "nonstationary" else "current",
    m_power = config$m_power, density_floor = config$density_floor,
    lw_interval = if (config$memory_mode == "adaptive") c(-0.49, 2.49) else c(-0.45, 0.49),
    diff_tuning_const = 0.4, max_diff_order = 3L)
  densities <- estimate_observed_densities(prepared$samples, setup, config)
  list(densities = densities, setup = setup, estimator = estimator, n = n)
}

estimate_memory <- function(data, config) {
  clr <- lapply(data$densities, clr_for_fit, config = data$estimator, setup = data$setup)
  rows <- lapply(names(clr), function(method) {
    Y <- clr[[method]]
    centered <- sweep(Y, 1L, rowMeans(Y), "-")
    Z <- sweep(centered, 1L, sqrt(data$setup$pweights), "*")
    lrc <- bartlett_long_run_covariance(Z, bw_power = config$lrc_bw_power)
    eigenfit <- eigen(lrc$Cmat, symmetric = TRUE)
    if (!is.finite(eigenfit$values[1L]) || eigenfit$values[1L] <= 0) stop("Degenerate leading direction in ", method, ".")
    V <- sweep(eigenfit$vectors[, 1L, drop = FALSE], 1L, sqrt(data$setup$pweights), "/")
    V[, 1L] <- V[, 1L] / weighted_norm(V[, 1L], data$setup$pweights)
    score <- as.numeric(crossprod(sweep(V, 1L, data$setup$pweights, "*"), centered))
    fit <- memory_score_fit(score, data$estimator)
    interval <- memory_confidence_interval(fit, config$confidence_level)
    fractional_fit <- fit$d_hat - fit$diff_order
    if (min(abs(fractional_fit - data$estimator$lw_interval)) < 0.015) {
      warning("Memory estimate for ", method, " is near its search boundary.", call. = FALSE)
    }
    data.frame(method = method, d_hat = fit$d_hat,
      ci_lower = unname(interval["lower"]), ci_upper = unname(interval["upper"]),
      stringsAsFactors = FALSE)
  })
  list(table = do.call(rbind, rows), clr = lapply(clr, t))
}

write_memory_table <- function(table, config) {
  write.csv(table, file.path(config$output_dir, "full_sample_memory.csv"), row.names = FALSE)
  number <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.3f", x))
  rows <- vapply(seq_len(nrow(table)), function(i) {
    x <- table[i, ]
    paste0(x$method, " & ", number(x$d_hat), " & [", number(x$ci_lower),
      ", ", number(x$ci_upper), "] \\\\")
  }, character(1))
  tex <- c("\\begin{table}[!htb]", "\\centering", "\\renewcommand{\\arraystretch}{1}",
    paste0("\\caption{Long memory estimates and ", 100 * config$confidence_level,
      "\\% confidence intervals for daily distributions of foreign exchange intraday returns.}\\label{table::real data}"),
    "\\fontsize{7.3pt}{9.4pt}\\selectfont", "\\setstretch{1.0}",
    "\\setlength{\\tabcolsep}{0.5in}", "\\begin{tabular}{@{}lrr@{}}", "\\toprule",
    "Method & $\\widehat d^{R,\\mathrm{est}}$ & Confidence interval \\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, file.path(config$output_dir, "full_sample_memory.tex"))
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
  list(statistic = statistic)
}

fd_bridge_reference <- function(config) {
  parameters <- list(replications = as.integer(config$test_replications),
    grid_size = as.integer(config$test_grid_size), seed = as.integer(config$test_seed),
    algorithm = "Brownian_bridge_mean_square_grid_v1")
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
  list(values = values)
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
  data.frame(series = series, method = method, d_used = d_used,
    statistic = fit$statistic, p_value = p_value, decision = decision,
    stringsAsFactors = FALSE)
}

fd_write_test_table <- function(table, config) {
  write.csv(table, file.path(config$output_dir, "d0_tests.csv"), row.names = FALSE)
  display <- lapply(seq_len(nrow(table)), function(i) {
    x <- table[i, ]
    p <- format.pval(x$p_value, digits = 3, eps = 2 / (config$test_replications + 1))
    if (startsWith(trimws(p), "<")) p <- paste0("$", trimws(p), "$")
    decision <- switch(x$decision, "Reject: d < 0" = "Reject: $d<0$",
      "Reject: d > 0" = "Reject: $d>0$", "Do not reject d = 0" = "Do not reject $d=0$")
    c(x$series, x$method, sprintf("%.3f", x$d_used), sprintf("%.4f", x$statistic), p, decision)
  })
  tex <- c("\\begin{table}[!htb]", "\\centering",
    "\\caption{Tests of $H_0:d_1=0$ before and after fractional differencing of the clr series.}",
    "\\begin{tabular}{llrrrl}", "\\toprule",
    "Series & Method & $d_{\\mathrm{used}}$ & Test statistic & $p$-value & Decision \\\\",
    "\\midrule", vapply(display, function(x) paste0(paste(x, collapse = " & "), " \\\\") , character(1)),
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, file.path(config$output_dir, "d0_tests.tex"))
}

run_real_data_analysis <- function(config = default_real_data_config()) {
  defaults <- default_real_data_config()
  unknown <- setdiff(names(config), names(defaults))
  if (length(unknown)) stop("Unknown configuration fields: ", paste(unknown, collapse = ", "))
  config <- modifyList(defaults, config, keep.null = TRUE)
  validate_real_config(config)
  if (!requireNamespace("FCVAR", quietly = TRUE)) stop('Install FCVAR first: install.packages("FCVAR")')
  data <- prepare_memory_data(config)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Full-sample memory estimation: ", data$n, " daily distributions.")
  memory <- estimate_memory(data, config)
  write_memory_table(memory$table, config)
  print(memory$table, row.names = FALSE)
  methods <- names(memory$clr)
  d_hats <- setNames(memory$table$d_hat[match(methods, memory$table$method)], methods)
  if (any(!is.finite(d_hats))) stop("Memory estimates must be finite for all four methods.")
  differenced <- lapply(methods, function(method) {
    x <- memory$clr[[method]]
    if (config$fractional_center == "mean") x <- sweep(x, 2L, colMeans(x), "-")
    z <- FCVAR::FracDiff(x = x, d = d_hats[[method]])
    if (!identical(dim(z), dim(x)) || any(!is.finite(z))) stop("Invalid FracDiff output for ", method)
    dimnames(z) <- dimnames(x)
    z
  })
  names(differenced) <- methods
  if (data$n - config$trim < 30L) stop("At least 30 curves must remain after trim.")
  keep <- seq.int(config$trim + 1L, data$n)
  reference <- fd_bridge_reference(config)
  rows <- list()
  for (series in c("Original", "Differenced")) {
    curves <- if (series == "Original") memory$clr else differenced
    for (method in methods) {
      fit <- fd_d0_statistic(curves[[method]][keep, , drop = FALSE], data$setup$grid,
        config$test_basis_dimension, config$test_bandwidth_power, config$test_bandwidth)
      rows[[length(rows) + 1L]] <- fd_test_row(fit, reference, series, method,
        if (series == "Original") 0 else d_hats[[method]], config)
    }
  }
  tests <- do.call(rbind, rows)
  rownames(tests) <- NULL
  fd_write_test_table(tests, config)
  print(tests, row.names = FALSE)
  message("Completed: ", normalizePath(config$output_dir))
  invisible(list(memory_table = memory$table, tests = tests))
}

real_data_command_line <- function(args = commandArgs(trailingOnly = TRUE)) {
  if ("--help" %in% args) {
    cat('Rscript "real_data(9).R" [--input=daily_returns.csv] [--output-dir=real_data_results]\n',
      "  [--grid-size=401] [--m-power=0.65] [--lrc-bw-power=0.33] [--confidence-level=0.95]\n",
      "  [--memory-mode=stationary|adaptive] [--support=LOWER,UPPER]\n",
      "  [--support-padding=0.05] [--support-calibration-trim=30] [--outside-support=error|truncate]\n",
      "  [--kernel=gaussian|student|epanechnikov] [--bandwidth-mode=daily|common]\n",
      "  [--bandwidth-return=VALUE] [--bw-adjust=1] [--tail-bw-multiplier=2]\n",
      "  [--tail-bounds=0.20,0.80] [--tail-transition=0.05] [--eps-const=0.10] [--eps-power=0.30]\n",
      "  [--density-floor=1e-300] [--kde-t-df=3] [--min-observations=10]\n",
      "  [--fractional-center=mean|none] [--trim=0] [--test-alpha=0.05]\n",
      "  [--test-basis-dimension=26] [--test-bandwidth-power=0.2] [--test-bandwidth=VALUE]\n",
      "  [--test-replications=200000] [--test-grid-size=1000] [--test-batch-size=1000] [--test-seed=137220731]\n",
      "Outputs: full_sample_memory.csv/.tex and d0_tests.csv/.tex\n", sep = "")
    return(invisible(NULL))
  }
  config <- default_real_data_config()
  fields <- setNames(names(config), gsub("_", "-", names(config), fixed = TRUE))
  fields <- c(fields[names(fields) != "input-file"], input = "input_file")
  strings <- c("input_file", "output_dir", "memory_mode", "outside_support",
    "kernel", "bandwidth_mode", "fractional_center")
  for (arg in args) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) < 2L || !(pair[1L] %in% names(fields))) stop("Unknown argument: ", arg)
    field <- unname(fields[pair[1L]])
    value <- paste(pair[-1L], collapse = "=")
    config[[field]] <- if (field %in% strings) value else as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
  }
  run_real_data_analysis(config)
}

if (sys.nframe() == 0L) real_data_command_line()
