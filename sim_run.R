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

make_orthonormal_basis <- function(grid, lambda, qweights, J = 5L) {
  pweights <- qweights * lambda
  basis_list <- list()
  candidate_id <- 1L
  while (length(basis_list) < J) {
    harmonic <- ceiling(candidate_id / 2)
    if (candidate_id %% 2L == 1L) {
      z <- cos(2 * pi * harmonic * grid)
    } else {
      z <- sin(2 * pi * harmonic * grid)
    }
    z <- z - sum(pweights * z)
    for (b in basis_list) {
      z <- z - weighted_inner(z, b, pweights) * b
    }
    nz <- weighted_norm(z, pweights)
    if (is.finite(nz) && nz > 1e-10) {
      basis_list[[length(basis_list) + 1L]] <- z / nz
    }
    candidate_id <- candidate_id + 1L
    if (candidate_id > 100L) stop("Could not construct enough basis functions.")
  }
  V <- do.call(cbind, basis_list)
  colnames(V) <- paste0("v", seq_len(J))
  V
}

fracint_weights <- function(d, n) {
  psi <- numeric(n)
  psi[1L] <- 1
  if (n >= 2L) {
    for (k in 2L:n) {
      psi[k] <- psi[k - 1L] * ((k - 2 + d) / (k - 1))
    }
  }
  psi
}

fft_causal_convolution <- function(x, h) {
  n_out <- length(x) + length(h) - 1L
  n_fft <- 2^(ceiling(log2(n_out)))
  fx <- fft(c(x, rep(0, n_fft - length(x))))
  fh <- fft(c(h, rep(0, n_fft - length(h))))
  z <- Re(fft(fx * fh, inverse = TRUE)) / n_fft
  z[seq_len(length(x))]
}

sim_arfima_1d1 <- function(n, d, phi = 0, theta = 0, sigma = 1, burn = 800L) {
  if (abs(phi) >= 1) stop("AR parameter phi must satisfy |phi| < 1.")
  total <- as.integer(n + burn)

  eps <- rnorm(total + 1L, sd = sigma)
  u <- eps[2L:(total + 1L)] + theta * eps[1L:total]

  psi <- fracint_weights(d, total)
  y <- fft_causal_convolution(u, psi)

  xi <- numeric(total)
  xi[1L] <- y[1L]
  if (total >= 2L) {
    for (tt in 2L:total) xi[tt] <- phi * xi[tt - 1L] + y[tt]
  }

  xi[(burn + 1L):total]
}

simulate_score_matrix <- function(T, d_vec, phi_vec, theta_vec, sigma_vec, burn = 800L) {
  J <- length(d_vec)
  phi_vec <- rep_len(phi_vec, J)
  theta_vec <- rep_len(theta_vec, J)
  sigma_vec <- rep_len(sigma_vec, J)

  Xi <- matrix(NA_real_, nrow = T, ncol = J)
  for (j in seq_len(J)) {
    Xi[, j] <- sim_arfima_1d1(
      n = T, d = d_vec[j], phi = phi_vec[j], theta = theta_vec[j],
      sigma = sigma_vec[j], burn = burn
    )
  }
  colnames(Xi) <- paste0("xi", seq_len(J))
  Xi
}

sample_from_grid_density <- function(f, grid, N) {
  p <- length(grid)
  if (length(f) != p) stop("f and grid must have the same length.")
  if (any(!is.finite(f)) || any(f < 0)) {
    stop("Density grid values must be finite and non-negative.")
  }

  h <- diff(grid)
  f0 <- f[1L:(p - 1L)]
  f1 <- f[2L:p]
  cell_mass <- 0.5 * (f0 + f1) * h
  total_mass <- sum(cell_mass)
  if (!is.finite(total_mass) || total_mass <= 0) {
    stop("Density has non-positive mass.")
  }

  cdf <- c(0, cumsum(cell_mass / total_mass))
  u <- runif(N)
  cell <- findInterval(u, cdf, all.inside = TRUE)
  cell <- pmin(cell, p - 1L)

  left_cdf <- cdf[cell]
  target_mass <- (u - left_cdf) * total_mass
  a <- f0[cell]
  b <- f1[cell] - f0[cell]
  hh <- h[cell]
  z <- target_mass / hh

  disc <- pmax(a^2 + 2 * b * z, 0)
  root_disc <- sqrt(disc)
  s <- ifelse(
    abs(b) < 1e-12,
    z / pmax(a, .Machine$double.eps),
    2 * z / pmax(a + root_disc, .Machine$double.eps)
  )
  s <- pmin(pmax(s, 0), 1)
  grid[cell] + hh * s
}

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

default_sim_config <- function() {
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(cores)) cores <- 1L
  list(
    T_values = c(200L, 500L, 1000L),
    N_values = c(100L, 200L, 400L),
    nrep = 200L,
    grid_size = 101L,
    J = 5L,
    phi_vec = c(0.25, 0.15, 0.10, 0, 0),
    theta_vec = c(-0.20, 0.10, 0, 0, 0),
    sigma_vec = c(0.34, 0.24, 0.17, 0.11, 0.08),
    bw_adjust = 1,
    density_floor = 1e-300,
    tail_bw_multiplier = 2,
    tail_bounds = c(0.20, 0.80),
    tail_transition = 0.05,
    eps_const = 0.10,
    eps_power = 0.30,
    m_power = 0.65,
    diff_tuning_const = 0.4,
    max_diff_order = 3L,
    lrc_bw_power = 0.33,
    nominal_coverage = 0.95,
    seed_base = 137220731L,
    seed_block = 10000L,
    reps_per_task = 20L,
    n_cores = max(4L, max(1L, as.integer(cores) - 1L)),
    save_rds = FALSE
  )
}

stationary_lw_experiments <- function() {
  data.frame(
    lw_experiment = c("fixed_interval", "true_d1_interval"),
    lower_intercept = c(-0.49, -0.49),
    lower_slope = c(0, 1),
    lw_upper = c(0.49, 0.49),
    interval_text = c(
      "The local Whittle search interval is $[-0.45,0.49]$. ",
      paste0("The local Whittle search interval is $(d_1-1/2,0.49]$, ",
        "where $d_1$ is the true memory parameter. ")
    ),
    stringsAsFactors = FALSE
  )
}

experiment_setup <- function(config) {
  grid <- seq(0, 1, length.out = config$grid_size)
  qweights <- trap_weights(grid)
  lambda <- rep(1, length(grid))
  V <- make_orthonormal_basis(grid, lambda, qweights, config$J)
  mu <- -config$tail_strength * (grid - 0.5)^2
  mu <- mu - sum(qweights * mu)
  list(grid = grid, qweights = qweights, lambda = lambda,
       pweights = qweights, V = V, mu = mu)
}

replication_seed <- function(T, d1, rep_id, config, N_index = 0L) {
  it <- match(T, config$T_values)
  id <- which(abs(config$d1_values - d1) < 1e-10)
  cell <- (it - 1L) * length(config$d1_values) + id - 1L
  as.integer(config$seed_base +
    (cell * (length(config$N_values) + 1L) + N_index) * config$seed_block + rep_id)
}

simulate_latent <- function(T, d1, rep_id, config, setup) {
  set.seed(replication_seed(T, d1, rep_id, config))
  d <- c(d1, 0.50 * d1, 0.25 * d1, 0, 0)
  Xi <- simulate_score_matrix(T, d, config$phi_vec, config$theta_vec,
                              config$sigma_vec, config$burn)
  Y <- sweep(setup$V %*% t(Xi), 1L, setup$mu, "+")
  list(density = inverse_clr_ref(Y, setup$lambda, setup$qweights),
       oracle_score = Xi[, 1L], scores = Xi, clr = Y)
}

reflected_kernel_matrix <- function(x, grid, h) {
  if (!is.finite(h) || h <= 0) stop("Invalid KDE bandwidth.")
  u <- outer(grid, x, "-")
  left <- outer(grid, -x, "-")
  right <- outer(grid, 2 - x, "-")
  cst <- 1 / (h * sqrt(2 * pi))
  K <- function(z) cst * exp(-0.5 * (z / h)^2)
  K(u) + K(left) + K(right)
}

estimate_densities <- function(latent_density, N, config, setup) {
  T <- ncol(latent_density)
  R0 <- R1 <- matrix(0, nrow(latent_density), T)
  tail_blend <- tail_score(setup$grid, config$tail_bounds, config$tail_transition)
  for (tt in seq_len(T)) {
    x <- sample_from_grid_density(latent_density[, tt], setup$grid, N)
    h <- config$bw_adjust * silverman_bw(x, setup$grid)
    K <- reflected_kernel_matrix(x, setup$grid, h)
    K_tail <- reflected_kernel_matrix(x, setup$grid, config$tail_bw_multiplier * h)
    f_base <- rowMeans(K)
    f_tail <- rowMeans(K_tail)
    R0[, tt] <- normalize_density(f_base, setup$qweights)
    R1[, tt] <- normalize_density((1 - tail_blend) * f_base + tail_blend * f_tail,
      setup$qweights)
  }
  epsilon <- min(0.25, config$eps_const * T^(-config$eps_power))
  mixture <- function(f) normalize_density(
    sweep((1 - epsilon) * f, 1L, epsilon * setup$lambda, "+"), setup$qweights)
  list(R0 = R0, R1 = R1, R2 = mixture(R0), R3 = mixture(R1))
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

long_run_covariance_pca_score <- function(Y, pweights, bw_power = 0.33, bandwidth = NULL) {
  Yc <- sweep(Y, 1L, rowMeans(Y), "-")
  Z <- sweep(Yc, 1L, sqrt(pweights), "*")
  lrc <- bartlett_long_run_covariance(Z, bw_power = bw_power, bandwidth = bandwidth)
  ev <- eigen(lrc$Cmat, symmetric = TRUE)
  v <- ev$vectors[, 1L] / sqrt(pweights)
  v <- v / sqrt(sum(pweights * v^2))

  score <- as.numeric(crossprod(v * pweights, Yc))
  score - mean(score)
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

stationary_memory_score_fit <- function(score, config, interval) {
  n <- length(score)
  m <- min(floor(n^config$m_power), floor((n - 1L) / 2L))
  fit <- local_whittle_fit(score, m, interval)
  fit$diff_order <- 0L
  fit$initial_d_hat <- fit$d_hat
  fit$m_used <- m
  fit
}

nonstationary_memory_score_fit <- function(score, config, interval) {
  n <- length(score)
  m <- min(floor(1 + n^config$m_power), floor((n - 1L) / 2L))
  fit <- local_whittle_fit(score, m, interval)
  initial_d_hat <- fit$d_hat
  diff_order <- 0L
  threshold <- 1 - config$diff_tuning_const / log(m)
  x <- as.numeric(score)
  while (fit$d_hat >= threshold) {
    if (diff_order >= config$max_diff_order) {
      stop("Adaptive differencing did not terminate by max_diff_order.")
    }
    x <- diff(x)
    diff_order <- diff_order + 1L
    fit <- local_whittle_fit(x, m, interval)
  }
  m_used <- min(m, floor((length(x) - 1L) / 2L))
  fit$d_hat <- fit$d_hat + diff_order
  fit$se <- 1 / (2 * sqrt(m_used))
  fit$diff_order <- diff_order
  fit$initial_d_hat <- initial_d_hat
  fit$m_used <- m_used
  fit
}

run_one_rep <- function(T, d1, rep_id, config, setup) {
  experiments <- config$lw_experiments
  experiments$lw_lower <- experiments$lower_intercept + experiments$lower_slope * d1
  fit_all_intervals <- function(score) {
    lapply(seq_len(nrow(experiments)), function(k) {
      interval <- c(experiments$lw_lower[k], experiments$lw_upper[k])
      config$fit_score(score, config, interval)
    })
  }
  latent <- simulate_latent(T, d1, rep_id, config, setup)
  oracle <- fit_all_intervals(latent$oracle_score)
  crit <- qnorm((1 + config$nominal_coverage) / 2)
  rows <- list()
  add <- function(N, method, direction, fits) {
    do.call(rbind, lapply(seq_along(fits), function(k) {
      fit <- fits[[k]]
      data.frame(T = T, N = N, d_true = d1, rep = rep_id,
        lw_experiment = experiments$lw_experiment[k],
        lw_lower = experiments$lw_lower[k], lw_upper = experiments$lw_upper[k],
        method = method, direction = direction, d_hat = fit$d_hat,
        error = fit$d_hat - d1,
        covered = abs(fit$d_hat - d1) <= crit * fit$se,
        diff_order = fit$diff_order, initial_d_hat = fit$initial_d_hat,
        lw_m = fit$m_used, observation_method = "iid", score_scale = 1,
        stringsAsFactors = FALSE)
    }))
  }
  for (N_index in seq_along(config$N_values)) {
    N <- config$N_values[N_index]
    set.seed(replication_seed(T, d1, rep_id, config, N_index))
    density <- estimate_densities(latent$density, N, config, setup)
    rows[[length(rows) + 1L]] <- add(N, "Oracle", "known_v1", oracle)
    for (method in names(density)) {
      Y <- clr_for_fit(density[[method]], config, setup)
      score_true <- as.numeric(crossprod(setup$V[, 1L] * setup$pweights, Y))
      score_estimated <- long_run_covariance_pca_score(Y, setup$pweights,
        bw_power = config$lrc_bw_power)
      rows[[length(rows) + 1L]] <- add(N, method, "known_v1",
        fit_all_intervals(score_true))
      rows[[length(rows) + 1L]] <- add(N, method, "estimated_v1",
        fit_all_intervals(score_estimated))
    }
  }
  do.call(rbind, rows)
}

run_task <- function(task, config, setup) {
  do.call(rbind, lapply(task$reps, function(r) {
    tryCatch(run_one_rep(task$T, task$d1, r, config, setup),
      error = function(e) stop("T=", task$T, ", d1=", task$d1,
        ", rep=", r, ": ", conditionMessage(e), call. = FALSE))
  }))
}

require_parallel_packages <- function() {
  missing <- c("foreach", "doSNOW")
  missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "),
      ". Install them with install.packages(c('foreach', 'doSNOW')).")
  }
  invisible(NULL)
}

run_simulation <- function(config) {
  require_parallel_packages()
  setup <- experiment_setup(config)
  cells <- expand.grid(T = config$T_values, d1 = config$d1_values)
  tasks <- list()
  for (i in seq_len(nrow(cells))) {
    for (first in seq.int(1L, config$nrep, by = config$reps_per_task)) {
      tasks[[length(tasks) + 1L]] <- list(T = cells$T[i], d1 = cells$d1[i],
        reps = first:min(first + config$reps_per_task - 1L, config$nrep))
    }
  }
  ncores <- min(config$n_cores, length(tasks))
  message("Starting doSNOW cluster with ", ncores, " worker(s); ",
    length(tasks), " tasks.")
  cl <- parallel::makePSOCKcluster(ncores)
  pb <- NULL
  on.exit({
    if (!is.null(pb)) try(close(pb), silent = TRUE)
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)
  doSNOW::registerDoSNOW(cl)
  env <- environment(run_simulation)
  functions <- Filter(function(n) is.function(get(n, envir = env)), ls(env, all.names = TRUE))
  parallel::clusterExport(cl, functions, envir = env)
  parallel::clusterExport(cl, c("config", "setup"), envir = environment())

  pb <- utils::txtProgressBar(min = 0, max = length(tasks), initial = 0, style = 3)
  progress <- function(completed) utils::setTxtProgressBar(pb, completed)

  opts <- list(progress = progress, preschedule = FALSE)
  `%dopar%` <- foreach::`%dopar%`
  pieces <- foreach::foreach(
    task = tasks,
    .options.snow = opts,
    .inorder = TRUE,
    .errorhandling = "stop"
  ) %dopar% {
    run_task(task, config, setup)
  }
  results <- do.call(rbind, pieces)
  results <- results[order(results$N, results$T, results$d_true,
                            results$rep, results$lw_experiment, results$method, results$direction), ]
  rownames(results) <- NULL
  summary <- summarise_simulation(results)
  expected <- length(config$T_values)*length(config$N_values)*length(config$d1_values)*9L*nrow(config$lw_experiments)
  if (nrow(summary) != expected || any(summary$nrep != config$nrep)) stop("Incomplete simulation grid.")
  list(config = config, results = results, summary = summary)
}

summarise_simulation <- function(results) {
  if (any(!is.finite(results$d_hat))) stop("Non-finite estimates cannot be omitted silently.")
  key <- interaction(results[c("lw_experiment", "N", "T", "d_true", "method", "direction")],
                     drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(results, key), function(g) {
    data.frame(lw_experiment = g$lw_experiment[1L],
      lw_lower = g$lw_lower[1L], lw_upper = g$lw_upper[1L],
      N = g$N[1L], T = g$T[1L], d_true = g$d_true[1L],
      method = g$method[1L], direction = g$direction[1L], nrep = nrow(g),
      bias = mean(g$error), rmse = sqrt(mean(g$error^2)),
      mae = mean(abs(g$error)), coverage = mean(g$covered),
      mean_diff_order = mean(g$diff_order),
      differenced_rate = mean(g$diff_order > 0L),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

aggregate_mc_performance <- function(summary) {
  key <- interaction(summary[c("lw_experiment", "N", "method", "direction")],
    drop = TRUE, lex.order = TRUE)
  do.call(rbind, lapply(split(summary, key), function(g) {
    data.frame(lw_experiment = g$lw_experiment[1L],
      N = g$N[1L], method = g$method[1L], direction = g$direction[1L],
      avg_abs_bias = mean(abs(g$bias)), avg_rmse = mean(g$rmse), avg_mae = mean(g$mae))
  }))
}

tex_number <- function(x, digits = 3L) {
  if (length(x) != 1L || !is.finite(x)) return("NA")
  sub("^0\\.", ".", formatC(x, digits = digits, format = "f"))
}

write_aggregate_table <- function(summary, config, path, lw_experiment) {
  a <- aggregate_mc_performance(summary)
  setting_title <- config$setting
  interval_text <- config$lw_experiments$interval_text[
    match(lw_experiment, config$lw_experiments$lw_experiment)]
  reconstruction_text <- paste0(config$baseline_text,
    "The ARFIMA scores are unscaled. Each latent density is estimated from $N$ ",
    "conditionally iid observations, sampled independently across periods, ",
    "using a reflected Gaussian KDE with a Silverman-type bandwidth. ",
    "$R_1$ blends unweighted kernel estimates at bandwidths $h$ and $")
  direction_text <- paste0("$v_1$ by Bartlett long-run covariance PCA from each method's level clr curves, with bandwidth ",
    "$\\lfloor n^{", format(config$lrc_bw_power, trim = TRUE), "}\\rfloor$. ",
    config$estimation_text, interval_text)
  values <- function(N, method, direction) {
    g <- a[a$N == N & a$method == method & a$direction == direction, ]
    if (!nrow(g)) return(rep("NA", 3L))
    vapply(g[c("avg_abs_bias", "avg_rmse", "avg_mae")], tex_number, character(1))
  }
  lines <- c("\\begin{table}[!htb]", "\\centering", "\\renewcommand{\\arraystretch}{0.8}",
    paste0("\\caption{", setting_title, ": aggregate Monte Carlo performance. ",
      "Each entry averages the cell-level ",
      "absolute bias, RMSE, or MAE equally over all reported $n$ and $d_1$ cells; each cell contains ",
      config$nrep, " replications. Panel A uses the known $v_1$, whereas Panel B estimates ",
      direction_text,
      reconstruction_text,
      format(config$tail_bw_multiplier, trim = TRUE), "h$ using a deterministic tail profile, ",
      "then normalizes the blended curve. $R_2$ applies the vanishing reference mixture to $R_0$, ",
      "and $R_3$ applies it to $R_1$. ",
      "Bias denotes average absolute bias.}"),
    paste0("\\label{tab_mc_aggregate_", gsub("--", "_", config$setting, fixed = TRUE),
      "_", lw_experiment, "}"), "\\fontsize{7.3pt}{9.4pt}\\selectfont",
    "\\setstretch{1.0}", "\\setlength{\\tabcolsep}{7pt}",
    "\\begin{tabular}{@{}clrrrrrr@{}}", "\\toprule",
    "& & \\multicolumn{3}{c}{Panel A: Known $v_1$} & \\multicolumn{3}{c}{Panel B: Estimated $v_1$} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
    "$N$ & Method & Bias & RMSE & MAE & Bias & RMSE & MAE \\\\", "\\midrule")
  Ns <- sort(unique(a$N))
  for (i in seq_along(Ns)) {
    N <- Ns[i]
    methods <- c("Oracle", "R0", "R1", "R2", "R3")
    for (j in seq_along(methods)) {
      method <- methods[j]
      label <- if (method == "Oracle") "Oracle" else paste0("$R_", substring(method, 2L), "$")
      row <- c(if (j == 1L) N else "", label,
               values(N, method, "known_v1"), values(N, method, "estimated_v1"))
      lines <- c(lines, paste0(paste(row, collapse = " & "), " \\\\"))
    }
    if (i < length(Ns)) lines <- c(lines, "\\addlinespace")
  }
  writeLines(c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"), path)
}

require_plot_packages <- function() {
  missing <- c("ggplot2", "patchwork")
  missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("The original plot style requires: ", paste(missing, collapse = ", "),
      ". Install these packages before generating reports.")
  }
  invisible(NULL)
}

performance_y_breaks <- function(limits) {

  ticks <- pretty(limits, n = 3L)
  ticks <- ticks[ticks >= limits[1L] & ticks <= limits[2L]]
  if (length(ticks) > 4L) {
    step <- ticks[2L] - ticks[1L]
    while (length(ticks) > 4L) {
      step <- 2 * step
      ticks <- seq(ceiling((limits[1L] - 1e-10) / step) * step,
                   floor((limits[2L] + 1e-10) / step) * step, by = step)
    }
  }
  ticks[abs(ticks) < 1e-12] <- 0
  ticks
}

performance_legend_guide <- function() {

  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::guide_legend(nrow = 1,
      theme = ggplot2::theme(legend.title.position = "left"))
  } else {
    ggplot2::guide_legend(nrow = 1, title.position = "left")
  }
}

plot_performance_panel <- function(summary, direction, metric = c("bias", "rmse", "coverage"),
                                   show_legend = TRUE, nominal = NULL) {
  metric <- match.arg(metric)
  methods <- c("R0", "R1", "R2", "R3")
  method_labels <- expression(R[0], R[1], R[2], R[3])
  colours <- c(R0 = "#1F77B4", R1 = "#D62728", R2 = "#2CA02C", R3 = "#9467BD")
  lty <- c(R0 = "solid", R1 = "dashed", R2 = "dotted", R3 = "dotdash")
  pch <- c(R0 = 16L, R1 = 17L, R2 = 15L, R3 = 3L)
  g <- summary[summary$direction == direction & summary$method %in% methods, ]
  if (!nrow(g)) stop("No performance rows for direction: ", direction)
  Ts <- sort(unique(g$T)); Ns <- sort(unique(g$N))
  g$method <- factor(g$method, levels = methods)
  g$N_label <- factor(paste0("N==", g$N), levels = paste0("N==", Ns))
  g$T_label <- factor(paste0("n==", g$T), levels = paste0("n==", Ts))
  g$plot_value <- g[[metric]]

  y <- summary[summary$method %in% methods, metric]

  reference <- if (metric == "bias") 0 else if (metric == "coverage") nominal else numeric(0)
  ylim <- range(c(reference, y), finite = TRUE)
  pad <- max(diff(ylim) * 0.06, 0.005)
  ylim <- ylim + c(-pad, pad)
  if (metric == "rmse") ylim[1L] <- max(0, ylim[1L])

  p <- ggplot2::ggplot(g, ggplot2::aes(
      x = d_true, y = plot_value, group = method,
      colour = method, linetype = method, shape = method)) +
    ggplot2::geom_line(linewidth = 0.35, alpha = 0.80, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.7, na.rm = TRUE) +
    ggplot2::facet_grid(N_label ~ T_label,
      labeller = ggplot2::labeller(N_label = ggplot2::label_parsed,
                                  T_label = ggplot2::label_parsed)) +
    ggplot2::scale_x_continuous(breaks = sort(unique(g$d_true)),
      labels = function(x) sprintf("%.2f", x)) +
    ggplot2::scale_y_continuous(breaks = performance_y_breaks,
      labels = function(x) sprintf("%.2f", x),
      minor_breaks = NULL, expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::scale_colour_manual(values = colours, breaks = methods, labels = method_labels) +
    ggplot2::scale_linetype_manual(values = lty, breaks = methods, labels = method_labels) +
    ggplot2::scale_shape_manual(values = pch, breaks = methods, labels = method_labels) +
    ggplot2::labs(title = NULL,
      x = if (metric == "bias") NULL else expression("True memory parameter " * d[1]),
      y = switch(metric, bias = "Bias", rmse = "RMSE", coverage = "Coverage"),
      colour = "Estimators", linetype = "Estimators", shape = "Estimators") +
    ggplot2::guides(
      colour = performance_legend_guide(),
      linetype = performance_legend_guide(),
      shape = performance_legend_guide()) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 11),
      legend.position = if (isTRUE(show_legend)) "top" else "none",
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      legend.direction = "horizontal",
      strip.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.55),
      strip.text = ggplot2::element_text(size = 12),
      panel.spacing = grid::unit(0.28, "lines"),
      panel.spacing.y = grid::unit(0.28, "lines"))
  if (metric %in% c("bias", "coverage")) {
    p <- p + ggplot2::geom_hline(yintercept = reference, colour = "black",
      linewidth = 0.4, linetype = "dashed")
  }
  p
}

make_performance_figure <- function(summary, direction) {
  require_plot_packages()

  p_bias <- plot_performance_panel(summary, direction, "bias", show_legend = TRUE)
  p_rmse <- plot_performance_panel(summary, direction, "rmse", show_legend = TRUE)
  figure <- patchwork::wrap_plots(p_bias, p_rmse, ncol = 1,
    heights = c(1, 1), guides = "collect")
  figure & ggplot2::theme(legend.position = "top")
}

make_coverage_figure <- function(summary, nominal) {
  require_plot_packages()

  p_known <- plot_performance_panel(summary, "known_v1", "coverage",
    show_legend = TRUE, nominal = nominal) +
    ggplot2::labs(x = NULL, y = expression("Coverage (known " * v[1] * ")"))
  p_estimated <- plot_performance_panel(summary, "estimated_v1", "coverage",
    show_legend = TRUE, nominal = nominal) +
    ggplot2::labs(y = expression("Coverage (estimated " * v[1] * ")"))
  figure <- patchwork::wrap_plots(p_known, p_estimated, ncol = 1,
    heights = c(1, 1), guides = "collect")
  figure & ggplot2::theme(legend.position = "top")
}

draw_performance_plot <- function(summary, direction) {
  figure <- make_performance_figure(summary, direction)
  print(figure)
  invisible(figure)
}

write_reports <- function(saved, output_dir = saved$config$output_dir) {
  config <- saved$config
  require_plot_packages()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (lw_experiment in config$lw_experiments$lw_experiment) {
    summary <- saved$summary[saved$summary$lw_experiment == lw_experiment, ]
    results <- saved$results[saved$results$lw_experiment == lw_experiment, ]
    experiment_dir <- file.path(output_dir, lw_experiment)
    dir.create(experiment_dir, recursive = TRUE, showWarnings = FALSE)
    write_aggregate_table(summary, config,
      file.path(experiment_dir, "aggregate_performance.tex"), lw_experiment)
    write.csv(results, file.path(experiment_dir, "simulation_results.csv"), row.names = FALSE)
    write.csv(summary, file.path(experiment_dir, "simulation_summary.csv"), row.names = FALSE)
    for (direction in c("known_v1", "estimated_v1")) {
      path <- file.path(experiment_dir, paste0("performance_", direction, ".png"))
      grDevices::png(path, width = 3520, height = 3200, res = 320, bg = "white")
      tryCatch(draw_performance_plot(summary, direction),
        finally = grDevices::dev.off())
    }
    path <- file.path(experiment_dir, "coverage_combined.png")
    grDevices::png(path, width = 3520, height = 3520, res = 320, bg = "white")
    tryCatch(print(make_coverage_figure(summary, config$nominal_coverage)),
      finally = grDevices::dev.off())
  }
  if (isTRUE(config$save_rds)) saveRDS(saved, file.path(output_dir, "simulation_complete.rds"))
  invisible(output_dir)
}

run_experiment <- function(config) {
  require_plot_packages()
  require_parallel_packages()
  message("Starting ", config$setting, ": ", config$nrep,
    " replications per cell; ", nrow(config$lw_experiments), " local Whittle interval experiment(s).")
  saved <- run_simulation(config)
  write_reports(saved)
  message("Tables, plots and CSV results written to ", normalizePath(config$output_dir))
  invisible(saved)
}

# stationary--unimodel ----
stationary_unimodel_config <- default_sim_config()
stationary_unimodel_config$setting <- "stationary--unimodel"
stationary_unimodel_config$d1_values <- seq(0, 0.45, by = 0.05)
stationary_unimodel_config$tail_strength <- 20
stationary_unimodel_config$burn <- 1200L
stationary_unimodel_config$fit_score <- stationary_memory_score_fit
stationary_unimodel_config$lw_experiments <- stationary_lw_experiments()
stationary_unimodel_config$baseline_text <- "The baseline density is proportional to $\\exp\\{-20(u-0.5)^2\\}$. "
stationary_unimodel_config$estimation_text <- "Local Whittle estimation uses undifferenced scores. "
stationary_unimodel_config$output_dir <- file.path("simulation_outputs", "stationary--unimodel")
stationary_unimodel_simulation <- run_experiment(stationary_unimodel_config)

# stationary--uniform ----
stationary_uniform_config <- default_sim_config()
stationary_uniform_config$setting <- "stationary--uniform"
stationary_uniform_config$d1_values <- seq(0, 0.45, by = 0.05)
stationary_uniform_config$tail_strength <- 0
stationary_uniform_config$burn <- 1200L
stationary_uniform_config$fit_score <- stationary_memory_score_fit
stationary_uniform_config$lw_experiments <- stationary_lw_experiments()
stationary_uniform_config$baseline_text <- "The baseline is $f_0(u)=1$. "
stationary_uniform_config$estimation_text <- "Local Whittle estimation uses undifferenced scores. "
stationary_uniform_config$output_dir <- file.path("simulation_outputs", "stationary--uniform")
stationary_uniform_simulation <- run_experiment(stationary_uniform_config)

# nonstationary--uniform ----
nonstationary_uniform_config <- default_sim_config()
nonstationary_uniform_config$setting <- "nonstationary--uniform"
nonstationary_uniform_config$d1_values <- c(0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90)
nonstationary_uniform_config$tail_strength <- 0
nonstationary_uniform_config$burn <- 0L
nonstationary_uniform_config$fit_score <- nonstationary_memory_score_fit
nonstationary_uniform_config$lw_experiments <- data.frame(
  lw_experiment = "fixed_interval", lower_intercept = -0.49, lower_slope = 0,
  lw_upper = 2.49, interval_text = "The local Whittle search interval is $[-0.49,2.49]$. ",
  stringsAsFactors = FALSE
)
nonstationary_uniform_config$baseline_text <- "The baseline is $f_0(u)=1$. "
nonstationary_uniform_config$estimation_text <- paste0(
  "Local Whittle estimation uses adaptive score differencing with threshold ",
  "$1-0.4/\\log(m)$, and adds the selected differencing order back to the final estimate. "
)
nonstationary_uniform_config$output_dir <- file.path("simulation_outputs", "nonstationary--uniform")
nonstationary_uniform_simulation <- run_experiment(nonstationary_uniform_config)
