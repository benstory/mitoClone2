# ============================================================================
# vireo.R - Vireo implementation for mitoClone2
# ============================================================================

# Normalize matrix rows or columns to sum to 1
normalize <- function(x, axis = 2) {
  if (axis == 2) {
    row_sums <- Matrix::rowSums(x)
    row_sums[row_sums == 0] <- 1
    return(x / row_sums)
  } else {
    col_sums <- Matrix::colSums(x)
    col_sums[col_sums == 0] <- 1
    return(t(t(x) / col_sums))
  }
}

# Amplify log likelihood by subtracting row maximum
loglik.amplify <- function(x) {
  row_max <- apply(x, 1, max)
  return(x - row_max)
}

# Calculate beta entropy or KL divergence
beta.entropy <- function(shapes, prior.shapes = NULL) {
  alpha <- shapes[, 1]
  beta <- shapes[, 2]

  if (is.null(prior.shapes)) {
    # Regular entropy
    entropy <- lbeta(alpha, beta) -
      (alpha - 1) * digamma(alpha) -
      (beta - 1) * digamma(beta) +
      (alpha + beta - 2) * digamma(alpha + beta)
  } else {
    # KL divergence
    alpha_q <- prior.shapes[, 1]
    beta_q <- prior.shapes[, 2]

    cross_entropy <- lbeta(alpha_q, beta_q) -
      (alpha_q - 1) * digamma(alpha) -
      (beta_q - 1) * digamma(beta) +
      (alpha_q + beta_q - 2) * digamma(alpha + beta)

    entropy <- lbeta(alpha, beta) -
      (alpha - 1) * digamma(alpha) -
      (beta - 1) * digamma(beta) +
      (alpha + beta - 2) * digamma(alpha + beta)

    entropy <- cross_entropy - entropy
  }

  return(sum(entropy))
}

# Update theta parameters
updateTheta <- function(model, ad, dp) {
  bd <- dp - ad
  s1.gt <- ad %*% model$id.prob
  s2.gt <- bd %*% model$id.prob

  theta.s1 <- model$theta.s1.prior
  theta.s2 <- model$theta.s2.prior

  for (ig in 1:model$n.gt) {
    theta.s1[, ig] <- theta.s1[, ig] + Matrix::rowSums(s1.gt * model$gt.prob[, , ig])
    theta.s2[, ig] <- theta.s2[, ig] + Matrix::rowSums(s2.gt * model$gt.prob[, , ig])
  }

  model$beta.mu <- theta.s1 / (theta.s1 + theta.s2)
  if (!model$fix.beta.sum) {
    model$beta.sum <- theta.s1 + theta.s2
  }

  return(model)
}

# Update ID probabilities
updateIDprob <- function(model, ad, dp) {
  bd <- dp - ad
  loglik.id <- matrix(0, nrow = ncol(ad), ncol = model$n.donor)

  theta.s1 <- model$beta.mu * model$beta.sum
  theta.s2 <- (1 - model$beta.mu) * model$beta.sum

  digamma1 <- digamma(theta.s1)
  digamma2 <- digamma(theta.s2)
  digammas <- digamma(theta.s1 + theta.s2)

  for (ig in 1:model$n.gt) {
    d1_mat <- matrix(rep(digamma1[, ig], model$n.donor),
      nrow = model$n.var, ncol = model$n.donor
    )
    d2_mat <- matrix(rep(digamma2[, ig], model$n.donor),
      nrow = model$n.var, ncol = model$n.donor
    )
    ds_mat <- matrix(rep(digammas[, ig], model$n.donor),
      nrow = model$n.var, ncol = model$n.donor
    )

    s1 <- t(ad) %*% (model$gt.prob[, , ig] * d1_mat)
    s2 <- t(bd) %*% (model$gt.prob[, , ig] * d2_mat)
    ss <- t(dp) %*% (model$gt.prob[, , ig] * ds_mat)
    loglik.id <- loglik.id + (s1 + s2 - ss)
  }

  model$id.prob <- normalize(exp(loglik.amplify(loglik.id + log(model$id.prior))))

  list(model = model, loglik.id = loglik.id)
}

# Update GT probabilities
updateGTprob <- function(model, ad, dp) {
  s1.gt <- ad %*% model$id.prob
  ss.gt <- dp %*% model$id.prob
  s2.gt <- ss.gt - s1.gt

  theta.s1 <- model$beta.mu * model$beta.sum
  theta.s2 <- (1 - model$beta.mu) * model$beta.sum

  digamma1 <- digamma(theta.s1)
  digamma2 <- digamma(theta.s2)
  digammas <- digamma(theta.s1 + theta.s2)

  loglik.gt <- array(0, dim = dim(model$gt.prior))

  for (ig in 1:model$n.gt) {
    d1_vec <- digamma1[, ig]
    d2_vec <- digamma2[, ig]
    ds_vec <- digammas[, ig]

    d1_mat <- matrix(rep(d1_vec, model$n.donor), nrow = model$n.var, ncol = model$n.donor)
    d2_mat <- matrix(rep(d2_vec, model$n.donor), nrow = model$n.var, ncol = model$n.donor)
    ds_mat <- matrix(rep(ds_vec, model$n.donor), nrow = model$n.var, ncol = model$n.donor)

    tmp <- s1.gt * d1_mat + s2.gt * d2_mat - ss.gt * ds_mat
    loglik.gt[, , ig] <- as.matrix(tmp)
  }

  ## Normalize along GT dimension
  for (i in 1:model$n.var) {
    for (j in 1:model$n.donor) {
      log_prior <- log(model$gt.prior[i, j, ])
      log_prob <- loglik.gt[i, j, ] + log_prior
      log_prob <- log_prob - max(log_prob)
      prob <- exp(log_prob)
      model$gt.prob[i, j, ] <- prob / sum(prob)
    }
  }

  return(model)
}

#' Initialize Vireo model
#' @param n.cell Number of cells
#' @param n.var Number of variants
#' @param n.donor Number of donors
#' @param n.gt Number of genotype states (default 3: 0,1,2)
#' @param learn.gt Whether to learn genotype probabilities
#' @param fix.beta.sum Whether to fix beta sum parameters
#' @param beta.mu.init Initial beta mu values
#' @param beta.sum.init Initial beta sum values
#' @param id.prob.init Initial ID probabilities
#' @param gt.prob.init Initial genotype probabilities
vireo <- function(n.cell, n.var, n.donor, n.gt = 3, learn.gt = TRUE,
                  fix.beta.sum = FALSE, beta.mu.init = NULL, beta.sum.init = NULL,
                  id.prob.init = NULL, gt.prob.init = NULL) {
  theta.len <- n.var

  if (is.null(beta.mu.init)) {
    beta.mu <- matrix(rep(seq(0.01, 0.99, length.out = n.gt), theta.len),
      nrow = theta.len, byrow = TRUE
    )
  } else {
    beta.mu <- beta.mu.init
  }

  if (is.null(beta.sum.init)) {
    beta.sum <- matrix(50, nrow = theta.len, ncol = n.gt)
  } else {
    beta.sum <- beta.sum.init
  }

  if (is.null(id.prob.init)) {
    id.prob <- normalize(matrix(stats::runif(n.cell * n.donor), nrow = n.cell))
  } else {
    id.prob <- normalize(id.prob.init)
  }

  if (is.null(gt.prob.init)) {
    gt.prob <- array(stats::runif(n.var * n.donor * n.gt),
      dim = c(n.var, n.donor, n.gt)
    )
    gt_sums <- apply(gt.prob, c(1, 2), sum)
    gt_sums_expanded <- array(gt_sums, dim = c(n.var, n.donor, n.gt))
    gt.prob <- gt.prob / gt_sums_expanded
  } else {
    gt.prob <- gt.prob.init
  }

  ## set the priors
  theta.s1.prior <- matrix(rep(seq(0.01, 0.99, length.out = n.gt), theta.len) * 50,
    nrow = theta.len, byrow = TRUE
  )
  theta.s2.prior <- matrix(rep((1 - seq(0.01, 0.99, length.out = n.gt)), theta.len) * 50,
    nrow = theta.len, byrow = TRUE
  )

  id.prior <- normalize(matrix(1, nrow = n.cell, ncol = n.donor))

  gt.prior <- array(1, dim = c(n.var, n.donor, n.gt))
  for (i in 1:n.var) {
    for (j in 1:n.donor) {
      gt.prior[i, j, ] <- gt.prior[i, j, ] / sum(gt.prior[i, j, ])
    }
  }

  list(
    n.cell = n.cell,
    n.var = n.var,
    n.donor = n.donor,
    n.gt = n.gt,
    learn.gt = learn.gt,
    fix.beta.sum = fix.beta.sum,
    beta.mu = beta.mu,
    beta.sum = beta.sum,
    id.prob = id.prob,
    gt.prob = gt.prob,
    theta.s1.prior = theta.s1.prior,
    theta.s2.prior = theta.s2.prior,
    id.prior = id.prior,
    gt.prior = gt.prior,
    elbo = numeric(0)
  )
}

# Calculate ELBO
get.elbo <- function(model, loglik.id = NULL, ad = NULL, dp = NULL) {
  if (is.null(loglik.id)) {
    result <- updateIDprob(model, ad, dp)
    loglik.id <- result$loglik.id
  }

  ## Expected log-likelihood term
  lb.p <- sum(loglik.id * model$id.prob)

  ## Convert to dense matrices for vectorized KL calculations
  id_prob_dense <- as.matrix(model$id.prob)
  id_prior_dense <- as.matrix(model$id.prior)

  ## ID KL divergence - vectorized
  valid_id <- id_prob_dense > 0 & id_prior_dense > 0
  kl.id <- sum(id_prob_dense[valid_id] * log(id_prob_dense[valid_id] / id_prior_dense[valid_id]))

  ## GT KL divergence - vectorized
  valid_gt <- model$gt.prob > 0 & model$gt.prior > 0
  kl.gt <- sum(model$gt.prob[valid_gt] * log(model$gt.prob[valid_gt] / model$gt.prior[valid_gt]))

  ## Beta KL divergence
  theta.s1 <- model$beta.mu * model$beta.sum
  theta.s2 <- (1 - model$beta.mu) * model$beta.sum

  shapes <- cbind(as.vector(theta.s1), as.vector(theta.s2))
  prior.shapes <- cbind(
    as.vector(model$theta.s1.prior),
    as.vector(model$theta.s2.prior)
  )

  kl.theta <- beta.entropy(shapes, prior.shapes)

  ## Return ELBO
  elbo <- lb.p - kl.id - kl.gt - kl.theta

  if (!is.finite(elbo)) {
    return(-Inf)
  }

  return(elbo)
}

# Fit the Vireo model
fit.vireo <- function(model, ad, dp, max.iter = 200, min.iter = 5,
                      epsilon.conv = 1e-2, delay.fit.theta = 0,
                      verbose = TRUE) {
  elbo <- numeric(max.iter)

  for (it in 1:max.iter) {
    ## Update theta if learning and past delay
    if (it > delay.fit.theta) {
      model <- updateTheta(model, ad, dp)
    }

    ## Update GT if learning
    if (model$learn.gt) {
      model <- updateGTprob(model, ad, dp)
    }

    ## Update ID probabilities
    result <- updateIDprob(model, ad, dp)
    model <- result$model
    loglik.id <- result$loglik.id

    ## Calculate ELBO
    elbo[it] <- get.elbo(model, loglik.id)

    ## Check convergence
    if (it > min.iter) {
      if (elbo[it] < elbo[it - 1] - 1e-6) {
        if (verbose) {
          cat("Warning: ELBO decreases!\n")
        }
      } else if (it == max.iter) {
        if (verbose) {
          cat("Warning: Did not converge!\n")
        }
      } else if (elbo[it] - elbo[it - 1] < epsilon.conv) {
        break
      }
    }

    if (verbose && it %% 20 == 0) {
      cat("Iteration", it, "ELBO:", elbo[it], "\n")
    }
  }

  model$elbo <- elbo[1:it]
  return(model)
}

#' Fit Vireo model with multiple initializations
#'
#' @param data A mitoClone2 data object containing M (ALT) and N (non-ALT) matrices
#' @param n.donor Number of donors to identify
#' @param n.gt Number of genotype states (default 3)
#' @param learn.gt Whether to learn genotype probabilities
#' @param n.init Number of random initializations
#' @param max.iter Maximum iterations per initialization
#' @param random.seed Random seed for reproducibility
#' @param verbose Print progress messages
#' @param ... Additional arguments passed to vireo.filter
#' @return Best fitted Vireo model
#' @examples load(system.file("extdata/LudwigFig7.Rda", package = "mitoClone2"))
#' test.data <- list(N = as.matrix(t(LudwigFig7@N)), M = as.matrix(t(LudwigFig7@M)))
#' vireoModel <- vireoFit(test.data, n.donor = 9, filter.variants = FALSE, min_cells_per_sample = 5)
#' @export
vireoFit <- function(data, n.donor, n.gt = 3, learn.gt = TRUE,
                     n.init = 10, max.iter = 200,
                     random.seed = NULL, verbose = TRUE, ...) {
  if (!is.null(random.seed)) {
    set.seed(random.seed)
  }

  ## extract and prepare data
  if (!all(c("M", "N") %in% names(data))) {
    stop("Data object must contain M (ALT) and N (non-ALT) matrices")
  }

  ## apply variant filtering if requested
  ad <- data$M
  dp <- data$M + data$N

  # set counts
  n.var <- nrow(ad)
  n.cell <- ncol(ad)

  if (verbose) {
    cat("Fitting Vireo model with:\n")
    cat("  Variants:", n.var, "\n")
    cat("  Cells:", n.cell, "\n")
    cat("  Donors:", n.donor, "\n")
  }

  best.elbo <- -Inf
  best.model <- NULL

  for (i in 1:n.init) {
    if (verbose) {
      cat("Initialization", i, "of", n.init, "\n")
    }

    model <- vireo(n.cell, n.var, n.donor, n.gt, learn.gt)
    model <- fit.vireo(model, ad, dp, max.iter = max.iter, verbose = verbose)

    final.elbo <- utils::tail(model$elbo, 1)
    if (final.elbo > best.elbo) {
      best.elbo <- final.elbo
      best.model <- model
    }
  }

  if (verbose) {
    cat("Best ELBO:", best.elbo, "\n")
  }

  return(best.model)
}

#' Predict cell assignments from fitted Vireo model
#'
#' @param model Fitted Vireo model
#' @param threshold Minimum probability threshold for assignment
#' @return Data frame with cell assignments and probabilities
#' @examples
#' load(system.file("extdata/LudwigFig7.Rda", package = "mitoClone2"))
#' test.data <- list(N = as.matrix(t(LudwigFig7@N)), M = as.matrix(t(LudwigFig7@M)))
#' vireoModel <- vireoFit(test.data, n.donor = 9, filter.variants = FALSE, min_cells_per_sample = 5)
#' cellAssignments <- predictCellAssignment(vireoModel, threshold = 0.9)
#' @export
predictCellAssignment <- function(model, threshold = 0.9) {
  max.prob <- apply(model$id.prob, 1, max)
  assignments <- apply(model$id.prob, 1, which.max)

  ## mark uncertain assignments
  assignments[max.prob < threshold] <- NA

  data.frame(
    cell = 1:model$n.cell,
    donor = assignments,
    prob = max.prob,
    stringsAsFactors = FALSE
  )
}
