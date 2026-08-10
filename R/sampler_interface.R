# Rcpp dispatch, progress handling, and multi-chain combination.

.blm_gibbs_rcpp <- function(y, x, residual_shape, residual_scale,
                            iterations, burnin, thin,
                            progress_callback = NULL,
                            block_id = NULL, block_model = 0L,
                            normal_shape = 2, normal_scale = 1,
                            pi_alpha = 1, pi_beta = 1,
                            spike_var_shape = 2, spike_var_scale = 1,
                            global_scale = 1, residual_var = NULL,
                            local_a = 1, local_b = 0.5,
                            multi_gamma = list(c(0, 0.01, 0.1, 1)),
                            multi_pi_alpha = list(rep(1, 4)),
                            multi_var_shape = 2, multi_var_scale = 1,
                            store_samples = TRUE,
                            store_coefficient_cov = TRUE,
                            effective_n = NULL, fit_intercept = TRUE,
                            intercept_x_mean = NULL,
                            intercept_y_mean = NULL,
                            XtX = NULL, XtX_center = NULL,
                            XtX_indices = NULL, XtX_types = NULL, Xty = NULL,
                            ld_blocks = NULL, ld_indices = NULL,
                            ld_scale = NULL,
                            eigen_X = NULL, eigen_y = NULL,
                            eigen_indices = NULL,
                            nthreads = 1L,
                            yty = NULL, center_observations = TRUE,
                            residual_sse_offset = 0, compute_pve = FALSE,
                            pve_type = c("standalone", "allocated")) {
  pve_controls <- .validate_pve_controls(compute_pve, pve_type)
  compute_pve <- pve_controls$compute_pve
  pve_type <- pve_controls$pve_type
  .validate_mcmc(iterations, burnin, thin)
  use_sufficient_statistics <- !is.null(XtX) || !is.null(eigen_X) ||
    !is.null(ld_blocks)
  if (is.null(block_id)) {
    block_id <- rep.int(
      1L, if (use_sufficient_statistics) length(Xty) else ncol(x)
    )
  }
  if (is.null(progress_callback)) {
    progress_callback <- function(amount, iteration) invisible(NULL)
  }
  if (is.null(effective_n)) effective_n <- length(y)
  if (is.null(intercept_x_mean)) intercept_x_mean <- colMeans(x)
  if (is.null(intercept_y_mean)) intercept_y_mean <- mean(y)
  common_arguments <- list(
    y = y,
    X = x,
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    progress_callback = progress_callback,
    block_id = block_id,
    block_model = block_model,
    normal_shape = normal_shape,
    normal_scale = normal_scale,
    pi_alpha = pi_alpha,
    pi_beta = pi_beta,
    spike_var_shape = spike_var_shape,
    spike_var_scale = spike_var_scale,
    global_scale = global_scale,
    local_a = local_a,
    local_b = local_b,
    multi_gamma_list = multi_gamma,
    multi_pi_alpha_list = multi_pi_alpha,
    multi_var_shape = multi_var_shape,
    multi_var_scale = multi_var_scale,
    learn_residual_var = is.null(residual_var),
    fixed_residual_var = if (is.null(residual_var)) 1 else residual_var,
    store_samples = store_samples,
    store_coefficient_cov = store_coefficient_cov,
    effective_n = effective_n,
    fit_intercept = fit_intercept,
    intercept_x_mean = intercept_x_mean,
    intercept_y_mean = intercept_y_mean,
    use_sufficient_statistics = use_sufficient_statistics,
    summary_Xty = if (use_sufficient_statistics) Xty else numeric(),
    summary_yty = if (use_sufficient_statistics && !is.null(yty)) yty else 0,
    center_observations = center_observations,
    residual_sse_offset = residual_sse_offset,
    compute_pve = compute_pve,
    pve_type_code = match(pve_type, c("standalone", "allocated")) - 1L
  )
  samples <- if (!is.null(ld_blocks)) {
    do.call(
      blm_gibbs_ld_rcpp_cpp,
      c(
        common_arguments[c(
          "residual_shape", "residual_scale", "iterations", "burnin", "thin",
          "progress_callback", "block_id", "block_model", "normal_shape",
          "normal_scale", "pi_alpha", "pi_beta", "spike_var_shape",
          "spike_var_scale", "global_scale", "local_a", "local_b",
          "multi_gamma_list", "multi_pi_alpha_list", "multi_var_shape",
          "multi_var_scale", "learn_residual_var", "fixed_residual_var",
          "store_samples", "store_coefficient_cov", "effective_n",
          "fit_intercept", "intercept_x_mean", "intercept_y_mean",
          "compute_pve", "pve_type_code"
        )],
        list(
          ld_blocks = ld_blocks,
          ld_indices = ld_indices,
          ld_scale = ld_scale,
          summary_Xty = Xty,
          summary_yty = if (is.null(yty)) 0 else yty,
          nthreads = nthreads
        )
      )
    )
  } else if (!is.null(eigen_X)) {
    do.call(
      blm_gibbs_eigen_block_rcpp_cpp,
      c(
        common_arguments[c(
          "residual_shape", "residual_scale", "iterations", "burnin", "thin",
          "progress_callback", "block_id", "block_model", "normal_shape",
          "normal_scale", "pi_alpha", "pi_beta", "spike_var_shape",
          "spike_var_scale", "global_scale", "local_a", "local_b",
          "multi_gamma_list", "multi_pi_alpha_list", "multi_var_shape",
          "multi_var_scale", "learn_residual_var", "fixed_residual_var",
          "store_samples", "store_coefficient_cov", "effective_n",
          "fit_intercept", "intercept_x_mean", "intercept_y_mean",
          "compute_pve", "pve_type_code"
        )],
        list(
          transformed_X = eigen_X,
          transformed_y = eigen_y,
          eigen_indices = eigen_indices,
          summary_Xty = Xty,
          summary_yty = if (is.null(yty)) 0 else yty,
          nthreads = nthreads
        )
      )
    )
  } else if (is.list(XtX)) {
    do.call(
      blm_gibbs_block_rcpp_cpp,
      c(
        common_arguments[c(
          "residual_shape", "residual_scale", "iterations", "burnin", "thin",
          "progress_callback", "block_id", "block_model", "normal_shape",
          "normal_scale", "pi_alpha", "pi_beta", "spike_var_shape",
          "spike_var_scale", "global_scale", "local_a", "local_b",
          "multi_gamma_list", "multi_pi_alpha_list", "multi_var_shape",
          "multi_var_scale", "learn_residual_var", "fixed_residual_var",
          "store_samples", "store_coefficient_cov", "effective_n",
          "fit_intercept", "intercept_x_mean", "intercept_y_mean",
          "compute_pve", "pve_type_code"
        )],
        list(
          summary_XtX = XtX,
          summary_indices = XtX_indices,
          summary_types = XtX_types,
          summary_center = XtX_center,
          summary_Xty = Xty,
          summary_yty = if (is.null(yty)) 0 else yty,
          nthreads = nthreads
        )
      )
    )
  } else if (inherits(XtX, "dgCMatrix")) {
    do.call(
      blm_gibbs_sparse_rcpp_cpp,
      c(
        common_arguments[c(
          "residual_shape", "residual_scale", "iterations", "burnin", "thin",
          "progress_callback", "block_id", "block_model", "normal_shape",
          "normal_scale", "pi_alpha", "pi_beta", "spike_var_shape",
          "spike_var_scale", "global_scale", "local_a", "local_b",
          "multi_gamma_list", "multi_pi_alpha_list", "multi_var_shape",
          "multi_var_scale", "learn_residual_var", "fixed_residual_var",
          "store_samples", "store_coefficient_cov", "effective_n",
          "fit_intercept", "intercept_x_mean", "intercept_y_mean",
          "compute_pve", "pve_type_code"
        )],
        list(
          summary_XtX = XtX,
          summary_center = XtX_center,
          summary_Xty = Xty,
          summary_yty = if (is.null(yty)) 0 else yty
        )
      )
    )
  } else {
    common_arguments$summary_XtX <- if (use_sufficient_statistics) {
      XtX
    } else {
      matrix(0, 0L, 0L)
    }
    do.call(blm_gibbs_rcpp_cpp, common_arguments)
  }
  if (store_samples) {
    predictor_names <- if (use_sufficient_statistics) {
      names(Xty)
    } else {
      colnames(x)
    }
    predictor_model <- block_model[block_id]
    colnames(samples$coefficient_samples) <- predictor_names
    if (!any(block_model == 0L)) {
      samples$normal_var_samples <- NULL
    }
    if (any(block_model == 1L)) {
      colnames(samples$inclusion_samples) <-
        predictor_names[predictor_model == 1L]
    } else {
      samples$inclusion_samples <- NULL
      samples$pi_samples <- NULL
      samples$slab_var_samples <- NULL
    }
    if (any(block_model == 2L)) {
      colnames(samples$local_var_samples) <-
        predictor_names[predictor_model == 2L]
    } else {
      samples$local_var_samples <- NULL
      samples$tau_sq_samples <- NULL
    }
    if (any(block_model == 3L)) {
      colnames(samples$multi_component_samples) <-
        predictor_names[predictor_model == 3L]
    } else {
      samples$multi_component_samples <- NULL
      samples$multi_pi_samples <- NULL
      samples$multi_var_samples <- NULL
    }
  }
  samples
}

.chain_progress_callback <- function(progressor, chain, nchains) {
  if (is.null(progressor)) {
    return(NULL)
  }
  force(progressor)
  force(chain)
  force(nchains)
  function(amount, iteration) {
    progressor(
      amount = amount,
      message = sprintf(
        "Chain %d/%d: iteration %d",
        chain,
        nchains,
        iteration
      )
    )
  }
}

.run_blm_chains <- function(sampler_arguments, version, nchains,
                            block_model, progressor = NULL) {
  if (nchains == 1L) {
    sampler <- if (version == "Rcpp") .blm_gibbs_rcpp else .blm_gibbs
    return(do.call(
      sampler,
      c(
        sampler_arguments,
        list(
          progress_callback = .chain_progress_callback(
            progressor,
            chain = 1L,
            nchains = 1L
          )
        )
      )
    ))
  }

  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  future::plan(future::multisession, workers = nchains)

  chain_futures <- lapply(seq_len(nchains), function(chain) {
    chain_progress <- .chain_progress_callback(progressor, chain, nchains)
    future::future({
      namespace <- asNamespace("BayesLinReg")
      chain_sampler <- if (version == "Rcpp") {
        get(".blm_gibbs_rcpp", envir = namespace)
      } else {
        get(".blm_gibbs", envir = namespace)
      }
      do.call(
        chain_sampler,
        c(
          sampler_arguments,
          list(
            progress_callback = chain_progress
          )
        )
      )
    }, seed = TRUE)
  })
  chain_samples <- lapply(chain_futures, future::value)
  .combine_blm_chains(
    chain_samples,
    block_model = block_model,
    store_samples = sampler_arguments$store_samples,
    store_coefficient_cov = sampler_arguments$store_coefficient_cov,
    compute_pve = sampler_arguments$compute_pve
  )
}

.combine_blm_chains <- function(chain_samples, block_model = 0L,
                                store_samples = TRUE,
                                store_coefficient_cov = TRUE,
                                compute_pve = FALSE) {
  if (!store_samples) {
    summary_names <- c(
      "number_of_draws", "coefficient_sum", "coefficient_sum_sq",
      "intercept_sum", "intercept_sum_sq", "residual_var_sum",
      "residual_var_sum_sq"
    )
    if (compute_pve) {
      summary_names <- c(
        summary_names, "block_pve_sum", "block_pve_sum_sq",
        "total_pve_sum", "total_pve_sum_sq", "cross_block_pve_sum",
        "cross_block_pve_sum_sq"
      )
    }
    if (any(block_model == 0L)) {
      summary_names <- c(summary_names, "normal_var_sum", "normal_var_sum_sq")
    }
    if (any(block_model == 1L)) {
      summary_names <- c(
        summary_names, "inclusion_sum", "pi_sum", "pi_sum_sq",
        "slab_var_sum", "slab_var_sum_sq"
      )
    }
    if (any(block_model == 2L)) {
      summary_names <- c(
        summary_names, "local_var_sum", "local_var_sum_sq",
        "tau_sq_sum", "tau_sq_sum_sq"
      )
    }
    if (any(block_model == 3L)) {
      summary_names <- c(
        summary_names, "multi_var_sum", "multi_var_sum_sq"
      )
    }
    combined <- stats::setNames(lapply(summary_names, function(name) {
      Reduce(`+`, lapply(chain_samples, `[[`, name))
    }), summary_names)
    if (store_coefficient_cov) {
      combined$coefficient_crossprod <- lapply(
        seq_along(block_model),
        function(block) {
          Reduce(`+`, lapply(chain_samples, function(samples) {
            samples$coefficient_crossprod[[block]]
          }))
        }
      )
    }
    if (any(block_model == 3L)) {
      list_names <- c(
        "multi_component_sum", "multi_pi_sum", "multi_pi_sum_sq"
      )
      for (name in list_names) {
        combined[[name]] <- lapply(seq_along(block_model), function(block) {
          values <- lapply(chain_samples, function(samples) {
            samples[[name]][[block]]
          })
          if (all(vapply(values, is.null, logical(1)))) NULL else Reduce(`+`, values)
        })
      }
    }
    return(combined)
  }
  number_of_draws <- vapply(
    chain_samples,
    function(samples) nrow(samples$coefficient_samples),
    integer(1)
  )
  combined <- list(
    coefficient_samples = do.call(
      rbind,
      lapply(chain_samples, `[[`, "coefficient_samples")
    ),
    intercept_samples = unlist(
      lapply(chain_samples, `[[`, "intercept_samples"),
      use.names = FALSE
    ),
    residual_var_samples = unlist(
      lapply(chain_samples, `[[`, "residual_var_samples"),
      use.names = FALSE
    ),
    chain_id = rep.int(seq_along(chain_samples), number_of_draws)
  )
  if (compute_pve) {
    combined$block_pve_samples <- do.call(
      rbind, lapply(chain_samples, `[[`, "block_pve_samples")
    )
    combined$total_pve_samples <- unlist(
      lapply(chain_samples, `[[`, "total_pve_samples"), use.names = FALSE
    )
    combined$cross_block_pve_samples <- unlist(
      lapply(chain_samples, `[[`, "cross_block_pve_samples"),
      use.names = FALSE
    )
  }
  if (any(block_model == 0L)) {
    combined$normal_var_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "normal_var_samples")
    )
  }
  if (any(block_model == 1L)) {
    combined$inclusion_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "inclusion_samples")
    )
    combined$pi_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "pi_samples")
    )
    combined$slab_var_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "slab_var_samples")
    )
  }
  if (any(block_model == 2L)) {
    combined$local_var_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "local_var_samples")
    )
    combined$tau_sq_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "tau_sq_samples")
    )
  }
  if (any(block_model == 3L)) {
    combined$multi_component_samples <- do.call(
      rbind, lapply(chain_samples, `[[`, "multi_component_samples")
    )
    combined$multi_pi_samples <- lapply(
      seq_along(block_model),
      function(block) {
        if (block_model[block] != 3L) return(NULL)
        do.call(rbind, lapply(chain_samples, function(samples) {
          samples$multi_pi_samples[[block]]
        }))
      }
    )
    combined$multi_var_samples <- do.call(
      rbind, lapply(chain_samples, `[[`, "multi_var_samples")
    )
  }
  combined
}
