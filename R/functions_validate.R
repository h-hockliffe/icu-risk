suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(survival)
  library(timeROC)
  library(riskRegression)
  library(splines)
})

options(brier_boot_B = 200)
if (!exists("rubin_pool", mode = "function")) {
  rubin_pool <- function(q, u = NULL) {
    ok <- is.finite(q)
    m <- sum(ok)
    if (m == 0) {
      return(list(qbar = NA_real_, tvar = NA_real_))
    }
    q <- q[ok]

    if (is.null(u)) {
      u <- rep(0, m)
    } else {
      u <- u[ok]
      u[!is.finite(u)] <- 0
    }

    qbar <- mean(q)
    wbar <- mean(u)
    bvar <- stats::var(q)

    tvar <- wbar + (1 + 1 / m) * bvar
    list(qbar = qbar, tvar = tvar)
  }
}

with_progress_if <- function(expr, n_steps = NULL) {
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress(expr)
  } else {
    force(expr)
  }
}

tick <- local({
  last_time <- Sys.time()
  max_gap <- getOption("validate_tick_secs", 60) 
  function(b, B) {
    beat <- max(1, floor(B / 50)) 
    now <- Sys.time()
    if (b == 1 ||
      b %% beat == 0 ||
      difftime(now, last_time, units = "secs") >= max_gap) {
      elapsed <- difftime(now, global_start, units = "mins")
      eta <- (elapsed / b) * (B - b)
      message(sprintf(
        "  … completed bootstrap %d / %d  (%.1f min elapsed, ETA %.1f min)",
        b, B, as.numeric(elapsed), as.numeric(eta)
      ))
      last_time <<- now
    }
  }
})

validate_models <- function(
    cp_list,
    bayes_list,
    cox_list,
    fg_list,
    B = getOption("val_B", 50),
    time_horizon = getOption("val_horizon", 30),
    min_events = 5,
    seed = 1000,
    seeds_file = file.path("_targets", "seeds_boot.rds"),
    brier_mode = c("temporal", "none", "both")
) {
  brier_mode <- match.arg(brier_mode)
  
  time_horizon <- sort(unique(time_horizon))
  if (!length(time_horizon))
    stop("validate_models(): `time_horizon` is empty.")
  
  if (length(cp_list) == 0L) {
    stop("cp_list is empty – nothing to validate.", call. = FALSE)
  }
  
  start <- Sys.time()
  message(
    format(start, "%Y-%m-%d %H:%M:%S"),
    " — validate_models(): starting (B = ", B,
    ", imputations = ", length(cp_list), ")"
  )
  
  if (!is.null(seeds_file) && file.exists(seeds_file)) {
    seeds_vec <- readRDS(seeds_file)
    if (length(seeds_vec) < B) {
      stop("seeds_file exists but has fewer than B seeds – please regenerate.",
           call. = FALSE
      )
    }
  } else {
    set.seed(seed)
    seeds_vec <- sample.int(1e9, B)
    if (!is.null(seeds_file)) {
      dir.create(dirname(seeds_file), recursive = TRUE, showWarnings = FALSE)
      saveRDS(seeds_vec, seeds_file)
      message("✓ Bootstrap seed vector saved to: ", seeds_file)
    }
  }
  
  compute_brier <- function(test_df, risk_INLA, risk_Cox, risk_FG, h) {
    objs <- list()
    if (any(is.finite(risk_INLA))) objs$INLA <- risk_INLA
    if (any(is.finite(risk_Cox)))  objs$Cox  <- risk_Cox
    if (any(is.finite(risk_FG)))   objs$FG   <- risk_FG
    if (length(objs) == 0L) {
      return(list(
        val = c(INLA = NA_real_, Cox = NA_real_, FG = NA_real_),
        se  = c(INLA = NA_real_, Cox = NA_real_, FG = NA_real_)
      ))
    }
    
    sc <- tryCatch(
      riskRegression::Score(
        object       = objs,
        formula      = prodlim::Hist(stop, status) ~ 1,
        data         = test_df,
        times        = h,
        cause        = 1,
        metrics      = "Brier",
        conf.int     = FALSE,  
        B            = 0,
        cens.model   = "km",
        split.method = "none"
      ),
      error = function(e) NULL
    )
    if (is.null(sc)) {
      return(list(
        val = c(INLA = NA_real_, Cox = NA_real_, FG = NA_real_),
        se  = c(INLA = NA_real_, Cox = NA_real_, FG = NA_real_)
      ))
    }
    
    tb   <- sc$Brier$score
    tb_h <- tb[tb$times == h & tb$type == "Brier", , drop = FALSE]
    
    get_val <- function(model) {
      if (model %in% tb_h$model) tb_h$Brier[tb_h$model == model][1] else NA_real_
    }
    
    list(
      val = c(INLA = get_val("INLA"), Cox = get_val("Cox"), FG = get_val("FG")),
      se  = c(INLA = NA_real_,        Cox = NA_real_,       FG = NA_real_)
    )
  }
  
  eval_split <- function(test, inla_fit, cox_fit, fg_fit) {
    test <- dplyr::distinct(test, patient, .keep_all = TRUE)
    
    if (nrow(test) == 0L) {
      return(
        expand.grid(
          horizon  = time_horizon,
          model    = c("INLA", "Cox", "FG")
        ) |>
          dplyr::mutate(
            c_index = NA_real_, c_index_se = NA_real_,
            intercept = NA_real_, intercept_se = NA_real_,
            slope = NA_real_, slope_se = NA_real_,
            brier = NA_real_, brier_se = NA_real_
          )
      )
    }
    
    if (sum(test$status == 1L, na.rm = TRUE) < min_events ||
        sum(test$status == 2L, na.rm = TRUE) < min_events) {
      return(
        expand.grid(
          horizon  = time_horizon,
          model    = c("INLA", "Cox", "FG")
        ) |>
          dplyr::mutate(
            c_index = NA_real_, c_index_se = NA_real_,
            intercept = NA_real_, intercept_se = NA_real_,
            slope = NA_real_, slope_se = NA_real_,
            brier = NA_real_, brier_se = NA_real_
          )
      )
    }
    
    score_one <- function(risk_vec, label, h, brier_val, brier_se) {
      ok <- is.finite(risk_vec)
      if (!any(ok)) {
        return(tibble::tibble(
          horizon = h,
          model = label,
          c_index = NA_real_, c_index_se = NA_real_,
          intercept = NA_real_, intercept_se = NA_real_,
          slope = NA_real_, slope_se = NA_real_,
          brier = brier_val, brier_se = brier_se
        ))
      }
      
      c_auc <- NA_real_; c_auc_se <- NA_real_
      tT   <- test$stop[ok]
      tSt  <- test$status[ok]
      tMrk <- risk_vec[ok]
      
      if (sum(tSt == 1L & tT <= h, na.rm = TRUE) > 0L) {
        roc <- tryCatch(
          timeROC::timeROC(
            T        = tT,
            delta    = tSt,
            marker   = tMrk,
            cause    = 1,
            weighting = "marginal",
            times     = h,
            iid       = FALSE           
          ),
          error = function(e) NULL
        )
        if (!is.null(roc)) {
          c_auc    <- suppressWarnings(as.numeric(roc$AUC[1]))
          c_auc_se <- NA_real_
          if (!is.finite(c_auc))  c_auc <- NA_real_
        }
      }
      
      obs_event <- as.integer((test$status == 1L) & (test$stop <= h))
      
      p  <- pmin(pmax(tMrk, 1e-6), 1 - 1e-6)
      lz <- qlogis(p)
      
      int_cal <- NA_real_; int_se <- NA_real_
      cil_fit <- tryCatch(
        stats::glm(obs_event[ok] ~ 1, offset = lz[ok], family = stats::binomial()),
        error = function(e) NULL
      )
      if (!is.null(cil_fit)) {
        cf <- stats::coef(cil_fit)
        vc <- stats::vcov(cil_fit)
        int_cal <- unname(if (!is.na(cf[1])) cf[1] else NA_real_)
        int_se  <- if (!anyNA(vc)) sqrt(vc[1,1]) else NA_real_
      }
      
      slope_val <- NA_real_; slope_se <- NA_real_
      sl_fit <- tryCatch(
        stats::glm(obs_event[ok] ~ lz[ok], family = stats::binomial()),
        error = function(e) NULL
      )
      if (!is.null(sl_fit)) {
        cf <- stats::coef(sl_fit)
        vc <- stats::vcov(sl_fit)
        slope_val <- unname(if (length(cf) >= 2) cf[2] else NA_real_)
        slope_se  <- if (!anyNA(vc) && ncol(vc) >= 2) sqrt(vc[2,2]) else NA_real_
      }
      
      tibble::tibble(
        horizon   = h,
        model     = label,
        c_index   = c_auc,
        c_index_se = c_auc_se,
        intercept = int_cal,
        intercept_se = int_se,
        slope       = slope_val,
        slope_se    = slope_se,
        brier       = brier_val,
        brier_se    = brier_se
      )
    }
    
    purrr::map_dfr(time_horizon, function(h) {
      risk_inla <- tryCatch(predict_inla_cif(inla_fit, test, h)$death,
                            error = function(e) rep(NA_real_, nrow(test)))
      risk_cox  <- tryCatch(predict_cox_cif(cox_fit,  test, h),
                            error = function(e) rep(NA_real_, nrow(test)))
      risk_fg   <- tryCatch(predict_fg_cif(fg_fit,   test, h),
                            error = function(e) rep(NA_real_, nrow(test)))
      
      brier_out <- compute_brier(test, risk_inla, risk_cox, risk_fg, h)
      
      dplyr::bind_rows(
        score_one(risk_inla, "INLA", h,
                  brier_val = brier_out$val["INLA"], brier_se = brier_out$se["INLA"]),
        score_one(risk_cox,  "Cox",  h,
                  brier_val = brier_out$val["Cox"],  brier_se = brier_out$se["Cox"]),
        score_one(risk_fg,   "FG",   h,
                  brier_val = brier_out$val["FG"],   brier_se = brier_out$se["FG"])
      )
    })
  }
  
  bootstrap_stats <- with_progress_if(
    {
      prog <- if (requireNamespace("progressr", quietly = TRUE)) {
        progressr::progressor(steps = B)
      } else {
        function(...) NULL
      }
      
      purrr::map_dfr(seq_len(B), function(b) {
        prog()
        set.seed(seeds_vec[b])
        
        res <- tryCatch(
          {
            all_blocks   <- unique(cp_list[[1]]$study_grp)
            train_blocks <- sample(all_blocks, length(all_blocks), replace = TRUE)
            
            imp_res <- purrr::imap_dfr(cp_list, function(dat, nm) {
              test_ids <- dat$patient[!dat$study_grp %in% train_blocks]
              eval_split(
                test     = dat[dat$patient %in% test_ids, ],
                inla_fit = bayes_list[[nm]],
                cox_fit  = cox_list[[nm]],
                fg_fit   = fg_list[[nm]]
              ) |>
                dplyr::mutate(.imp = nm)
            })
            
            stopifnot(all(c("horizon","c_index","c_index_se","intercept","intercept_se",
                            "slope","slope_se","brier","brier_se") %in% names(imp_res)))
            
            pooled <- imp_res |>
              dplyr::group_by(model, horizon) |>
              dplyr::summarise(
                across(
                  c(c_index, intercept, slope, brier,
                    c_index_se, intercept_se, slope_se, brier_se),
                  ~ list(.),
                  .names = "{.col}"
                ),
                .groups = "drop"
              ) |>
              dplyr::rowwise() |>
              dplyr::mutate(
                c_index   = rubin_pool(unlist(c_index),   unlist(c_index_se)^2)$qbar,
                intercept = rubin_pool(unlist(intercept), unlist(intercept_se)^2)$qbar,
                slope     = rubin_pool(unlist(slope),     unlist(slope_se)^2)$qbar,
                brier     = if (identical(brier_mode, "both")) {
                  b <- unlist(brier); se <- unlist(brier_se)
                  if (any(is.finite(se))) rubin_pool(b, se^2)$qbar else mean(b, na.rm = TRUE)
                } else NA_real_
              ) |>
              dplyr::ungroup() |>
              dplyr::transmute(
                model, horizon, c_index, intercept, slope,
                brier = if (identical(brier_mode, "both")) brier else NA_real_,
                split = paste0("boot", b)
              )
            
            pooled
          },
          error = function(e) {
            warning("Bootstrap ", b, " failed: ", e$message)
            tibble::tibble(
              model = rep(c("INLA", "Cox", "FG"), each = length(time_horizon)),
              horizon = rep(time_horizon, times = 3),
              c_index = NA_real_,
              intercept = NA_real_,
              slope = NA_real_,
              brier = if (identical(brier_mode, "both")) NA_real_ else NA_real_,
              split = paste0("boot", b)
            )
          }
        )
        
        res
      })
    },
    n_steps = B
  )
  
  temporal_stats <- {
    imp_res <- purrr::imap_dfr(cp_list, function(dat, nm) {
      hold_ids <- dat$patient[lubridate::year(dat$icu_admit) >= 2022]
      eval_split(
        test     = dat[dat$patient %in% hold_ids, ],
        inla_fit = bayes_list[[nm]],
        cox_fit  = cox_list[[nm]],
        fg_fit   = fg_list[[nm]]
      ) |>
        dplyr::mutate(.imp = nm)
    })
    
    stopifnot(all(c("horizon","c_index","c_index_se","intercept","intercept_se",
                    "slope","slope_se","brier","brier_se") %in% names(imp_res)))
    
    imp_res |>
      dplyr::group_by(model, horizon) |>
      dplyr::summarise(
        across(
          c(c_index, intercept, slope, brier, c_index_se, intercept_se, slope_se, brier_se),
          ~ list(.),
          .names = "{.col}"
        ),
        .groups = "drop"
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        c_index   = rubin_pool(unlist(c_index),   unlist(c_index_se)^2)$qbar,
        intercept = rubin_pool(unlist(intercept), unlist(intercept_se)^2)$qbar,
        slope     = rubin_pool(unlist(slope),     unlist(slope_se)^2)$qbar,
        brier     = if (identical(brier_mode, "none")) {
          NA_real_
        } else {
          b <- unlist(brier); se <- unlist(brier_se)
          if (any(is.finite(se))) rubin_pool(b, se^2)$qbar else mean(b, na.rm = TRUE)
        }
      ) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        model, horizon, c_index, intercept, slope, brier,
        split = "temporal"
      )
  }
  
  out <- list(
    bootstrap = bootstrap_stats,
    temporal  = temporal_stats
  )
  attr(out, "bootstrap_seeds") <- seeds_vec
  out
}


cox_ph_diagnostics <- function(cox_list, transform = "rank") {
  empty <- tibble::tibble(
    .imp  = integer(),
    term  = character(),
    rho   = numeric(),
    chisq = numeric(),
    p     = numeric()
  )
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  pick_col <- function(mat, opts) {
    cn <- tolower(colnames(mat) %||% character())
    idx <- which(cn %in% tolower(opts))[1]
    if (!length(idx) || is.na(idx)) return(rep(NA_real_, nrow(mat)))
    as.numeric(mat[, idx])
  }
  
  get_model <- function(obj) {
    if (inherits(obj, "coxph")) return(obj)
    if (is.list(obj) && inherits(obj$death, "coxph")) return(obj$death)
    NULL
  }
  
  purrr::imap_dfr(cox_list, function(obj, nm) {
    mdl <- get_model(obj)
    if (is.null(mdl)) return(empty)
    
    z <- NULL
    
    ph_form <- attr(mdl, "ph_form")
    ph_df   <- attr(mdl, "ph_df")
    if (!is.null(ph_form) && !is.null(ph_df)) {
      fit <- try(
        survival::coxph(ph_form, data = ph_df, ties = "breslow",
                        x = FALSE, y = FALSE, model = TRUE),
        silent = TRUE
      )
      if (!inherits(fit, "try-error")) {
        z <- try(survival::cox.zph(fit, transform = transform), silent = TRUE)
      }
    }
    
    if (is.null(z) || inherits(z, "try-error")) {
      Yc <- attr(mdl, "ph_y_counting")
      X  <- mdl$x
      if (!is.null(Yc) && !is.null(X)) {
        nx <- NROW(X); ny <- NROW(Yc); n <- min(nx, ny)
        if (is.finite(n) && n > 0L) {
          if (ny != n) Yc <- Yc[seq_len(n), , drop = FALSE]
          if (nx != n) X  <- X [seq_len(n), , drop = FALSE]
          
          df <- cbind.data.frame(
            start = Yc[, 1],
            stop  = Yc[, 2],
            event = Yc[, 3],
            as.data.frame(X, check.names = FALSE)
          )
          
          orig   <- attr(mdl, "ph_x_colnames") %||% colnames(X)
          tt_col <- if (!is.null(orig)) grepl("(?i)^tt\\s*\\(", orig) else rep(FALSE, ncol(X))
          
          covars <- setdiff(names(df), c("start", "stop", "event"))
          covars <- covars[!tt_col]
          if (length(covars)) {
            keep <- vapply(df[covars], function(v) {
              v <- suppressWarnings(as.numeric(v))
              any(is.finite(v)) && stats::sd(v, na.rm = TRUE) > 0
            }, logical(1))
            covars <- covars[keep]
          }
          
          if (length(covars)) {
            formB <- stats::as.formula(
              paste("survival::Surv(start, stop, event) ~", paste(covars, collapse = " + "))
            )
            fitB <- try(
              survival::coxph(formB,
                              data  = df[, c("start","stop","event", covars), drop = FALSE],
                              ties  = "breslow",
                              x     = FALSE, y = FALSE, model = TRUE),
              silent = TRUE
            )
            if (!inherits(fitB, "try-error")) {
              z <- try(survival::cox.zph(fitB, transform = transform), silent = TRUE)
            }
          }
        }
      }
    }
    
    if (is.null(z) || inherits(z, "try-error")) return(empty)
    
    tt <- as.matrix(z$table)
    if (is.null(colnames(tt))) colnames(tt) <- paste0("V", seq_len(ncol(tt)))
    rn <- rownames(tt) %||% paste0("term", seq_len(nrow(tt)))
    
    drop_pat  <- "(?i)^(frailty|cluster|strata)\\b|^tt\\s*\\("
    keep_rows <- !grepl(drop_pat, rn, perl = TRUE)
    tt <- tt[keep_rows, , drop = FALSE]
    rn <- rn[keep_rows]
    if (!nrow(tt)) return(empty)
    
    imp_id <- NA_integer_
    if (!is.null(ph_df) && ".imp" %in% names(ph_df)) {
      imp_id <- as.integer(unique(ph_df$.imp)[1])
    } else {
      maybe <- suppressWarnings(as.integer(gsub(".*_(\\d+)$", "\\1", nm)))
      if (is.finite(maybe)) imp_id <- maybe
    }
    
    tibble::tibble(
      .imp  = as.integer(imp_id),
      term  = rn,
      rho   = pick_col(tt, c("rho", "cor", "corr")),
      chisq = pick_col(tt, c("chisq", "chi.sq", "chi-square")),
      p     = pick_col(tt, c("p", "pvalue", "p.value", "pr(>|z|)"))
    )
  })
}
