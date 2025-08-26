suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

as_list <- function(x) if (is.list(x)) x else list(x)
get_cox_death <- function(x) {
  if (inherits(x, "coxph")) return(x)
  if (is.list(x) && !is.null(x$death) && inherits(x$death, "coxph")) return(x$death)
  NULL
}

safe_num <- function(expr, default = NA_real_) {
  out <- try(expr, silent = TRUE)
  if (inherits(out, "try-error")) return(default)
  v <- suppressWarnings(as.numeric(out))
  if (!length(v)) return(default)
  v[1]
}

extract_model_metrics <- function(bayes_fit, cox_fit, fg_fit, strict = FALSE, context_label = NULL) {
  bayes_list <- as_list(bayes_fit)
  cox_list   <- as_list(cox_fit)
  fg_list    <- as_list(fg_fit)
  
  # ---------- INLA ----------
  bay_tbl <- purrr::map_dfr(bayes_list, function(obj) {
    mdl_d <- tryCatch(obj$death$model,     error = function(e) NULL)
    mdl_s <- tryCatch(obj$discharge$model, error = function(e) NULL)
    if (is.null(mdl_d) && is.null(mdl_s)) return(NULL)
    
    get_dic  <- function(mdl) if (!is.null(mdl) && !is.null(mdl$dic))  safe_num(mdl$dic$dic)  else NA_real_
    get_waic <- function(mdl) if (!is.null(mdl) && !is.null(mdl$waic)) safe_num(mdl$waic$waic) else NA_real_
    
    tau_row <- NULL
    if (!is.null(mdl_d) && !is.null(mdl_d$summary.hyperpar)) {
      ix <- grepl("Precision for study_grp", rownames(mdl_d$summary.hyperpar))
      if (any(ix)) tau_row <- mdl_d$summary.hyperpar[ix, , drop = FALSE]
    }
    
    if (!is.null(tau_row)) {
      var_hat <- safe_num(1 / tau_row["mean"])
      lwr     <- safe_num(1 / tau_row["0.975quant"])
      upr     <- safe_num(1 / tau_row["0.025quant"])
    } else {
      var_hat <- lwr <- upr <- NA_real_
    }
    
    tibble::tibble(
      DIC        = sum(get_dic(mdl_d), get_dic(mdl_s), na.rm = TRUE),
      WAIC       = sum(get_waic(mdl_d), get_waic(mdl_s), na.rm = TRUE),
      centre_var = var_hat,
      centre_lwr = lwr,
      centre_upr = upr
    )
  })
  have_inla <- nrow(bay_tbl) > 0
  
  bay_pooled <- if (!have_inla) NULL else tibble::tibble(
    model   = "Bayesian (INLA)",
    logLik  = NA_real_,
    AIC     = NA_real_,
    BIC     = NA_real_,
    DIC     = if (any(is.finite(bay_tbl$DIC)))  mean(bay_tbl$DIC,  na.rm = TRUE) else NA_real_,
    WAIC    = if (any(is.finite(bay_tbl$WAIC))) mean(bay_tbl$WAIC, na.rm = TRUE) else NA_real_,
    centre_var = if (any(is.finite(bay_tbl$centre_var)) &&
                     any(is.finite(bay_tbl$centre_lwr)) &&
                     any(is.finite(bay_tbl$centre_upr))) {
      sprintf("%.3f (%.3f–%.3f)",
              mean(bay_tbl$centre_var, na.rm = TRUE),
              mean(bay_tbl$centre_lwr, na.rm = TRUE),
              mean(bay_tbl$centre_upr, na.rm = TRUE))
    } else NA_character_,
    ridge_lambda = NA_character_
  )
  
  cox_tbl <- purrr::map_dfr(cox_list, function(obj) {
    mdl <- get_cox_death(obj)
    if (is.null(mdl)) return(NULL)
    tibble::tibble(
      logLik     = safe_num(logLik(mdl)),
      AIC        = safe_num(AIC(mdl)),
      BIC        = safe_num(BIC(mdl)),
      centre_var = {
        theta_val <- if (!is.null(mdl$theta)) safe_num(mdl$theta[1]) else NA_real_
        if (is.finite(theta_val)) theta_val^2 else NA_real_
      }
    )
  })
  have_cox <- nrow(cox_tbl) > 0
  
  cox_pooled <- if (!have_cox) NULL else tibble::tibble(
    model   = "Cox frailty",
    logLik  = if (any(is.finite(cox_tbl$logLik))) mean(cox_tbl$logLik, na.rm = TRUE) else NA_real_,
    AIC     = if (any(is.finite(cox_tbl$AIC)))    mean(cox_tbl$AIC,    na.rm = TRUE) else NA_real_,
    BIC     = if (any(is.finite(cox_tbl$BIC)))    mean(cox_tbl$BIC,    na.rm = TRUE) else NA_real_,
    DIC     = NA_real_,
    WAIC    = NA_real_,
    centre_var = if (any(is.finite(cox_tbl$centre_var))) sprintf(
      "%.3f (%.3f–%.3f)",
      mean(cox_tbl$centre_var, na.rm = TRUE),
      stats::quantile(cox_tbl$centre_var, 0.025, na.rm = TRUE),
      stats::quantile(cox_tbl$centre_var, 0.975, na.rm = TRUE)
    ) else NA_character_,
    ridge_lambda = NA_character_
  )
  
  fg_tbl <- purrr::map_dfr(fg_list, function(obj) {
    if (is.null(obj)) {
      return(tibble::tibble(engine = NA_character_,
                            logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, lambda = NA_real_))
    }
    if (inherits(obj, "FGR")) {
      ll  <- safe_num(obj$logLik)
      k   <- safe_num(length(stats::coef(obj)))
      n   <- safe_num(obj$n)
      lam <- safe_num(obj$lambda)
      tibble::tibble(
        engine = "FGR",
        logLik = ll,
        AIC    = if (is.finite(ll) && is.finite(k)) -2 * ll + 2 * k else NA_real_,
        BIC    = if (is.finite(ll) && is.finite(k) && is.finite(n)) -2 * ll + log(n) * k else NA_real_,
        lambda = lam
      )
    } else if (inherits(obj, "crr")) {
      ll <- safe_num(obj$loglik[2])
      k  <- safe_num(length(obj$coef))
      n  <- safe_num(obj$n)
      tibble::tibble(
        engine = "crr",
        logLik = ll,
        AIC    = if (is.finite(ll) && is.finite(k)) -2 * ll + 2 * k else NA_real_,
        BIC    = if (is.finite(ll) && is.finite(k) && is.finite(n)) -2 * ll + log(n) * k else NA_real_,
        lambda = NA_real_
      )
    } else {
      tibble::tibble(engine = paste(class(obj), collapse = "/"),
                     logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, lambda = NA_real_)
    }
  })
  
  have_fg_any <- nrow(fg_tbl) > 0 && any(is.finite(fg_tbl$logLik) | is.finite(fg_tbl$AIC) | is.finite(fg_tbl$BIC))
  
  if (!have_fg_any) {
    if (strict) {
      stop(
        paste0(
          "Fine–Gray models were not fitted in any imputation. ",
          "Typical reason: no censoring (non-identifiable). ",
          if (!is.null(context_label)) paste0("[context: ", context_label, "]") else ""
        ),
        call. = FALSE
      )
    }
    fg_rows <- tibble::tibble(
      model = "Fine–Gray (cluster): not fitted",
      logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
      DIC = NA_real_, WAIC = NA_real_,
      centre_var = NA_character_,
      ridge_lambda = NA_character_
    )
  } else {
    fg_rows <- list()
    if (any(fg_tbl$engine == "FGR", na.rm = TRUE)) {
      sub <- dplyr::filter(fg_tbl, engine == "FGR")
      fg_rows[["FGR"]] <- tibble::tibble(
        model = "Fine–Gray ridge",
        logLik = if (any(is.finite(sub$logLik))) mean(sub$logLik, na.rm = TRUE) else NA_real_,
        AIC    = if (any(is.finite(sub$AIC)))    mean(sub$AIC,    na.rm = TRUE) else NA_real_,
        BIC    = if (any(is.finite(sub$BIC)))    mean(sub$BIC,    na.rm = TRUE) else NA_real_,
        DIC    = NA_real_,
        WAIC   = NA_real_,
        centre_var   = NA_character_,
        ridge_lambda = if (any(is.finite(sub$lambda))) sprintf("%.3f", mean(sub$lambda, na.rm = TRUE)) else NA_character_
      )
    }
    if (any(fg_tbl$engine == "crr", na.rm = TRUE)) {
      sub <- dplyr::filter(fg_tbl, engine == "crr")
      fg_rows[["crr"]] <- tibble::tibble(
        model = "Fine–Gray (crr)",
        logLik = if (any(is.finite(sub$logLik))) mean(sub$logLik, na.rm = TRUE) else NA_real_,
        AIC    = if (any(is.finite(sub$AIC)))    mean(sub$AIC,    na.rm = TRUE) else NA_real_,
        BIC    = if (any(is.finite(sub$BIC)))    mean(sub$BIC,    na.rm = TRUE) else NA_real_,
        DIC    = NA_real_,
        WAIC   = NA_real_,
        centre_var   = NA_character_,
        ridge_lambda = NA_character_
      )
    }
    fg_rows <- dplyr::bind_rows(fg_rows)
  }
  
  levels_order <- c(
    "Bayesian (INLA)",
    "Cox frailty",
    "Fine–Gray ridge",
    "Fine–Gray (crr)",
    "Fine–Gray (cluster)",
    "Fine–Gray (cluster): not fitted"
  )
  
  pooled <- dplyr::bind_rows(bay_pooled, cox_pooled, fg_rows) %>%
    dplyr::mutate(model = factor(model, levels = levels_order)) %>%
    dplyr::arrange(model) %>%
    dplyr::mutate(model = as.character(model))
  
  bay_imp <- dplyr::bind_rows(bayes_list %>% purrr::imap(function(obj, i) {
    mdl_d <- tryCatch(obj$death$model,     error = function(e) NULL)
    mdl_s <- tryCatch(obj$discharge$model, error = function(e) NULL)
    tibble::tibble(
      .imp = i,
      family = "INLA",
      DIC  = if (!is.null(mdl_d) && !is.null(mdl_s))
        sum(suppressWarnings(mdl_d$dic$dic), suppressWarnings(mdl_s$dic$dic), na.rm = TRUE) else NA_real_,
      WAIC = if (!is.null(mdl_d) && !is.null(mdl_s))
        sum(suppressWarnings(mdl_d$waic$waic), suppressWarnings(mdl_s$waic$waic), na.rm = TRUE) else NA_real_
    )
  }))
  
  cox_imp <- dplyr::bind_rows(cox_list %>% purrr::imap(function(obj, i) {
    mdl <- get_cox_death(obj)
    tibble::tibble(
      .imp = i, family = "Cox frailty",
      logLik = if (!is.null(mdl)) suppressWarnings(as.numeric(logLik(mdl))) else NA_real_,
      AIC    = if (!is.null(mdl)) suppressWarnings(as.numeric(AIC(mdl)))    else NA_real_,
      BIC    = if (!is.null(mdl)) suppressWarnings(as.numeric(BIC(mdl)))    else NA_real_
    )
  }))
  
  fg_imp <- dplyr::bind_rows(fg_list %>% purrr::imap(function(obj, i) {
    if (is.null(obj)) {
      tibble::tibble(.imp = i, family = "FG", engine = NA_character_,
                     logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, lambda = NA_real_)
    } else if (inherits(obj, "FGR")) {
      tibble::tibble(.imp = i, family = "FG", engine = "FGR",
                     logLik = suppressWarnings(as.numeric(obj$logLik)),
                     AIC = NA_real_, BIC = NA_real_, lambda = suppressWarnings(as.numeric(obj$lambda)))
    } else if (inherits(obj, "crr")) {
      ll <- try(obj$loglik[2], silent = TRUE); ll <- if (inherits(ll,"try-error")) NA_real_ else as.numeric(ll)
      k  <- try(length(obj$coef), silent = TRUE); k <- if (inherits(k,"try-error")) NA_real_ else as.numeric(k)
      n  <- try(obj$n,           silent = TRUE); n <- if (inherits(n,"try-error")) NA_real_ else as.numeric(n)
      tibble::tibble(.imp = i, family = "FG", engine = "crr",
                     logLik = ll,
                     AIC = if (is.finite(ll) && is.finite(k)) -2*ll + 2*k else NA_real_,
                     BIC = if (is.finite(ll) && is.finite(k) && is.finite(n)) -2*ll + log(n)*k else NA_real_,
                     lambda = NA_real_)
    } else {
      tibble::tibble(.imp = i, family = "FG", engine = paste(class(obj), collapse = "/"),
                     logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, lambda = NA_real_)
    }
  }))
  
  per_imputation <- dplyr::bind_rows(bay_imp, cox_imp, fg_imp)
  
  diagnostics <- list(
    context = context_label,
    n_imputations = length(bayes_list),
    fg_counts = dplyr::count(fg_tbl, .data$engine),
    reason_fg_unfitted = if (!have_fg_any)
      "No censoring in distinct-patient clustered data; FG not identifiable."
    else NA_character_,
    timestamp = Sys.time()
  )
  
  out <- list(
    pooled = pooled,
    per_imputation = per_imputation,
    diagnostics = diagnostics
  )
  class(out) <- c("metrics_bundle", class(out))
  out
}

auc_simple <- function(score, y) {
  ok <- is.finite(score) & is.finite(y)
  y <- y[ok]; score <- score[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

brier_simple <- function(score, y) {
  ok <- is.finite(score) & is.finite(y)
  mean((score[ok] - y[ok])^2)
}

compute_raw_metrics <- function(cp_list, bayes_list, cox_list, fg_list,
                                horizons = c(0.5,1,3,7,15,30,60,90)) {
  stopifnot(length(cp_list) == length(bayes_list),
            length(cp_list) == length(cox_list),
            length(cp_list) == length(fg_list))
  
  one_imp_h <- function(i, h) {
    cp <- dplyr::distinct(cp_list[[i]], patient, .keep_all = TRUE)
    y  <- as.integer(cp$status == 1L & cp$stop <= h)
    
    pin <- predict_inla_cif(bayes_list[[i]], cp, h)$death
    pcx <- predict_cox_cif( cox_list[[i]],   cp, h)
    pfg <- predict_fg_cif(  fg_list[[i]],    cp, h)
    
    tibble::tibble(
      .imp   = i,
      Horizon= h,
      model  = c("INLA","Cox","FG"),
      AUC    = c(auc_simple(pin,y), auc_simple(pcx,y), auc_simple(pfg,y)),
      Brier  = c(brier_simple(pin,y), brier_simple(pcx,y), brier_simple(pfg,y))
    )
  }
  
  purrr::map_dfr(seq_along(cp_list), function(i) {
    purrr::map_dfr(horizons, function(h) one_imp_h(i, h))
  })
}

compute_pooled_metrics <- function(raw_metrics) {
  raw_metrics |>
    dplyr::group_by(model, Horizon) |>
    dplyr::summarise(
      AUC   = mean(AUC,   na.rm = TRUE),
      Brier = mean(Brier, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_cindex_td <- function(cp_list, bayes_list, cox_list, fg_list,
                              horizons = TIME_HORIZONS) {

    
    auc_simple <- function(score, y) {
      ok <- is.finite(score) & is.finite(y)
      y <- y[ok]; score <- score[ok]
      n1 <- sum(y == 1); n0 <- sum(y == 0)
      if (n1 == 0 || n0 == 0) return(NA_real_)
      r <- rank(score, ties.method = "average")
      (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    }
    
    c_td_one <- function(marker, cp, h) {
      ok <- is.finite(marker) & is.finite(cp$stop) & is.finite(cp$status)
      if (sum(ok) < 5) return(NA_real_)
      T      <- as.numeric(cp$stop[ok])
      delta  <- as.integer(cp$status[ok])   
      marker <- as.numeric(marker[ok])
      
      y_h <- as.integer(delta == 1L & T <= h)
      if (sum(y_h) == 0L || sum(1 - y_h) == 0L) return(NA_real_)
      
      out <- try(
        timeROC::timeROC(
          T = T, delta = delta, cause = 1, marker = marker,
          weighting = "marginal", times = h, iid = FALSE
        ),
        silent = TRUE
      )
      auc_td <- if (inherits(out, "try-error")) NA_real_ else suppressWarnings(as.numeric(out$AUC[1]))
      if (!length(auc_td) || !is.finite(auc_td)) auc_simple(marker, y_h) else auc_td
    }
    
    cindex_td_imp <- function(i, hvec = horizons) {
      cp <- dplyr::distinct(cp_list[[i]], patient, .keep_all = TRUE) |>
        dplyr::mutate(stop = as.numeric(stop), status = as.integer(status))
      
      purrr::map_dfr(hvec, function(h) {
        pin <- try(predict_inla_cif(bayes_list[[i]], cp, h)$death, silent = TRUE)
        pcx <- try(predict_cox_cif( cox_list[[i]],   cp, h),      silent = TRUE)
        pfg <- try(predict_fg_cif(  fg_list[[i]],    cp, h),      silent = TRUE)
        
        if (inherits(pin, "try-error")) pin <- rep(NA_real_, nrow(cp))
        if (inherits(pcx, "try-error")) pcx <- rep(NA_real_, nrow(cp))
        if (inherits(pfg, "try-error")) pfg <- rep(NA_real_, nrow(cp))
        
        tibble::tibble(
          Horizon = as.numeric(h),
          INLA = c_td_one(pin, cp, h),
          Cox  = c_td_one(pcx, cp, h),
          FG   = c_td_one(pfg,  cp, h)
        ) |>
          tidyr::pivot_longer(-Horizon, names_to = "model", values_to = "C_td")
      })
    }
    
    raw <- purrr::map_dfr(seq_along(cp_list), ~cindex_td_imp(.x, horizons))
    
    if (!nrow(raw)) {
      return(list(
        raw    = tibble::tibble(Horizon = numeric(), model = character(), C_td = numeric()),
        pooled = tibble::tibble(Horizon = numeric(), model = character(), C_td = numeric())
      ))
    }
    
    raw$model <- factor(raw$model, levels = c("INLA", "Cox", "FG"))
    
    pooled <- raw |>
      dplyr::group_by(Horizon, model) |>
      dplyr::summarise(C_td = mean(C_td, na.rm = TRUE), .groups = "drop")
    
    list(raw = raw, pooled = pooled)
}


compute_group_metrics <- function(cp_list, bayes_list, cox_list, fg_list,
                                  horizons = TIME_HORIZONS,
                                  min_n = 20) {
  auc_simple   <- function(s,y){ ok <- is.finite(s)&is.finite(y); y<-y[ok]; s<-s[ok]
  n1<-sum(y==1); n0<-sum(y==0); if(n1==0||n0==0) return(NA_real_)
  r<-rank(s,ties.method="average"); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0)
  }
  brier_simple <- function(s,y){ ok <- is.finite(s)&is.finite(y); mean((s[ok]-y[ok])^2) }
  
  one_imp_one_group <- function(i, g) {
    dat_all <- dplyr::distinct(cp_list[[i]], patient, .keep_all = TRUE)
    sub     <- dplyr::filter(dat_all, study_grp == g)
    if (nrow(sub) < min_n) return(tibble::tibble())  # skip tiny groups
    
    purrr::map_dfr(horizons, function(h) {
      pin <- try(predict_inla_cif(bayes_list[[i]], sub, h)$death, silent = TRUE)
      pcx <- try(predict_cox_cif( cox_list[[i]],   sub, h),      silent = TRUE)
      pfg <- try(predict_fg_cif(  fg_list[[i]],    sub, h),      silent = TRUE)
      if (inherits(pin,"try-error")) pin <- rep(NA_real_, nrow(sub))
      if (inherits(pcx,"try-error")) pcx <- rep(NA_real_, nrow(sub))
      if (inherits(pfg,"try-error")) pfg <- rep(NA_real_, nrow(sub))
      
      y <- as.integer(sub$status == 1L & sub$stop <= h)
      
      tibble::tibble(
        study_grp = g, .imp = i, Horizon = h,
        INLA_AUC = auc_simple(pin, y), Cox_AUC = auc_simple(pcx, y), FG_AUC = auc_simple(pfg, y),
        INLA_Brier = brier_simple(pin, y), Cox_Brier = brier_simple(pcx, y), FG_Brier = brier_simple(pfg, y),
        n = nrow(sub)
      )
    })
  }
  
  raw <- purrr::map_dfr(seq_along(cp_list), function(i) {
    grps <- unique(dplyr::distinct(cp_list[[i]], patient, .keep_all = TRUE)$study_grp)
    purrr::map_dfr(grps, ~one_imp_one_group(i, .x))
  })
  
  long_auc <- raw |>
    tidyr::pivot_longer(c(INLA_AUC, Cox_AUC, FG_AUC),
                        names_to = "model", values_to = "AUC") |>
    dplyr::mutate(model = sub("_AUC$", "", model))
  long_bri <- raw |>
    tidyr::pivot_longer(c(INLA_Brier, Cox_Brier, FG_Brier),
                        names_to = "model", values_to = "Brier") |>
    dplyr::mutate(model = sub("_Brier$", "", model))
  long <- dplyr::inner_join(
    dplyr::select(long_auc, study_grp, Horizon, .imp, model, AUC, n),
    dplyr::select(long_bri, study_grp, Horizon, .imp, model, Brier),
    by = c("study_grp","Horizon",".imp","model")
  )
  
  pooled <- long |>
    dplyr::group_by(study_grp, Horizon, model) |>
    dplyr::summarise(
      AUC   = mean(AUC,   na.rm = TRUE),
      Brier = mean(Brier, na.rm = TRUE),
      n     = max(n, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(raw = long, pooled = pooled)
}

  