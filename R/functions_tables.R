suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})


as_list <- function(x) if (is.list(x)) x else list(x)

rubin_pool <- function(q, u) {
  ok <- is.finite(q) & is.finite(u)
  m <- sum(ok)
  if (m == 0) {
    return(list(qbar = NA_real_, tvar = NA_real_))
  }
  q <- q[ok]
  u <- u[ok]
  qbar <- mean(q)
  wbar <- mean(u)
  bvar <- stats::var(q)
  tvar <- wbar + (1 + 1 / m) * bvar
  list(qbar = qbar, tvar = tvar)
}

fmt_ci <- function(est, lwr, upr) {
  sprintf("%.2f (%.2f–%.2f)", est, lwr, upr)
}

extract_death_model <- function(obj, pattern) {
  if (inherits(obj, pattern)) {
    obj
  } else if (is.list(obj) && !is.null(obj$death)) {
    obj$death
  } else {
    NULL
  }
}

nice_term <- function(x) {
  recode(x,
         `log_egfr`            = "log(eGFR)",
         `log(egfr)`           = "log(eGFR)",
         `log_bmi`             = "log(BMI)",
         `log(bmi)`            = "log(BMI)",
         `sexVrouw`            = "Female vs Male",
         `sex`                 = "Female vs Male",
         `cr_resid`            = "Creatinine residual (⊥ eGFR)",
         `mega_pool`           = "Merged-small-sites stratum",
         .default = x
  )
}

make_coef_table <- function(inla_obj, cox_obj, fg_obj) {
  gather_coefs <- function(obj_list, extractor) {
    dplyr::bind_rows(
      lapply(seq_along(obj_list), function(i) {
        mdl <- extractor(obj_list[[i]])
        if (is.null(mdl)) return(NULL)
        
        beta_vec <- tryCatch(stats::coef(mdl), error = function(e) NULL)
        vcov_m   <- tryCatch(stats::vcov(mdl), error = function(e) NULL)
        if (!is.null(beta_vec) && !is.null(vcov_m)) {
          se_vec <- sqrt(diag(vcov_m))
          return(tibble::tibble(
            .imp = i,
            term = names(beta_vec),
            beta = as.numeric(beta_vec),
            var  = as.numeric(se_vec)^2
          ))
        }
        NULL
      })
    )
  }
  
  pool_family <- function(df, label) {
    if (nrow(df) == 0) return(NULL)
    
    pooled <- df %>%
      dplyr::group_by(term) %>%
      dplyr::summarise(pool = list(rubin_pool(beta, var)), .groups = "drop") %>%
      dplyr::mutate(
        beta = vapply(pool, `[[`, numeric(1), "qbar"),
        se   = sqrt(vapply(pool, `[[`, numeric(1), "tvar")),
        hr   = fmt_ci(
          exp(beta),
          exp(beta - 1.96 * se),
          exp(beta + 1.96 * se)
        ),
        p_raw = 2 * stats::pnorm(-abs(beta / pmax(se, .Machine$double.eps)))
      ) %>%
      dplyr::select(
        term,
        !!paste0(label, "_hr") := hr,
        !!paste0(label, "_p")  := p_raw
      )
    
    pooled
  }
  
  inla_list <- as_list(inla_obj)
  inla_df <- dplyr::bind_rows(
    lapply(seq_along(inla_list), function(i) {
      mdl <- inla_list[[i]]
      if (is.list(mdl) && !is.null(mdl$death) && !is.null(mdl$death$model)) mdl <- mdl$death$model
      if (is.null(mdl) || !inherits(mdl, "inla")) return(NULL)
      sf <- mdl$summary.fixed
      if (is.null(sf) || !nrow(sf)) return(NULL)
      tibble::tibble(
        .imp = i,
        term = rownames(sf),
        beta = as.numeric(sf$mean),
        var  = as.numeric(sf$sd)^2
      )
    })
  )
  
  cox_df <- gather_coefs(as_list(cox_obj), function(x) {
    if (inherits(x, "coxph")) return(x)
    if (is.list(x) && !is.null(x$death) && inherits(x$death, "coxph")) return(x$death)
    NULL
  })
  
  fg_list <- as_list(fg_obj)
  
  fg_engines <- vapply(fg_list, function(o) {
    if (is.null(o)) return(NA_character_)
    if (inherits(o, "FGR")) return("FG ridge")
    if (inherits(o, "crr")) return("FG (crr)")
    NA_character_
  }, character(1))
  
  fg_label <- {
    labs <- unique(na.omit(fg_engines))
    if (length(labs) == 1) labs else "FG"  
  }
  
  fg_df <- gather_coefs(fg_list, function(x) x)  
  
  inla_tab <- pool_family(inla_df, "INLA")
  cox_tab  <- pool_family(cox_df,  "Cox")
  fg_tab   <- pool_family(fg_df,   fg_label)
  
  wide_tbl <- list(inla_tab, cox_tab, fg_tab) %>%
    purrr::compact() %>%
    purrr::reduce(dplyr::full_join, by = "term")
  
  if (is.null(wide_tbl) || !nrow(wide_tbl)) return(tibble::tibble())
  
  p_long <- wide_tbl %>%
    tidyr::pivot_longer(
      cols = tidyselect::ends_with("_p"),
      names_to = "family",
      values_to = "p_raw"
    )
  adj_vec <- stats::p.adjust(p_long$p_raw, method = "holm")
  p_long  <- dplyr::mutate(p_long, p_adj = adj_vec)
  
  p_adj_wide <- p_long %>%
    dplyr::select(term, family, p_adj) %>%
    tidyr::pivot_wider(
      names_from  = family,
      values_from = p_adj,
      names_glue  = "{family}_adj"
    )
  
  out <- wide_tbl %>%
    dplyr::full_join(p_adj_wide, by = "term") %>%
    dplyr::mutate(term = nice_term(term)) %>%
    dplyr::arrange(term)
  
  fmt_p <- function(p) ifelse(is.na(p), "", sprintf("%.2g", p))
  out <- out %>% dplyr::mutate(dplyr::across(tidyselect::ends_with("_adj"), fmt_p))
  
  out
}
