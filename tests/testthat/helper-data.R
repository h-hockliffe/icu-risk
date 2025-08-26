toy_cp <- data.frame(
  patient    = 1:2,
  study_grp  = factor(c("A", "A")),
  start      = 0,
  stop       = c(5, 7),
  status     = c(1, 0),
  age        = c(60, 70),
  bmi        = c(26, 22),
  egfr       = c(75, 90),
  sex        = factor(c("Man", "Vrouw"), levels = c("Man", "Vrouw"))
)

make_cp_toy <- function(n = 120, seed = 1L) {
  set.seed(seed)
  id <- seq_len(n)
  age <- round(rnorm(n, 68, 10))
  bmi <- pmax(15, rnorm(n, 27, 5))
  egfr <- pmin(120, pmax(5, rlnorm(n, log(60), 0.6)))
  sex <- factor(sample(c("Man", "Vrouw"), n, TRUE))
  study_grp <- factor(sample(c("CentreA", "CentreB", "Merged-small-sites"), n, TRUE, prob = c(.2, .28, .52)))
  t <- rexp(n, rate = 1/12)             # ~12-day mean
  cause <- sample(c(1L, 2L, 0L), n, TRUE, prob = c(.12, .65, .23))
  cause[t > 90] <- 0L
  data.frame(
    .imp = 1L,
    patient = paste0("p", id),
    study_grp = study_grp,
    age = age,
    bmi = bmi,
    egfr = egfr,
    sex = sex,
    icu_admit = as.POSIXct("2023-01-01", tz = "UTC"),
    icu_disc  = as.POSIXct("2023-01-01", tz = "UTC") + as.difftime(pmin(90, t), units = "days"),
    start = 0, stop = pmin(90, t),
    status = cause
  )
}

fake_inla_death_fit <- function() {
  rn <- c("log_egfr","log_bmi","sexVrouw","age_spline1","age_spline2","age_spline3")
  sf <- data.frame(
    mean = c(-0.40, 0.10, 0.05, 0.02, -0.01, 0.00),
    sd   = c(0.08,  0.05, 0.06, 0.04,  0.04, 0.04),
    `0.025quant` = c(-0.56, 0.00, -0.07, -0.05, -0.09, -0.08),
    `0.975quant` = c(-0.24, 0.20,  0.17,  0.09,  0.07,  0.08),
    row.names = rn,
    check.names = FALSE
  )
  list(
    model = list(summary.fixed = sf),
    age_mean = 68,
    age_knots_int = c(-11, 0, 7),
    age_knots_bd  = c(-30, 25),
    cuts = c(0, 10, 20, 40, 60, 90)
  )
}

fake_mcmc_samp <- function(n = 2000, seed = 1L) {
  set.seed(seed)
  mat <- cbind(
    b_egfr = rnorm(n, -0.39, 0.08),
    b_bmi  = rnorm(n,  0.11, 0.05),
    b_sex  = rnorm(n,  0.05, 0.06),
    b_a1   = rnorm(n,  0.02, 0.04),
    b_a2   = rnorm(n, -0.01, 0.04),
    b_a3   = rnorm(n,  0.00, 0.04)
  )
  coda::mcmc.list(coda::mcmc(mat))
}
