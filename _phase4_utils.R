# STAT3888 | Gordon | Phase 4
# Shared reproducible helpers used by three R Markdown reports.
# Run each Rmd from its location in the code/ folder (RStudio Knit default).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(knitr)
})
options(scipen = 5, digits = 5)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE,
                      fig.width = 8.5, fig.height = 4.8,
                      fig.align = 'center', dpi = 140)

phase4_root <- function() {
  if (dir.exists('../inputs')) return(normalizePath('..'))
  if (dir.exists('inputs')) return(normalizePath('.'))
  stop('Cannot find inputs/. Knit the Rmd files in the code/ folder of the Phase4_PUFA package.')
}
phase4_dir <- phase4_root()
input_dir <- file.path(phase4_dir, 'inputs')
figure_dir <- file.path(phase4_dir, 'figures')
table_dir <- file.path(phase4_dir, 'tables')
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

locate_file <- function(stem, ext) {
  exact <- file.path(input_dir, paste0(stem, '.', ext))
  if (file.exists(exact)) return(exact)
  f <- list.files(input_dir, full.names = TRUE)
  eligible <- f[grepl(paste0('^', stem, '(\\([0-9]+\\))*\\.', ext, '$'), basename(f), ignore.case = TRUE)]
  if (length(eligible) != 1L) stop('Cannot uniquely locate file: ',stem,'.',ext, ' in inputs/')
  eligible[[1]]
}

phase4_load <- function(reps = FALSE) {
  nn <- as.data.table(readRDS(locate_file('NNPAS_PUFA_cholesterol_study', 'rds')))
  coh <- as.data.table(readRDS(locate_file('phase2_clean_cohort_with_raw_audit', 'rds')))
  if (anyDuplicated(nn$ABSPID) || anyDuplicated(coh$ABSPID)) stop('Duplicate ABSPID in RDS sources')
  want <- c('ABSPID', 'seifa_quintile', 'remoteness', 'FASTSTAD',
            'smoking3', 'mvpa_equiv_min_week', 'alcohol_g_day1')
  stopifnot(all(want %in% names(coh)))
  d <- merge(nn, coh[, ..want], by = 'ABSPID', all.x = TRUE, sort = FALSE)
  stopifnot(nrow(d) == nrow(nn), !anyNA(d$ABSPID))
  repcols <- sprintf('RPWGT%02d', 1:60)
  bio_cols <- c('ABSPID', 'NHMSPERW', if (reps) repcols)
  path <- locate_file('AHSnpa11bb', 'csv')
  bio <- fread(path, select = bio_cols, showProgress = FALSE)
  stopifnot(!anyDuplicated(bio$ABSPID))
  d <- merge(d, bio, by = 'ABSPID', all.x = TRUE, sort = FALSE)
  stopifnot(nrow(d) == nrow(nn), all(!is.na(d$NHMSPERW)))
  if (any(abs(d$NHMSPERW - d$biomedical_weight) > 0.1, na.rm = TRUE))
    stop('Biomedical weight mismatch between Phase 2 and raw file')
  d[, exposure := fifelse(pufa_sfa_group == 'High: 1.0-1.5', 1L,
                     fifelse(pufa_sfa_group == 'Low: <0.45', 0L, NA_integer_))]
  list(data = d, replicate_columns = repcols)
}

phase4_primary <- function(d, repcols = character()) {
  for (nm in c('age', 'biomedical_weight', repcols))
    if (nm %in% names(d)) d[, (nm) := as.numeric(get(nm))]
  primary_covars <- c('age', 'sex', 'seifa_quintile', 'remoteness')
  needs <- c('exposure', 'total_chol_elevated', 'biomedical_weight', primary_covars, repcols)
  x <- copy(d[!is.na(exposure) & complete.cases(d[, ..needs]) & biomedical_weight > 0])
  x[, `:=`(sex = droplevels(factor(sex)),
           seifa_quintile = droplevels(factor(seifa_quintile)),
           remoteness = droplevels(factor(remoteness)))]
  if (any(table(x$exposure) < 20)) stop('At least one exposure group has fewer than 20 observations')
  x
}

wmean <- function(x, wt) {
  ok <- is.finite(x) & is.finite(wt) & wt > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], wt[ok])
}
ESS <- function(w) sum(w)^2 / sum(w^2)

phase4_estimate <- function(x, wt, covars, outcome = 'total_chol_elevated', clip = 0) {
  # x must be complete for the requested outcome and covariates.
  stopifnot(nrow(x) == length(wt), all(is.finite(wt)), all(wt >= 0), sum(wt) > 0)
  z <- copy(x)
  ww <- wt / (sum(wt) / sum(wt > 0))
  z[, .model_weight := ww]
  fprop <- reformulate(covars, response = 'exposure')
  fy <- reformulate(c('exposure', covars), response = outcome)
  psfit <- suppressWarnings(glm(fprop, data = z, weights = .model_weight,
                                family = quasibinomial(), control = glm.control(maxit = 100)))
  yfit <- suppressWarnings(glm(fy, data = z, weights = .model_weight,
                               family = quasibinomial(), control = glm.control(maxit = 100)))
  if (!psfit$converged || !yfit$converged) stop('GLM did not converge')
  ps0 <- as.numeric(predict(psfit, newdata = z, type = 'response'))
  if (any(!is.finite(ps0)) || any(ps0 <= 0 | ps0 >= 1)) stop('Invalid/extreme propensity score')
  ps <- if (clip > 0) pmin(1 - clip, pmax(clip, ps0)) else ps0
  zlow <- copy(z); zlow[, exposure := 0L]
  zhigh <- copy(z); zhigh[, exposure := 1L]
  m0 <- as.numeric(predict(yfit, newdata = zlow, type = 'response'))
  m1 <- as.numeric(predict(yfit, newdata = zhigh, type = 'response'))
  aa <- z$exposure; yy <- z[[outcome]]
  w1 <- wt * aa / ps
  w0 <- wt * (1 - aa) / (1 - ps)
  estimates <- rbindlist(list(
    data.table(method = 'G-computation', low = wmean(m0,wt), high = wmean(m1,wt)),
    data.table(method = 'IPW', low = wmean(yy,w0), high = wmean(yy,w1)),
    data.table(method = 'AIPW',
               low = wmean(m0 + (1-aa)/(1-ps)*(yy-m0), wt),
               high = wmean(m1 + aa/ps*(yy-m1), wt))
  ))
  estimates[, `:=`(RD = high - low, RR = high / low)]
  list(estimates = estimates, ps = ps0, ESS_high = ESS(w1),
       ESS_low = ESS(w0), psfit = psfit, outcome_fit = yfit)
}

phase4_replicate_rd <- function(x, repcols, covars) {
  main <- phase4_estimate(x, x$biomedical_weight, covars)$estimates
  rr <- matrix(NA_real_, nrow = length(repcols), ncol = nrow(main))
  failures <- character()
  for (j in seq_along(repcols)) {
    result <- tryCatch(phase4_estimate(x, as.numeric(x[[repcols[j]]]), covars), error = function(e) e)
    if (inherits(result, 'error')) failures <- c(failures, paste(repcols[j], conditionMessage(result)))
    else rr[j, ] <- result$estimates$RD
  }
  if (length(failures)) {
    warning('Replication failed: ', paste(failures, collapse = '; '))
    main[, `:=`(SE_RD = NA_real_, CI_low = NA_real_, CI_high = NA_real_)]
  } else {
    se <- sqrt((59/60) * colSums(sweep(rr, 2, main$RD, '-')^2))
    crit <- qt(.975, df = 59)
    main[, `:=`(SE_RD = se, CI_low = RD - crit * se, CI_high = RD + crit * se)]
  }
  attr(main, 'failed_replicates') <- failures
  main
}

save_table <- function(x, name) fwrite(as.data.table(x), file.path(table_dir, name), na = '')
save_figure <- function(p, name, width = 8.5, height = 5) {
  ggsave(file.path(figure_dir, name), p, width = width, height = height, dpi = 170)
}
format_table <- function(x, caption) knitr::kable(as.data.frame(x), caption = caption, digits = 4)
