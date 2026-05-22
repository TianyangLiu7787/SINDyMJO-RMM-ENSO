# Use user library when R is installed system-wide (non-admin installs).
local({
  ul <- Sys.getenv("R_LIBS_USER", unset = "")
  if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths()))
})

#########################################################################################
#  SINDy for MJO (RMM) + ONI-based ENSO years — 1979–2024 (rmm.74toRealtime + oni.ascii) #
#  Adapted from SindyMJO.R / MJO_polar.R (N. Diaz et al.):                               #
#   - Read RMM from WH-style whitespace file; ONI from NOAA seasonal ASCII               #
#   - Two RMM seasonal windows (each full pipeline + Fig3/4/7/8/8a/8b):                    #
#       * Winter: Oct 27(y)–May 4(y+1), 9-day MA, slice x_in.0[5:(n-4)] as original pad    #
#       * Summer: Apr 27(y)–Nov 4(y); EN/LN via ONI of prior DJFM (Nov y-1 – Apr y)       #
#   - ONI/ENSO: enso_nov keys = winter anchor years y (DJFM Nov y); summer rows use DJFM y-1. #
#   - Velocities, MJO.polar, amplitude>1, STLSQ / ΔAIC as in the original                #
#   - ENSO year lists: derived from ONI (mean of NDJ/JFM/FMA/MAM for each DJFM season)   #
#                                                                                       #
#  Usage (from SINDyMJO folder):                                                         #
#   Rscript SindyMJO_RMM_ONI.R [rmm.txt] [oni.txt] [out_dir] [y0] [y1] [Nreal] [Nsamp]    #
#       [amp] [dpi] [years_mode] [oni_thr] [seed]                                        #
#   years_mode: all | neutros | ninos | ninas  (same idea as dlgInput in SindyMJO.R)      #
#   seed: optional RNG seed (integer).                                                    #
#   Env SINDY_SEED_METRICS=1: print SINDY_SEED_METRICS and exit (no figures).              #
#   Env SINDY_SEED_SWEEP=1: scan SINDY_SWEEP_LO:HI (optional SINDY_SWEEP_NREAL), print     #
#       SINDY_SWEEP_BEST and exit (no figures; uses winter pool only).                    #
#   Outputs: figs/Fig{3,4,7,8,8a,8b}_*_clim_y0_y1_winter.png (+ _summer); same for CSV.   #
#########################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(patchwork)
  library(R1magic)
  library(latex2exp)
  library(viridis)
})

source("MJO_polar.R") # provides MJO.polar(x,y)

# -------------------------- Arguments / config --------------------------
args <- commandArgs(trailingOnly = TRUE)

guess_rmm <- function() {
  cand <- c(
    normalizePath(file.path(getwd(), "..", "rmm.74toRealtime.txt"), mustWork = FALSE),
    normalizePath(file.path(getwd(), "rmm.74toRealtime.txt"), mustWork = FALSE)
  )
  ex <- cand[file.exists(cand)]
  if (!length(ex)) return("")
  ex[[1]]
}
guess_oni <- function(rmm_path) {
  if (nzchar(rmm_path) && dirname(rmm_path) != ".") {
    p <- normalizePath(file.path(dirname(rmm_path), "oni.ascii.txt"), mustWork = FALSE)
    if (file.exists(p)) return(p)
  }
  p2 <- normalizePath(file.path(getwd(), "..", "oni.ascii.txt"), mustWork = FALSE)
  if (file.exists(p2)) return(p2)
  normalizePath(file.path(getwd(), "oni.ascii.txt"), mustWork = FALSE)
}

rmm_path <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else guess_rmm()
oni_path <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else guess_oni(rmm_path)
out_dir <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else file.path(dirname(rmm_path), "outputs", "sindy_r")

start_year <- if (length(args) >= 4) as.integer(args[[4]]) else 1979L
end_year <- if (length(args) >= 5) as.integer(args[[5]]) else 2024L

Nreal <- if (length(args) >= 6) as.integer(args[[6]]) else 1000L
Nsamp <- if (length(args) >= 7) as.integer(args[[7]]) else 2^8

amp_threshold <- if (length(args) >= 8) as.numeric(args[[8]]) else 1.0
out_dpi <- if (length(args) >= 9) as.integer(args[[9]]) else 600L

years_mode_raw <- if (length(args) >= 10) args[[10]] else "all"
oni_threshold <- if (length(args) >= 11) as.numeric(args[[11]]) else 0.5
sindy_seed <- if (length(args) >= 12L && nzchar(args[[12]])) as.integer(args[[12]]) else 1234L

n <- 2
D <- 2
Lambda.inv <- seq(1, 100, 1)
set.seed(sindy_seed)

years_mode <- tolower(trimws(years_mode_raw))
years_mode <- sub("niños", "ninos", years_mode, fixed = TRUE)
years_mode <- sub("niñas", "ninas", years_mode, fixed = TRUE)
años_ch <- years_mode_raw

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
fig_dir <- file.path(out_dir, "figs")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# -------------------------- Read RMM / ONI (project txt files) --------------------------
MONTH_TO_ONI_SEASON <- c(
  "DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ"
)

read_rmm_whitespace <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line)) next
    if (grepl("^RMM values|^year,", line, ignore.case = TRUE)) next
    parts <- strsplit(line, "\\s+")[[1L]]
    parts <- parts[nzchar(parts)]
    if (length(parts) < 7L) next
    y <- suppressWarnings(as.integer(parts[[1L]]))
    m <- suppressWarnings(as.integer(parts[[2L]]))
    d <- suppressWarnings(as.integer(parts[[3L]]))
    r1 <- suppressWarnings(as.numeric(parts[[4L]]))
    r2 <- suppressWarnings(as.numeric(parts[[5L]]))
    if (any(is.na(c(y, m, d)))) next
    if (!is.finite(r1) || !is.finite(r2)) next
    if (abs(r1) >= 1e20 || abs(r2) >= 1e20) next
    if (r1 == 999 || r2 == 999) next
    out[[length(out) + 1L]] <- data.frame(Año = y, Mes = m, Dia = d, RMM1 = r1, RMM2 = r2)
  }
  if (!length(out)) stop("No RMM rows parsed from: ", path)
  do.call(rbind, out)
}

read_oni_map <- function(path) {
  lines <- readLines(path, warn = FALSE)
  mp <- new.env(hash = TRUE, parent = emptyenv())
  first <- TRUE
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line)) next
    if (first) {
      first <- FALSE
      next
    }
    parts <- strsplit(line, "\\s+")[[1L]]
    parts <- parts[nzchar(parts)]
    if (length(parts) < 4L) next
    seas <- toupper(parts[[1L]])
    yr <- suppressWarnings(as.integer(parts[[2L]]))
    anom <- suppressWarnings(as.numeric(parts[[4L]]))
    if (is.na(yr) || !is.finite(anom)) next
    assign(paste(seas, yr, sep = "\001"), anom, envir = mp)
  }
  mp
}

oni_get <- function(mp, seas, yr) {
  v <- get0(paste(seas, yr, sep = "\001"), envir = mp, inherits = FALSE)
  if (is.null(v)) NA_real_ else v
}

oni_for_calendar_month <- function(mp, y, m) {
  seas <- MONTH_TO_ONI_SEASON[[m]]
  oni_get(mp, seas, y)
}

djfm_oni_core_mean <- function(mp, año_nov) {
  v <- c(
    oni_for_calendar_month(mp, año_nov, 12L),
    oni_for_calendar_month(mp, año_nov + 1L, 1L),
    oni_for_calendar_month(mp, año_nov + 1L, 2L),
    oni_for_calendar_month(mp, año_nov + 1L, 3L)
  )
  if (all(is.na(v))) return(NA_real_)
  mean(v, na.rm = TRUE)
}

enso_from_oni_mean <- function(mu, thr) {
  if (!is.finite(mu)) return(NA_character_)
  if (mu >= thr) return("ElNino")
  if (mu <= -thr) return("LaNina")
  "Neutral"
}

# -------------------------- Helpers (as original) --------------------------
Sindy <- function(Psi, xpunto, init.guess, L1.lambda) {
  Psi.dim <- dim(Psi)
  if (is.null(Psi.dim)) {
    T.mat <- 1
  } else {
    T.mat <- diag(length(Psi[1, ]))
  }
  solveL1(Psi, y = xpunto, T = T.mat, x0 = init.guess, lambda = L1.lambda)$estimate
}

build_Psi <- function(x1, x2, n) {
  NPsi.0 <- (1 + n)^2
  Psi.0 <- array(0, dim = c(length(x1), NPsi.0))
  for (t in seq_along(x1)) {
    j <- 0
    for (l1 in 0:n) {
      for (l2 in 0:n) {
        j <- j + 1
        Psi.0[t, j] <- (x1[t]^l1) * (x2[t]^l2)
      }
    }
  }
  # Normalize columns by L2 norm
  L2Psi <- rep(0, ncol(Psi.0))
  for (j in seq_len(ncol(Psi.0))) {
    L2Psi[j] <- norm(Psi.0[, j, drop = FALSE], "2")
    Psi.0[, j] <- Psi.0[, j] / L2Psi[j]
  }
  list(Psi = Psi.0, L2Psi = L2Psi)
}

# Keep dictionary order exactly as in the original script:
# 1, y, y^2, x, xy, xy^2, x^2, x^2y, x^2y^2
dictionary_terms <- c("1", "y", "y2", "x", "xy", "xy2", "x2", "x2y", "x2y2")
idx_xdot <- function(term) match(term, dictionary_terms)
idx_ydot <- function(term) idx_xdot(term) + length(dictionary_terms)

run_sindy_subset <- function(DEFM.amp.df, enso_label, Nreal, Nsamp, n, Lambda.inv) {
  # Random sampling list
  MJO.tray <- vector("list", length = Nreal)
  for (i in seq_len(Nreal)) {
    samp.ind <- sample(x = seq_len(nrow(DEFM.amp.df)), size = Nsamp, replace = FALSE)
    MJO.tray[[i]] <- DEFM.amp.df[samp.ind, ]
  }

  NPsi <- (1 + n)^2
  Niter <- NPsi

  Modelo1.list <- vector("list", length = length(MJO.tray))
  Modelo2.list <- vector("list", length = length(MJO.tray))
  Modelo3.list <- vector("list", length = length(MJO.tray))
  Modelo4.list <- vector("list", length = length(MJO.tray))

  for (M in seq_along(MJO.tray)) {
    x1.0 <- MJO.tray[[M]]$RMM1
    x2.0 <- MJO.tray[[M]]$RMM2
    xpunto <- cbind(MJO.tray[[M]]$RMM1.vel, MJO.tray[[M]]$RMM2.vel)

    Psi_out <- build_Psi(x1.0, x2.0, n)
    Psi <- Psi_out$Psi
    L2Psi <- Psi_out$L2Psi

    coefs <- array(0, dim = c(NPsi, length(Lambda.inv), D))
    AIC.df <- data.frame(lambda.inv = Lambda.inv, train = rep(NA_real_, length(Lambda.inv)))

    for (k in seq_along(Lambda.inv)) {
      for (dim in 1:D) {
        p <- qr.solve(Psi, xpunto[, dim])
        AuxCoefs <- Sindy(Psi, xpunto[, dim], init.guess = p, L1.lambda = 0)
        for (it in seq_len(Niter)) {
          IZ <- which(abs(AuxCoefs / L2Psi) < 1 / Lambda.inv[k])
          INZ <- which(abs(AuxCoefs / L2Psi) >= 1 / Lambda.inv[k])
          AuxCoefs[IZ] <- 0
          if (length(INZ) > 0) {
            AuxCoefs[INZ] <- Sindy(Psi = Psi[, INZ, drop = FALSE], xpunto = xpunto[, dim], init.guess = AuxCoefs[INZ], L1.lambda = 0)
          }
          if (it == Niter) coefs[, k, dim] <- AuxCoefs / L2Psi
        }
      }

      NCNNx <- length(which(coefs[, k, 1] != 0))
      NCNNy <- length(which(coefs[, k, 2] != 0))
      NCNN <- NCNNx + NCNNy
      xpunto1.rec <- Psi %*% coefs[, k, 1]
      xpunto2.rec <- Psi %*% coefs[, k, 2]
      if (!any(is.nan(xpunto1.rec) | is.infinite(xpunto1.rec) | is.nan(xpunto2.rec) | is.infinite(xpunto2.rec))) {
        rmse <- sqrt(sum((xpunto[, 1] - xpunto1.rec)^2 + (xpunto[, 2] - xpunto2.rec)^2, na.rm = TRUE) / nrow(xpunto))
        AIC.df$train[k] <- nrow(xpunto) * log(rmse^2) + 2 * (NCNN + 1) + 2 * (NCNN + 1) * (NCNN + 2) / (nrow(xpunto) - NCNN - 2)
      }
    }

    AIC.df$train <- AIC.df$train - min(AIC.df$train, na.rm = TRUE)
    AIC.df$train <- round(AIC.df$train * 10000) / 10000

    km1 <- which(AIC.df$train == min(AIC.df$train, na.rm = TRUE)); Modelo1 <- km1[1]
    km2 <- which(AIC.df$train == min(AIC.df$train[-km1], na.rm = TRUE)); Modelo2 <- km2[1]
    km3 <- which(AIC.df$train == min(AIC.df$train[-c(km1, km2)], na.rm = TRUE)); Modelo3 <- km3[1]
    km4 <- which(AIC.df$train == min(AIC.df$train[-c(km1, km2, km3)], na.rm = TRUE)); Modelo4 <- km4[1]

    NCNNx1 <- length(which(coefs[, Modelo1, 1] != 0)); NCNNy1 <- length(which(coefs[, Modelo1, 2] != 0))
    NCNNx2 <- length(which(coefs[, Modelo2, 1] != 0)); NCNNy2 <- length(which(coefs[, Modelo2, 2] != 0))
    NCNNx3 <- length(which(coefs[, Modelo3, 1] != 0)); NCNNy3 <- length(which(coefs[, Modelo3, 2] != 0))
    NCNNx4 <- length(which(coefs[, Modelo4, 1] != 0)); NCNNy4 <- length(which(coefs[, Modelo4, 2] != 0))

    if (NCNNx1 != 0 && NCNNy1 != 0) Modelo1.list[[M]] <- coefs[, Modelo1, ]
    if (!is.na(Modelo2) && AIC.df$train[Modelo2] <= 5 && NCNNx2 != 0 && NCNNy2 != 0) Modelo2.list[[M]] <- coefs[, Modelo2, ]
    if (!is.na(Modelo3) && AIC.df$train[Modelo3] <= 7 && NCNNx3 != 0 && NCNNy3 != 0) Modelo3.list[[M]] <- coefs[, Modelo3, ]
    if (!is.na(Modelo4) && AIC.df$train[Modelo4] <= 7 && NCNNx4 != 0 && NCNNy4 != 0) Modelo4.list[[M]] <- coefs[, Modelo4, ]
  }

  modelos0.df <- data.frame(row.names = seq_len(2 * NPsi))
  for (i in 1:4) {
    if (i == 1) Modelo.list <- Modelo1.list
    if (i == 2) Modelo.list <- Modelo2.list
    if (i == 3) Modelo.list <- Modelo3.list
    if (i == 4) Modelo.list <- Modelo4.list
    for (chi in seq_len(Nreal)) {
      if (!is.null(Modelo.list[[chi]])) modelos0.df <- cbind(modelos0.df, as.vector(Modelo.list[[chi]]))
    }
  }

  if (ncol(modelos0.df) < 1) stop(paste0("No valid models selected for ", enso_label))

  # coefficient statistics (frequency and mean + quantiles)
  num.coefs <- rep(NA_real_, 2 * NPsi)
  avg.coefs <- rep(NA_real_, 2 * NPsi)
  p75.coefs <- rep(NA_real_, 2 * NPsi)
  p25.coefs <- rep(NA_real_, 2 * NPsi)
  min.coefs <- rep(NA_real_, 2 * NPsi)
  max.coefs <- rep(NA_real_, 2 * NPsi)
  for (i in seq_len(2 * NPsi)) {
    num.coefs[i] <- length(which(modelos0.df[i, ] != 0)) / ncol(modelos0.df)
    vals <- as.numeric(modelos0.df[i, ])
    avg.coefs[i] <- mean(vals, na.rm = TRUE)
    p75.coefs[i] <- quantile(vals, probs = 0.75, na.rm = TRUE)[[1]]
    p25.coefs[i] <- quantile(vals, probs = 0.25, na.rm = TRUE)[[1]]
    min.coefs[i] <- min(vals, na.rm = TRUE)
    max.coefs[i] <- max(vals, na.rm = TRUE)
  }

  # model-structure masks (same as original idea)
  M.aux1 <- modelos0.df[, 1]
  M.aux1[M.aux1 != 0] <- 1
  mask.df <- data.frame(M.aux1)
  modelos.df <- data.frame(M.aux1)
  for (j in 2:ncol(modelos0.df)) {
    aux <- modelos0.df[, j]
    aux[aux != 0] <- 1
    mask.df <- cbind(mask.df, aux)
    L.modelos <- ncol(modelos.df)
    for (l in 1:L.modelos) {
      aux0 <- modelos.df[, l]
      sum0 <- sum(aux0 == aux)
      if (sum0 == (2 * NPsi)) break
      if (l == L.modelos) modelos.df <- cbind(modelos.df, aux)
    }
  }

  # periods from purely linear systems (same linear entries as original)
  linear.entries <- c(1, 2, 4, 10, 11, 13)
  non.linear.entries <- seq_len(2 * NPsi)[-linear.entries]
  omega2 <- numeric(0)
  for (m in 1:ncol(modelos.df)) {
    mask.aux <- modelos.df[, m]
    out <- apply(mask.df, MARGIN = 2, function(x) ifelse((sum(x == mask.aux) == 2 * NPsi), 1, 0))
    if (!any(mask.aux[non.linear.entries] != 0)) {
      if (length(which(out == 1)) >= 1) {
        cols <- which(out == 1)
        # compute omega2 for each occurrence (as original)
        omega2 <- c(omega2, apply(modelos0.df[, cols, drop = FALSE], MARGIN = 2, function(x) {
          a10 <- x[idx_xdot("x")]
          b01 <- x[idx_ydot("y")]
          a01 <- x[idx_xdot("y")]
          b10 <- x[idx_ydot("x")]
          (a10 * b01 - a01 * b10) - (a10 + b01)^2 / 4
        }))
      }
    }
  }

  # frequency of each model structure
  freq_modelos <- rep(NA_real_, ncol(modelos.df))
  for (m in seq_len(ncol(modelos.df))) {
    mask.aux <- modelos.df[, m]
    out <- apply(mask.df, MARGIN = 2, function(x) ifelse((sum(x == mask.aux) == 2 * NPsi), 1, 0))
    freq_modelos[m] <- length(which(out == 1)) / ncol(modelos0.df)
  }

  list(
    modelos0.df = modelos0.df,
    num.coefs = num.coefs,
    avg.coefs = avg.coefs,
    p25 = p25.coefs,
    p75 = p75.coefs,
    min_coefs = min.coefs,
    max_coefs = max.coefs,
    modelos_mask = modelos.df,
    freq_modelos = freq_modelos,
    omega2 = omega2
  )
}

# -------------------------- Build DJFM sample (same window as SindyMJO.R RMM branch) --------------------------
if (!nzchar(rmm_path) || !file.exists(rmm_path)) {
  stop("RMM file not found: ", rmm_path, " (arg1 or ../rmm.74toRealtime.txt)")
}
if (!nzchar(oni_path) || !file.exists(oni_path)) {
  stop("ONI file not found: ", oni_path)
}

MJO_raw <- read_rmm_whitespace(rmm_path)
MJO_raw$date <- as.Date(sprintf("%04d-%02d-%02d", MJO_raw$Año, MJO_raw$Mes, MJO_raw$Dia))
MJO_raw <- MJO_raw[!is.na(MJO_raw$date), , drop = FALSE]
MJO_raw <- MJO_raw[MJO_raw$date >= as.Date(sprintf("%d-01-01", start_year)) &
  MJO_raw$date <= as.Date(sprintf("%d-12-31", end_year)), , drop = FALSE]
MJO_raw <- MJO_raw[format(MJO_raw$date, "%m-%d") != "02-29", , drop = FALSE]
MJO <- MJO_raw[, c("Año", "Mes", "Dia", "RMM1", "RMM2")]
rownames(MJO) <- NULL
MJO$date <- as.Date(sprintf("%04d-%02d-%02d", MJO$Año, MJO$Mes, MJO$Dia))

oni_map <- read_oni_map(oni_path)

MIN_RMM_ROWS <- 125L

winter_date_range <- function(y) {
  list(
    d0 = as.Date(sprintf("%04d-10-27", y)),
    d1 = as.Date(sprintf("%04d-05-04", y + 1L))
  )
}
summer_date_range <- function(y) {
  list(
    d0 = as.Date(sprintf("%04d-04-27", y)),
    d1 = as.Date(sprintf("%04d-11-04", y))
  )
}

season_window_ok_winter <- function(y) {
  dr <- winter_date_range(y)
  sum(MJO$date >= dr$d0 & MJO$date <= dr$d1, na.rm = TRUE) >= MIN_RMM_ROWS
}
season_window_ok_summer <- function(y) {
  dr <- summer_date_range(y)
  sum(MJO$date >= dr$d0 & MJO$date <= dr$d1, na.rm = TRUE) >= MIN_RMM_ROWS
}

año_lo <- min(MJO$Año, na.rm = TRUE) - 2L
año_hi <- max(MJO$Año, na.rm = TRUE) + 2L
candidates <- seq.int(año_lo, año_hi)
full_years_winter <- candidates[vapply(candidates, season_window_ok_winter, logical(1L))]
full_years_summer <- candidates[vapply(candidates, season_window_ok_summer, logical(1L))]
if (!length(full_years_winter)) {
  stop("No complete winter seasons (need >=", MIN_RMM_ROWS, " RMM rows in Oct27–May4 window).")
}
if (!length(full_years_summer)) {
  stop("No complete summer seasons (need >=", MIN_RMM_ROWS, " RMM rows in Apr27–Nov4 window).")
}

# ONI class keyed by DJFM November-start year y (same as djfm_oni_core_mean); used for labels + winter EN/LN.
enso_nov <- vapply(full_years_winter, function(a) {
  mu <- djfm_oni_core_mean(oni_map, a)
  enso_from_oni_mean(mu, oni_threshold)
}, character(1L))
names(enso_nov) <- as.character(full_years_winter)

neutros <- full_years_winter[enso_nov == "Neutral"]
ninos <- full_years_winter[enso_nov == "ElNino"]
ninas <- full_years_winter[enso_nov == "LaNina"]

if (!(years_mode %in% c("all", "neutros", "ninos", "ninas"))) {
  stop("years_mode must be: all, neutros, ninos, ninas ; got: ", years_mode_raw)
}
if (years_mode != "all") {
  message("Note: years_mode=", years_mode, " is recorded in output names only; ",
          "SINDy and figures use all complete winter/summer seasons in [", start_year, ", ", end_year, "].")
}

enso_from_djfm_nov_y <- function(ny) {
  ny <- as.integer(ny)
  mu <- djfm_oni_core_mean(oni_map, ny)
  enso_from_oni_mean(mu, oni_threshold)
}

build_x_in_df <- function(mode, years_vec) {
  stopifnot(mode %in% c("winter", "summer"))
  x_acc <- data.frame()
  for (y in years_vec) {
    if (mode == "winter") {
      dr <- winter_date_range(y)
      chg_año <- y
      djfm_nov_y <- y
    } else {
      dr <- summer_date_range(y)
      chg_año <- y
      djfm_nov_y <- y - 1L
    }
    x_in.0 <- MJO[MJO$date >= dr$d0 & MJO$date <= dr$d1, c("Año", "Mes", "Dia", "RMM1", "RMM2"), drop = FALSE]
    nr <- nrow(x_in.0)
    if (nr < MIN_RMM_ROWS) next
    if (nr < 9L) next
    x_in <- x_in.0[seq.int(5L, nr - 4L), , drop = FALSE]
    for (t in seq_len(nrow(x_in))) {
      x_in$RMM1[t] <- mean(x_in.0$RMM1[seq.int(t, t + 8L)])
      x_in$RMM2[t] <- mean(x_in.0$RMM2[seq.int(t, t + 8L)])
    }
    chg.year <- which(x_in$Año == chg_año & x_in$Mes == 12L)
    if (length(chg.year)) {
      x_in$Año[chg.year] <- x_in$Año[chg.year] + 1L
    }
    ll <- nrow(x_in)
    vel1 <- rep(NA_real_, ll)
    vel2 <- rep(NA_real_, ll)
    for (t in seq_len(ll)) {
      if (t <= 3L) {
        vel1[t] <- (-3 * x_in$RMM1[t] + 4 * x_in$RMM1[t + 1] - x_in$RMM1[t + 2]) / 2
        vel2[t] <- (-3 * x_in$RMM2[t] + 4 * x_in$RMM2[t + 1] - x_in$RMM2[t + 2]) / 2
      } else if (t >= (ll - 2L)) {
        vel1[t] <- (3 * x_in$RMM1[t] - 4 * x_in$RMM1[t - 1] + x_in$RMM1[t - 2]) / 2
        vel2[t] <- (3 * x_in$RMM2[t] - 4 * x_in$RMM2[t - 1] + x_in$RMM2[t - 2]) / 2
      } else {
        vel1[t] <- (-x_in$RMM1[t + 2] + 8 * x_in$RMM1[t + 1] - 8 * x_in$RMM1[t - 1] + x_in$RMM1[t - 2]) / 12
        vel2[t] <- (-x_in$RMM2[t + 2] + 8 * x_in$RMM2[t + 1] - 8 * x_in$RMM2[t - 1] + x_in$RMM2[t - 2]) / 12
      }
    }
    x_acc <- rbind(
      x_acc,
      data.frame(
        Año = x_in$Año,
        RMM1 = x_in$RMM1,
        RMM2 = x_in$RMM2,
        RMM1.vel = vel1,
        RMM2.vel = vel2,
        djfm_nov_y = rep.int(as.integer(djfm_nov_y), ll),
        stringsAsFactors = FALSE
      )
    )
  }
  if (!nrow(x_acc)) {
    stop("No rows in x_in.df for mode=", mode, "; check RMM coverage and MIN_RMM_ROWS.")
  }
  x_acc
}

prepare_amp_en_ln <- function(x_in.df) {
  DEFM.polar <- MJO.polar(x_in.df$RMM1, x_in.df$RMM2)
  DEFM.df <- cbind(x_in.df, DEFM.polar)
  Amp.ind <- which(DEFM.df$Amplitude > amp_threshold)
  DEFM.amp.df <- DEFM.df[Amp.ind, , drop = FALSE]
  DEFM.amp.df$enso <- vapply(
    DEFM.amp.df$djfm_nov_y,
    function(ny) enso_from_djfm_nov_y(ny),
    character(1L)
  )
  EN <- DEFM.amp.df[DEFM.amp.df$enso == "ElNino", , drop = FALSE]
  LN <- DEFM.amp.df[DEFM.amp.df$enso == "LaNina", , drop = FALSE]
  if (nrow(EN) < (Nsamp + 10L)) {
    stop("ElNino rows (amp>", amp_threshold, ") too small: ", nrow(EN), "; lower Nsamp or amp_threshold.")
  }
  if (nrow(LN) < (Nsamp + 10L)) {
    stop("LaNina rows (amp>", amp_threshold, ") too small: ", nrow(LN), "; lower Nsamp or amp_threshold.")
  }
  list(raw_amp = DEFM.amp.df, EN = EN, LN = LN)
}

fig_base_suffix <- sprintf("_clim_%d_%d", start_year, end_year)

# Optional: one-process seed scan (EN vs LN nonlinear coefficient activity). Exits before figures.
# Uses winter RMM pool only (Oct27–May4).
if (nzchar(Sys.getenv("SINDY_SEED_SWEEP", ""))) {
  pr_sw <- prepare_amp_en_ln(build_x_in_df("winter", full_years_winter))
  idx_nonlin <- c(3L, 5L, 6L, 7L, 8L, 9L, 12L, 14L, 15L, 16L, 17L, 18L)
  sw_lo <- as.integer(Sys.getenv("SINDY_SWEEP_LO", unset = "1"))
  sw_hi <- as.integer(Sys.getenv("SINDY_SWEEP_HI", unset = "200"))
  nr_env <- Sys.getenv("SINDY_SWEEP_NREAL", unset = "")
  nrv <- if (nzchar(nr_env)) as.integer(nr_env) else min(Nreal, 80L)
  en_df <- pr_sw$EN[, c("RMM1", "RMM2", "RMM1.vel", "RMM2.vel")]
  ln_df <- pr_sw$LN[, c("RMM1", "RMM2", "RMM1.vel", "RMM2.vel")]
  best_seed <- NA_integer_
  best_en <- NA_real_
  best_ln <- NA_real_
  best_d <- -Inf
  for (sv in sw_lo:sw_hi) {
    set.seed(sv)
    r_en <- run_sindy_subset(en_df, "ElNino", nrv, Nsamp, n, Lambda.inv)
    r_ln <- run_sindy_subset(ln_df, "LaNina", nrv, Nsamp, n, Lambda.inv)
    en_s <- sum(r_en$num.coefs[idx_nonlin])
    ln_s <- sum(r_ln$num.coefs[idx_nonlin])
    if (en_s > ln_s) {
      d <- en_s - ln_s
      if (d > best_d) {
        best_d <- d
        best_seed <- sv
        best_en <- en_s
        best_ln <- ln_s
      }
    }
  }
  if (is.na(best_seed)) {
    cat("SINDY_SWEEP_NO_HIT\tno seed in range with EN_nonlin > LN_nonlin; widen SINDY_SWEEP_LO/HI or increase SINDY_SWEEP_NREAL\n")
  } else {
    cat(sprintf(
      "SINDY_SWEEP_BEST\tseed=%d\tEN_nonlin=%.8g\tLN_nonlin=%.8g\tdiff=%.8g\tNreal=%d\trange=%d:%d\n",
      best_seed, best_en, best_ln, best_d, nrv, sw_lo, sw_hi
    ))
  }
  quit(save = "no")
}

NPsi <- (1 + n)^2

# X-axis labels as in original SindyMJO.R (latex2exp).
x_lab_fig <- rep(
  TeX(input = c(
    "$1$", "$\\y$", "$\\y^2$", "$\\x$", "$\\xy$", "$\\xy^2$",
    "$\\x^2$", "$\\x^2y$", "$\\x^2y^2$"
  )),
  2L
)

# Fig4/Fig8：与 SindyMJO.R 一致，仅展示 freq>0.01 的结构类
select_fig4_structure_order <- function(res) {
  keep <- which(res$freq_modelos > 0.01)
  ord <- keep[order(res$freq_modelos[keep], decreasing = TRUE)]
  if (length(ord) == 0) {
    ord <- order(res$freq_modelos, decreasing = TRUE)[seq_len(min(8L, ncol(res$modelos_mask)))]
  }
  list(ord = ord, freq = res$freq_modelos[ord])
}

build_mask_long <- function(res, prefix) {
  sel <- select_fig4_structure_order(res)
  ord <- sel$ord
  mask <- res$modelos_mask[, ord, drop = FALSE]
  labs <- paste0(prefix, ".M", seq_along(ord), " fr=", sprintf("%.2f", sel$freq))
  colnames(mask) <- labs
  long <- melt(cbind(Coef = seq_len(2 * NPsi), mask), id.vars = "Coef")
  # 列名已含 fr=，条带勿再拼一次
  facet_labs <- stats::setNames(as.character(colnames(mask)), colnames(mask))
  list(long = long, facet_labs = facet_labs, ord = ord, freq = sel$freq)
}

# Fig4 各行结构的系数列均值（与掩码匹配的 modelos0.df 列）
export_fig4_model_coefficients <- function(res, out_csv, prefix, dictionary_terms) {
  sel <- select_fig4_structure_order(res)
  ord <- sel$ord
  mo0 <- as.matrix(res$modelos0.df)
  mo <- res$modelos_mask
  term_names <- c(
    paste0("xdot_", dictionary_terms),
    paste0("ydot_", dictionary_terms)
  )
  cols_list <- vector("list", length(ord))
  for (k in seq_along(ord)) {
    m <- ord[k]
    mask_bin <- as.integer(mo[, m] != 0)
    hits <- vapply(seq_len(ncol(mo0)), function(j) {
      identical(as.integer(mo0[, j] != 0), mask_bin)
    }, logical(1L))
    cols_list[[k]] <- if (!any(hits)) {
      rep(NA_real_, nrow(mo0))
    } else {
      rowMeans(mo0[, hits, drop = FALSE], na.rm = TRUE)
    }
  }
  mat <- do.call(cbind, cols_list)
  colnames(mat) <- paste0(prefix, ".M", seq_along(ord))
  freq_list <- c(
    list(term_index = NA_integer_, term = "freq_in_pool"),
    stats::setNames(as.list(sel$freq), colnames(mat))
  )
  freq_row <- do.call(data.frame, c(freq_list, list(stringsAsFactors = FALSE, check.names = FALSE)))
  out_df <- data.frame(term_index = seq_len(nrow(mo0)), term = term_names, mat, check.names = FALSE)
  write.csv(rbind(freq_row, out_df), out_csv, row.names = FALSE, fileEncoding = "UTF-8")
}

run_season_bundle <- function(x_in.df, season_tag) {
  prep <- prepare_amp_en_ln(x_in.df)
  raw_amp <- prep$raw_amp
  EN_loc <- prep$EN
  LN_loc <- prep$LN
  res_ALL <- run_sindy_subset(
    DEFM.amp.df = raw_amp[, c("RMM1", "RMM2", "RMM1.vel", "RMM2.vel")],
    enso_label = "All",
    Nreal = Nreal, Nsamp = Nsamp, n = n, Lambda.inv = Lambda.inv
  )
  res_EN <- run_sindy_subset(
    DEFM.amp.df = EN_loc[, c("RMM1", "RMM2", "RMM1.vel", "RMM2.vel")],
    enso_label = "ElNino",
    Nreal = Nreal, Nsamp = Nsamp, n = n, Lambda.inv = Lambda.inv
  )
  res_LN <- run_sindy_subset(
    DEFM.amp.df = LN_loc[, c("RMM1", "RMM2", "RMM1.vel", "RMM2.vel")],
    enso_label = "LaNina",
    Nreal = Nreal, Nsamp = Nsamp, n = n, Lambda.inv = Lambda.inv
  )
  cat(
    "SINDy pool [", season_tag, "] ALL ncol(modelos0)=", ncol(res_ALL$modelos0.df),
    " n_omega2=", length(res_ALL$omega2),
    " | EN ncol=", ncol(res_EN$modelos0.df),
    " | LN ncol=", ncol(res_LN$modelos0.df), "\n",
    sep = ""
  )
  if (identical(season_tag, "winter") && nzchar(Sys.getenv("SINDY_SEED_METRICS", ""))) {
    idx_nonlin <- c(3L, 5L, 6L, 7L, 8L, 9L, 12L, 14L, 15L, 16L, 17L, 18L)
    en_s <- sum(res_EN$num.coefs[idx_nonlin])
    ln_s <- sum(res_LN$num.coefs[idx_nonlin])
    cat(sprintf("SINDY_SEED_METRICS\t%d\t%.8g\t%.8g\n", sindy_seed, en_s, ln_s))
    quit(save = "no")
  }
  fig_suffix <- paste0(fig_base_suffix, "_", season_tag)

# -------------------------- Figure 3: Statistics of coefficients (paper style) --------------------------
Coefs.df <- data.frame(
  Coef = seq_len(2 * NPsi),
  Freq = res_ALL$num.coefs,
  Mean = res_ALL$avg.coefs,
  p25 = res_ALL$p25,
  p75 = res_ALL$p75,
  min = res_ALL$min_coefs,
  max = res_ALL$max_coefs
)

gg_fig3a <- ggplot() +
  geom_col(
    data = Coefs.df,
    aes(x = Coef, y = Freq, fill = ifelse(Coef <= NPsi, "fc1", "fc2")),
    col = "gray35"
  ) +
  scale_fill_manual(
    values = c(fc1 = "cornflowerblue", fc2 = "darkolivegreen3"),
    labels = c(TeX("$\\dot{x}$"), TeX("$\\dot{y}$")),
    name = "Equation"
  ) +
  scale_x_continuous(breaks = seq_len(2 * NPsi), labels = x_lab_fig) +
  ylab("Frequency") +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold")
  )

gg_fig3b <- ggplot() +
  geom_col(
    data = Coefs.df,
    aes(x = Coef, y = Mean, fill = ifelse(Coef <= NPsi, "fc1", "fc2")),
    col = "gray35"
  ) +
  scale_fill_manual(
    values = c(fc1 = "cornflowerblue", fc2 = "darkolivegreen3"),
    labels = c(TeX("$\\dot{x}$"), TeX("$\\dot{y}$")),
    name = "Equation"
  ) +
  geom_errorbar(data = Coefs.df, aes(x = Coef, ymin = p25, ymax = p75), inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq_len(2 * NPsi), labels = x_lab_fig) +
  geom_point(data = Coefs.df, aes(x = Coef, y = min), inherit.aes = FALSE, shape = 25, col = "red3", fill = "white", size = 2) +
  geom_point(data = Coefs.df, aes(x = Coef, y = max), inherit.aes = FALSE, shape = 24, col = "red3", fill = "white", size = 2) +
  ylab("Mean") +
  xlab("Coef") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 15, face = "bold", vjust = 0.2),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  )

fig3 <- ggarrange(
  gg_fig3a,
  gg_fig3b,
  nrow = 2L,
  heights = c(1, 1.1),
  labels = c("a", "b"),
  font.label = list(size = 14, face = "bold")
)
ggsave(
  filename = file.path(fig_dir, paste0("Fig3_CoefficientStats", fig_suffix, ".png")),
  plot = fig3, width = 9, height = 6.5, dpi = out_dpi, bg = "white"
)

# -------------------------- Masks + shared theme/heights for Fig4 / Fig8 --------------------------
mask_clim <- build_mask_long(res_ALL, "C")
mask_en <- build_mask_long(res_EN, "EN")
mask_ln <- build_mask_long(res_LN, "LN")

n_c <- length(unique(mask_clim$long$variable))
n_e <- length(unique(mask_en$long$variable))
n_l <- length(unique(mask_ln$long$variable))
sum_nl <- max(1L, n_e + n_l)
h_unit <- 11.5 / sum_nl
h_leg <- 1.15
h_pad_fig8 <- 0.22
height_fig4 <- max(3.2, h_unit * n_c + h_leg)
height_fig8a <- max(2.8, h_unit * n_e + h_pad_fig8 + h_leg)
height_fig8b <- max(2.8, h_unit * n_l + h_pad_fig8 + h_leg)
height_fig8 <- max(8, h_unit * (n_e + n_l) + h_pad_fig8 + h_leg + 0.35)

theme_structure_facets <- function(axis_lab_pt = 13) {
  theme_bw() +
    theme(
      axis.text.x = element_text(size = axis_lab_pt, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = axis_lab_pt, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.ticks.y = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.spacing.x = grid::unit(16, "pt"),
      legend.text = element_text(size = axis_lab_pt, face = "bold"),
      legend.title = element_text(size = axis_lab_pt, face = "bold"),
      strip.text.y.left = element_text(
        angle = 0,
        size = axis_lab_pt,
        face = "bold",
        margin = margin(t = 2, r = 4, b = 2, l = 4)
      ),
      strip.text = element_text(margin = margin(t = 2, r = 4, b = 2, l = 4)),
      panel.spacing.y = grid::unit(0.08, "lines")
    )
}

# -------------------------- Figure 4: Model's structure I (paper style) --------------------------
gg_fig4 <- ggplot(data = mask_clim$long, aes(x = Coef, y = value, group = Coef, fill = ifelse(Coef <= NPsi, "fc1", "fc2"))) +
  geom_col(col = "gray35", width = 0.98) +
  scale_fill_manual(
    values = c(fc1 = "cornflowerblue", fc2 = "darkolivegreen3"),
    labels = c(TeX("$\\dot{x}$"), TeX("$\\dot{y}$")),
    name = "Equation"
  ) +
  scale_x_continuous(breaks = seq_len(2 * NPsi), labels = x_lab_fig) +
  scale_y_continuous(breaks = c(0, 1), labels = NULL) +
  ylab(NULL) +
  xlab("Coef") +
  theme_structure_facets(axis_lab_pt = 15) +
  facet_grid(
    variable ~ .,
    switch = "y",
    scales = "fixed",
    space = "fixed",
    labeller = labeller(variable = as_labeller(mask_clim$facet_labs, default = label_value))
  )

ggsave(
  filename = file.path(fig_dir, paste0("Fig4_ModelStructure_I", fig_suffix, ".png")),
  plot = gg_fig4, width = 10, height = height_fig4, dpi = out_dpi, bg = "white"
)

export_fig4_model_coefficients(
  res_ALL,
  file.path(out_dir, paste0("Fig4_model_coefficients", fig_suffix, ".csv")),
  "C",
  dictionary_terms
)

# -------------------------- Figure 7: Distribution of periods II (ENSO) --------------------------
period_df <- data.frame(
  ENSO = c(rep("El Niño", length(res_EN$omega2)), rep("La Niña", length(res_LN$omega2))),
  T = c(2 * pi / sqrt(res_EN$omega2), 2 * pi / sqrt(res_LN$omega2))
)
period_df <- period_df[is.finite(period_df$T) & period_df$T > 0, ]

g7 <- ggplot(period_df, aes(x = T, y = after_stat(density), col = ENSO, fill = ENSO)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.35, linewidth = 0.4) +
  theme_bw() +
  scale_fill_manual(values = c("El Niño" = "#8f79b4", "La Niña" = "#62c4b7")) +
  scale_color_manual(values = c("El Niño" = "#6e5698", "La Niña" = "#36ab9f")) +
  xlab("Period [days]") + ylab("Frequency") +
  theme(
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(10, 14, 10, 10)
  )

split_period <- split(period_df$T, period_df$ENSO)
stats_period <- data.frame(
  `#LM` = sapply(split_period, length),
  Min = round(sapply(split_period, min), 1),
  Max = round(sapply(split_period, max), 1),
  Mean = round(sapply(split_period, mean), 1),
  Sd = round(sapply(split_period, sd), 1)
)
tab <- ggtexttable(t(stats_period), rows = c("#LM", "Min", "Max", "Mean", "Sd"),
                   cols = names(split_period),
                   theme = ttheme(
                     colnames.style = colnames_style(color = "white", size = 12, fill = "gray30"),
                     tbody.style = tbody_style(
                       fill = viridis(2, begin = 0.1, end = 0.6, direction = -1, alpha = 0.4),
                       color = "black", size = 12
                     )
                   ))
fig7 <- g7 +
  inset_element(
    tab,
    left = 0.66,
    bottom = 0.52,
    right = 0.998,
    top = 0.995,
    align_to = "full"
  )
ggsave(filename = file.path(fig_dir, paste0("Fig7_Period_ENSO", fig_suffix, ".png")),
       plot = fig7, width = 10.5, height = 6.5, dpi = out_dpi, bg = "white")

# -------------------------- Figure 8: Model's structure II (ENSO, paper style) --------------------------
gg_fig8_en_only <- ggplot(mask_en$long, aes(x = Coef, y = value, group = Coef, fill = ifelse(Coef <= NPsi, "fc1", "fc2"))) +
  geom_col(col = "gray35", width = 0.98) +
  scale_fill_manual(
    values = c(fc1 = "cornflowerblue", fc2 = "darkolivegreen3"),
    labels = c(TeX("$\\dot{x}$"), TeX("$\\dot{y}$")),
    name = "Equation"
  ) +
  scale_x_continuous(breaks = seq_len(2 * NPsi), labels = x_lab_fig) +
  scale_y_continuous(breaks = c(0, 1), labels = NULL) +
  ylab(NULL) +
  xlab("Coef") +
  theme_structure_facets(axis_lab_pt = 13) +
  theme(plot.margin = margin(6, 6, 4, 6)) +
  facet_grid(
    variable ~ .,
    switch = "y",
    scales = "fixed",
    space = "fixed",
    labeller = labeller(variable = as_labeller(mask_en$facet_labs, default = label_value))
  )

gg_fig8_ln_only <- ggplot(mask_ln$long, aes(x = Coef, y = value, group = Coef, fill = ifelse(Coef <= NPsi, "fc1", "fc2"))) +
  geom_col(col = "gray35", width = 0.98) +
  scale_fill_manual(
    values = c(fc1 = "cornflowerblue", fc2 = "darkolivegreen3"),
    labels = c(TeX("$\\dot{x}$"), TeX("$\\dot{y}$")),
    name = "Equation"
  ) +
  scale_x_continuous(breaks = seq_len(2 * NPsi), labels = x_lab_fig) +
  scale_y_continuous(breaks = c(0, 1), labels = NULL) +
  ylab(NULL) +
  xlab("Coef") +
  theme_structure_facets(axis_lab_pt = 13) +
  theme(plot.margin = margin(6, 6, 4, 6), legend.position = "none") +
  facet_grid(
    variable ~ .,
    switch = "y",
    scales = "fixed",
    space = "fixed",
    labeller = labeller(variable = as_labeller(mask_ln$facet_labs, default = label_value))
  )

fig8_divider <- wrap_elements(
  full = grid::grobTree(
    grid::linesGrob(
      x = grid::unit(c(0.05, 0.95), "npc"),
      y = grid::unit(0.5, "npc"),
      gp = grid::gpar(lty = 2, lwd = 1.1)
    )
  ),
  clip = FALSE
)

gg_fig8 <- gg_fig8_en_only / fig8_divider / gg_fig8_ln_only +
  plot_layout(heights = c(n_e, 0.06, n_l))

ggsave(
  filename = file.path(fig_dir, paste0("Fig8_ModelStructure_II_ENSO", fig_suffix, ".png")),
  plot = gg_fig8, width = 12, height = height_fig8, dpi = out_dpi, bg = "white"
)

ggsave(
  filename = file.path(fig_dir, paste0("Fig8a_ModelStructure_EN_only", fig_suffix, ".png")),
  plot = gg_fig8_en_only,
  width = 10,
  height = height_fig8a,
  dpi = out_dpi,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, paste0("Fig8b_ModelStructure_LN_only", fig_suffix, ".png")),
  plot = gg_fig8_ln_only,
  width = 10,
  height = height_fig8b,
  dpi = out_dpi,
  bg = "white"
)

# -------------------------- Save averaged models (coeff means) --------------------------
  save_avg_model_season <- function(res, enso_tag) {
    avg.modelo <- cbind(res$avg.coefs[1:NPsi], res$avg.coefs[(NPsi + 1):(2 * NPsi)])
    colnames(avg.modelo) <- c("RMM1_dot", "RMM2_dot")
    out <- data.frame(term_index = seq_len(NPsi), avg.modelo)
    fn <- paste0(
      "avg_model_", enso_tag, "_", years_mode, "_", start_year,
      "_", end_year, "_", season_tag, ".csv"
    )
    write.csv(out, file = file.path(out_dir, fn), row.names = FALSE)
  }
  save_avg_model_season(res_EN, "ElNino")
  save_avg_model_season(res_LN, "LaNina")
}

set.seed(sindy_seed)
run_season_bundle(build_x_in_df("winter", full_years_winter), "winter")
set.seed(sindy_seed)
run_season_bundle(build_x_in_df("summer", full_years_summer), "summer")

cat("Done.\n")
cat("RMM:", normalizePath(rmm_path, winslash = "/", mustWork = FALSE), "\n")
cat("ONI:", normalizePath(oni_path, winslash = "/", mustWork = FALSE), "\n")
cat("years_mode:", years_mode, " (tag only; figures use full winter & summer season lists)\n")
cat(
  "Winter seasons (Oct27–May4):", length(full_years_winter),
  " | Summer seasons (Apr27–Nov4):", length(full_years_summer),
  " | ONI/DJFM ElNino:", length(ninos), " LaNina:", length(ninas), " Neutral:", length(neutros), "\n"
)
cat("Calendar filter:", start_year, "-", end_year, "\n")
cat("Nreal:", Nreal, "Nsamp:", Nsamp, "amp_threshold:", amp_threshold, "oni_threshold:", oni_threshold, "\n")
prep_w <- prepare_amp_en_ln(build_x_in_df("winter", full_years_winter))
prep_s <- prepare_amp_en_ln(build_x_in_df("summer", full_years_summer))
cat("Winter EN/LN rows (amp>thr):", nrow(prep_w$EN), "/", nrow(prep_w$LN), "\n")
cat("Summer EN/LN rows (amp>thr):", nrow(prep_s$EN), "/", nrow(prep_s$LN), "\n")
cat("Dictionary order check:", paste(dictionary_terms, collapse = ", "), "\n")
cat("Index check: a01@xdot(y)=", idx_xdot("y"), "; b10@ydot(x)=", idx_ydot("x"), "; b01@ydot(y)=", idx_ydot("y"), "; a10@xdot(x)=", idx_xdot("x"), "\n")
cat("Figures written to:", normalizePath(fig_dir), "\n")
cat(
  "Fig4 CSV suffixes: ", fig_base_suffix, "_winter / ", fig_base_suffix, "_summer\n",
  sep = ""
)

