###############################################################################
#  oreo self-test  --  verifies the SPP implementation against exact results
#
#  Usage (after installing oreo):
#      source(system.file("tests", "oreo_selftest.R", package = "oreo"))
#  or simply
#      source("oreo_selftest.R")
#
#  Only base R + oreo are needed. The script prints a PASS/FAIL table and
#  ends with "ALL TESTS PASSED" or a list of failed checks.
#
#  What is checked
#   T1  Linear viscoelastic (LVE) benchmark, rpp_fft():  G't = G', G''t = G'' (Pa),
#       dG't/dt = dG''t/dt = speed = 0, sigma_d = 0, delta_t = atan(G''/G').
#   T2  Same benchmark with norm_rate = FALSE (oreo 1.0 convention):
#       G''t = G''/omega and sigma_d still = 0 (internal consistency).
#   T3  Nonlinear signal (3rd, 5th, 7th harmonics), rpp_fft() vs an independent
#       analytical Frenet-Serret reference -- this exercises the third-derivative
#       branch (dG't/dt, dG''t/dt, speed) that was wrong in oreo 1.0.
#   T4  LVE benchmark with the numerical-differentiation route, Rpp_num().
#   T5  Nonlinear signal: Rpp_num() agrees with rpp_fft().
#   T6  Smoke test on the example dataset shipped with the package.
#   T7  Phase-shifted strain, gamma = g0 sin(w t + phi), rpp_fft(): the internal
#       phase alignment (Delta) must recover the analytical reference.
#       phi = 2.2 rad exercises the "Delta + pi" branch; phi = pi/2 is the
#       common "strain starts as a cosine" case.
#   T8  Same phase-shifted signal with Rpp_num() (no internal shift).
#   T9  Three cycles (p = 3) plus phase shift, rpp_fft().
#   T10 Three cycles plus phase shift, Rpp_num() (looped and standard modes).
###############################################################################

suppressPackageStartupMessages(library(oreo))

cat("\n=== oreo self-test ===\n")
cat("oreo version :", as.character(utils::packageVersion("oreo")), "\n")
cat("R version    :", R.version.string, "\n\n")

has_norm <- "norm_rate" %in% names(formals(oreo::rpp_fft))
if (!has_norm) {
  cat("NOTE: this oreo build has no 'norm_rate' argument (oreo 1.0).\n",
      "     T1, T3, T4, T5 are EXPECTED TO FAIL on this version.\n\n", sep = "")
}

## ---------------------------------------------------------------- helpers ----
results <- data.frame(test = character(), check = character(),
                      value = character(), tol = character(),
                      status = character(), stringsAsFactors = FALSE)

check <- function(test, what, err, tol) {
  ok <- is.finite(err) && err <= tol
  results[nrow(results) + 1, ] <<- list(test, what, formatC(err, format = "e", digits = 2),
                                        formatC(tol, format = "e", digits = 0),
                                        if (ok) "PASS" else "FAIL")
  invisible(ok)
}
relerr <- function(x, ref) max(abs(x - ref)) / max(abs(ref))
abserr <- function(x, scale) max(abs(x)) / scale

call_fft <- function(tw, rw, L, w, M, norm) {
  if (has_norm) rpp_fft(tw, rw, L = L, omega = w, M = M, p = 1, norm_rate = norm)
  else          rpp_fft(tw, rw, L = L, omega = w, M = M, p = 1)
}
call_num <- function(tw, rw, L, k, mode, norm) {
  if (has_norm) Rpp_num(tw, rw, L = L, k = k, num_mode = mode, norm_rate = norm)
  else          Rpp_num(tw, rw, L = L, k = k, num_mode = mode)
}

## Independent analytical SPP reference, normalized coordinate [g, gdot/w, s]
## Input: stress = sum_n (a_n sin(n w t) + b_n cos(n w t)), strain = g0 sin(w t)
spp_reference <- function(t, w, g0, n, a, b) {
  S  <- function(d) {                     # d-th time derivative of stress
    out <- 0
    for (i in seq_along(n)) {
      ph <- n[i] * w * t
      out <- out + (n[i] * w)^d * (a[i] * sin(ph + d * pi / 2) + b[i] * cos(ph + d * pi / 2))
    }
    out
  }
  G  <- function(d) g0 * w^d * sin(w * t + d * pi / 2)          # strain
  R  <- function(d) g0 * w^(d + 1) * cos(w * t + d * pi / 2) / w # rate / w
  rd   <- cbind(G(1), R(1), S(1))
  rdd  <- cbind(G(2), R(2), S(2))
  rddd <- cbind(G(3), R(3), S(3))
  cr   <- cbind(rd[,2]*rdd[,3] - rd[,3]*rdd[,2],
                rd[,3]*rdd[,1] - rd[,1]*rdd[,3],
                rd[,1]*rdd[,2] - rd[,2]*rdd[,1])
  Gp  <- -cr[,1] / cr[,3]
  Gpp <- -cr[,2] / cr[,3]
  tr  <- rowSums(rddd * cr)
  Gpd  <- -rd[,2] * tr / cr[,3]^2
  Gppd <-  rd[,1] * tr / cr[,3]^2
  list(Gp = Gp, Gpp = Gpp, Gpd = Gpd, Gppd = Gppd, speed = sqrt(Gpd^2 + Gppd^2),
       disp = S(0) - (Gp * G(0) + Gpp * R(0)))
}

## ------------------------------------------------------------ test data ----
w  <- 3.16; g0 <- 0.1; L <- 1024
tt <- (0:(L - 1)) * 2 * pi / w / L                  # exactly one cycle
strain <- g0 * sin(w * tt)
rate   <- g0 * w * cos(w * tt)

## LVE: G' = 1000 Pa, G'' = 250 Pa
Gp0 <- 1000; Gpp0 <- 250
stress_lve <- Gp0 * strain + Gpp0 * rate / w
rw_lve <- data.frame(strain, rate, stress_lve)

## Nonlinear: odd harmonics up to 7
nh <- c(1, 3, 5, 7)
a  <- g0 * c(1000, -120,  30, -8)                  # sin coefficients
b  <- g0 * c( 250,   60, -15,  4)                  # cos coefficients
stress_nl <- rowSums(sapply(seq_along(nh), function(i)
  a[i] * sin(nh[i] * w * tt) + b[i] * cos(nh[i] * w * tt)))
rw_nl <- data.frame(strain, rate, stress_nl)
ref   <- spp_reference(tt, w, g0, nh, a, b)

## ------------------------------------------------------------------ T1 ----
o  <- call_fft(tt, rw_lve, L, w, 15, TRUE)$spp_data_out
sc <- Gp0 * w                                     # natural scale for dG/dt
check("T1 fft LVE", "Gp_t = 1000",           relerr(o$Gp_t,  rep(Gp0,  L)), 1e-8)
check("T1 fft LVE", "Gpp_t = 250 (Pa)",      relerr(o$Gpp_t, rep(Gpp0, L)), 1e-8)
check("T1 fft LVE", "delta_t = atan(0.25)",  relerr(o$delta_t, rep(atan(Gpp0/Gp0), L)), 1e-8)
check("T1 fft LVE", "Gp_t_dot = 0",          abserr(o$Gp_t_dot,  sc), 1e-8)
check("T1 fft LVE", "Gpp_t_dot = 0",         abserr(o$Gpp_t_dot, sc), 1e-8)
check("T1 fft LVE", "G_speed = 0",           abserr(o$G_speed,   sc), 1e-8)
check("T1 fft LVE", "disp_stress = 0",       abserr(o$disp_stress, max(abs(stress_lve))), 1e-8)
check("T1 fft LVE", "delta_t_dot = 0",       abserr(o$delta_t_dot, 1), 1e-8)

## ------------------------------------------------------------------ T2 ----
o2 <- call_fft(tt, rw_lve, L, w, 15, FALSE)$spp_data_out
check("T2 fft raw-rate", "Gpp_t = 250/omega",      relerr(o2$Gpp_t, rep(Gpp0 / w, L)), 1e-8)
check("T2 fft raw-rate", "disp_stress = 0",        abserr(o2$disp_stress, max(abs(stress_lve))), 1e-8)
check("T2 fft raw-rate", "Gpp_t(raw)*omega = Gpp_t(norm)", relerr(o2$Gpp_t * w, o$Gpp_t), 1e-8)

## ------------------------------------------------------------------ T3 ----
o3 <- call_fft(tt, rw_nl, L, w, 15, TRUE)$spp_data_out
check("T3 fft nonlinear", "Gp_t vs analytic",      relerr(o3$Gp_t,  ref$Gp),   1e-6)
check("T3 fft nonlinear", "Gpp_t vs analytic",     relerr(o3$Gpp_t, ref$Gpp),  1e-6)
check("T3 fft nonlinear", "Gp_t_dot vs analytic",  relerr(o3$Gp_t_dot,  ref$Gpd),  1e-6)
check("T3 fft nonlinear", "Gpp_t_dot vs analytic", relerr(o3$Gpp_t_dot, ref$Gppd), 1e-6)
check("T3 fft nonlinear", "G_speed vs analytic",   relerr(o3$G_speed,   ref$speed), 1e-6)
check("T3 fft nonlinear", "disp_stress vs analytic", relerr(o3$disp_stress, ref$disp), 1e-6)

## ------------------------------------------------------------------ T4 ----
o4 <- call_num(tt, rw_lve, L, 1, 2, TRUE)$spp_data_out
check("T4 num LVE", "Gp_t = 1000",       relerr(o4$Gp_t,  rep(Gp0,  L)), 1e-4)
check("T4 num LVE", "Gpp_t = 250 (Pa)",  relerr(o4$Gpp_t, rep(Gpp0, L)), 1e-4)
check("T4 num LVE", "Gp_t_dot ~ 0",      abserr(o4$Gp_t_dot,  sc), 1e-3)
check("T4 num LVE", "disp_stress ~ 0",   abserr(o4$disp_stress, max(abs(stress_lve))), 1e-4)

## ------------------------------------------------------------------ T5 ----
o5 <- call_num(tt, rw_nl, L, 1, 2, TRUE)$spp_data_out
check("T5 num vs fft", "Gp_t",      relerr(o5$Gp_t,      o3$Gp_t),      1e-3)
check("T5 num vs fft", "Gpp_t",     relerr(o5$Gpp_t,     o3$Gpp_t),     1e-3)
check("T5 num vs fft", "Gp_t_dot",  relerr(o5$Gp_t_dot,  o3$Gp_t_dot),  1e-2)
check("T5 num vs fft", "Gpp_t_dot", relerr(o5$Gpp_t_dot, o3$Gpp_t_dot), 1e-2)

## Helper: sample the nonlinear signal at shifted time tau = t + phi/w
## (strain = g0 sin(w t + phi)); returns the input frame and the reference.
make_signal <- function(t, phi) {
  tau <- t + phi / w
  st  <- rowSums(sapply(seq_along(nh), function(i)
           a[i] * sin(nh[i] * w * tau) + b[i] * cos(nh[i] * w * tau)))
  list(rw = data.frame(g0 * sin(w * tau), g0 * w * cos(w * tau), st),
       ref_in = spp_reference(tau, w, g0, nh, a, b))
}
## rpp_fft() returns one cycle re-aligned so that strain = g0 sin(w tau),
## on tau = (0:(L-1)) * 2*pi/w/L, whatever the number of cycles in the input.
ref_fft <- ref

## ------------------------------------------------------------------ T7 ----
phi <- 2.2
s7  <- make_signal(tt, phi)
o7  <- call_fft(tt, s7$rw, L, w, 15, TRUE)$spp_data_out
check("T7 fft phase 2.2", "strain realigned to sine", relerr(o7$strain, g0 * sin(w * tt)), 1e-8)
check("T7 fft phase 2.2", "Gp_t vs analytic",      relerr(o7$Gp_t,      ref_fft$Gp),    1e-6)
check("T7 fft phase 2.2", "Gpp_t vs analytic",     relerr(o7$Gpp_t,     ref_fft$Gpp),   1e-6)
check("T7 fft phase 2.2", "Gp_t_dot vs analytic",  relerr(o7$Gp_t_dot,  ref_fft$Gpd),   1e-6)
check("T7 fft phase 2.2", "Gpp_t_dot vs analytic", relerr(o7$Gpp_t_dot, ref_fft$Gppd),  1e-6)
check("T7 fft phase 2.2", "G_speed vs analytic",   relerr(o7$G_speed,   ref_fft$speed), 1e-6)
check("T7 fft phase 2.2", "disp_stress vs analytic", relerr(o7$disp_stress, ref_fft$disp), 1e-6)

s7c <- make_signal(tt, pi / 2)
o7c <- call_fft(tt, s7c$rw, L, w, 15, TRUE)$spp_data_out
check("T7 fft cosine", "strain realigned to sine", relerr(o7c$strain, g0 * sin(w * tt)), 1e-8)
check("T7 fft cosine", "Gp_t vs analytic",       relerr(o7c$Gp_t,     ref_fft$Gp),  1e-6)
check("T7 fft cosine", "Gp_t_dot vs analytic",   relerr(o7c$Gp_t_dot, ref_fft$Gpd), 1e-6)

## ------------------------------------------------------------------ T8 ----
o8 <- call_num(tt, s7$rw, L, 1, 2, TRUE)$spp_data_out
r8 <- s7$ref_in                                  # reference at the input times
check("T8 num phase 2.2", "Gp_t vs analytic",      relerr(o8$Gp_t,      r8$Gp),   1e-4)
check("T8 num phase 2.2", "Gpp_t vs analytic",     relerr(o8$Gpp_t,     r8$Gpp),  1e-4)
check("T8 num phase 2.2", "Gp_t_dot vs analytic",  relerr(o8$Gp_t_dot,  r8$Gpd),  1e-2)
check("T8 num phase 2.2", "Gpp_t_dot vs analytic", relerr(o8$Gpp_t_dot, r8$Gppd), 1e-2)
check("T8 num phase 2.2", "disp_stress vs analytic", relerr(o8$disp_stress, r8$disp), 1e-3)

## ------------------------------------------------------------------ T9 ----
pc  <- 3; Lp <- pc * L
ttp <- (0:(Lp - 1)) * 2 * pi / w / L              # three full cycles
s9  <- make_signal(ttp, phi)
o9  <- (if (has_norm) rpp_fft(ttp, s9$rw, L = Lp, omega = w, M = 15, p = pc, norm_rate = TRUE)
        else rpp_fft(ttp, s9$rw, L = Lp, omega = w, M = 15, p = pc))$spp_data_out
## output: Lp points over ONE cycle -> reference on the matching fine grid
tau9 <- (0:(Lp - 1)) * 2 * pi / w / Lp
r9   <- spp_reference(tau9, w, g0, nh, a, b)
check("T9 fft 3 cycles", "strain realigned to sine", relerr(o9$strain, g0 * sin(w * tau9)), 1e-8)
check("T9 fft 3 cycles", "Gp_t vs analytic",      relerr(o9$Gp_t,      r9$Gp),    1e-6)
check("T9 fft 3 cycles", "Gpp_t vs analytic",     relerr(o9$Gpp_t,     r9$Gpp),   1e-6)
check("T9 fft 3 cycles", "Gp_t_dot vs analytic",  relerr(o9$Gp_t_dot,  r9$Gpd),   1e-6)
check("T9 fft 3 cycles", "Gpp_t_dot vs analytic", relerr(o9$Gpp_t_dot, r9$Gppd),  1e-6)
check("T9 fft 3 cycles", "disp_stress vs analytic", relerr(o9$disp_stress, r9$disp), 1e-6)

## ----------------------------------------------------------------- T10 ----
r10 <- s9$ref_in
o10 <- call_num(ttp, s9$rw, Lp, 1, 2, TRUE)$spp_data_out
check("T10 num 3 cyc loop", "Gp_t vs analytic",      relerr(o10$Gp_t,      r10$Gp),   1e-4)
check("T10 num 3 cyc loop", "Gpp_t vs analytic",     relerr(o10$Gpp_t,     r10$Gpp),  1e-4)
check("T10 num 3 cyc loop", "Gp_t_dot vs analytic",  relerr(o10$Gp_t_dot,  r10$Gpd),  1e-2)
check("T10 num 3 cyc loop", "Gpp_t_dot vs analytic", relerr(o10$Gpp_t_dot, r10$Gppd), 1e-2)
o10s <- call_num(ttp, s9$rw, Lp, 1, 1, TRUE)$spp_data_out
inn  <- 10:(Lp - 10)                              # standard mode: skip one-sided ends
check("T10 num 3 cyc std", "Gp_t vs analytic (interior)",     relerr(o10s$Gp_t[inn],     r10$Gp[inn]),   1e-4)
check("T10 num 3 cyc std", "Gpp_t vs analytic (interior)",    relerr(o10s$Gpp_t[inn],    r10$Gpp[inn]),  1e-4)
check("T10 num 3 cyc std", "Gp_t_dot vs analytic (interior)", relerr(o10s$Gp_t_dot[inn], r10$Gpd[inn]),  1e-2)

## ------------------------------------------------------------------ T6 ----
ok6 <- tryCatch({
  e <- new.env(); utils::data("mydata", package = "oreo", envir = e)
  df <- rpp_read2(e$mydata, selected = c(2, 3, 4, 0, 0, 1, 0, 0))
  rw <- data.frame(df$strain, df$strain_rate, df$stress)
  f6 <- rpp_fft(df$raw_time, rw, L = 1024, omega = 3.16, M = 15, p = 1)$spp_data_out
  n6 <- Rpp_num(df$raw_time, rw, L = 1024, k = 8, num_mode = 1)$spp_data_out
  all(is.finite(as.matrix(f6))) && all(is.finite(as.matrix(n6)))
}, error = function(e) { cat("T6 error:", conditionMessage(e), "\n"); FALSE })
check("T6 example data", "rpp_fft + Rpp_num run, all finite", if (isTRUE(ok6)) 0 else Inf, 0)

## --------------------------------------------------------------- report ----
cat(sprintf("%-18s %-34s %-10s %-7s %s\n", "test", "check", "error", "tol", "status"))
cat(strrep("-", 78), "\n")
for (i in seq_len(nrow(results)))
  cat(sprintf("%-18s %-34s %-10s %-7s %s\n", results$test[i], results$check[i],
              results$value[i], results$tol[i], results$status[i]))
cat(strrep("-", 78), "\n")
nf <- sum(results$status == "FAIL")
if (nf == 0) {
  cat(sprintf("ALL TESTS PASSED (%d checks)\n\n", nrow(results)))
} else {
  cat(sprintf("%d of %d CHECKS FAILED:\n", nf, nrow(results)))
  print(results[results$status == "FAIL", c("test", "check")], row.names = FALSE)
  cat("\n")
}
invisible(nf == 0)
