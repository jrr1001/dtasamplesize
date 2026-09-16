# CONFIRMATORY property tests, PACKAGE side (PT1..PT6 of the locked grid).
suppressPackageStartupMessages({library(dtasamplesize); library(jsonlite)})
options(dtasamplesize.warn_small_B = FALSE)
grid <- fromJSON("LOCKED_confirmatory_grid.json", simplifyDataFrame = FALSE)
SC <- setNames(grid$scenarios, vapply(grid$scenarios, `[[`, "", "id"))
PT <- grid$property_tests

run <- function(s, target = s$target, dmult = 1, method = "exact",
                seed = 2026, B = 1000, N_range = NULL) {
  warns <- character(0)
  r <- withCallingHandlers(
    bam_sample_size(prior_se = unlist(s$prior_se),
                    prior_sp = unlist(s$prior_sp),
                    prior_prev = unlist(s$prior_prev),
                    delta_se = s$delta_se * dmult,
                    delta_sp = s$delta_sp * dmult,
                    target_assurance = target, alpha_ci = s$level,
                    n_range = 20:200, B = B, seed = seed,
                    N_range = if (is.null(N_range)) s$N_lo:s$N_hi else N_range,
                    method = method),
    warning = function(w) { warns <<- c(warns, conditionMessage(w))
                            invokeRestart("muffleWarning") })
  # "not found in range" -> +Inf, per the locked not_found_rule
  notfound <- any(grepl("achieved the target", warns))
  list(N = if (notfound) Inf else r$N_total, A = r$joint_assurance,
       raw_N = r$N_total, warns = warns)
}

pass <- list()
say <- function(id, ok, msg) {
  pass[[id]] <<- ok
  cat(sprintf("  [%s] %-4s %s\n", if (ok) "PASS" else "FAIL", id, msg))
}

cat("== PT1  larger target_assurance never reduces N (package) ==\n")
for (id in unlist(PT$PT1_target_monotone$package_scenarios)) {
  s <- SC[[id]]
  Ns <- vapply(unlist(PT$PT1_target_monotone$targets),
               function(t) run(s, target = t)$N, 0)
  say(paste0("PT1-", id), all(diff(Ns) >= 0),
      sprintf("targets %s -> N %s", paste(unlist(PT$PT1_target_monotone$targets),
              collapse = "/"), paste(Ns, collapse = "/")))
}

cat("== PT2  stricter width margin never reduces N (package) ==\n")
for (id in unlist(PT$PT2_delta_monotone$package_scenarios)) {
  s <- SC[[id]]
  Ns <- vapply(unlist(PT$PT2_delta_monotone$delta_multipliers),
               function(m) run(s, dmult = m)$N, 0)
  say(paste0("PT2-", id), all(diff(Ns) >= 0),
      sprintf("delta x%s -> N %s",
              paste(unlist(PT$PT2_delta_monotone$delta_multipliers),
                    collapse = "/x"), paste(Ns, collapse = "/")))
}

cat("== PT3  infeasible scenario is not passed off as valid ==\n")
s <- SC[["C12"]]
r <- run(s)
say("PT3-C12", length(r$warns) > 0 && r$A < s$target,
    sprintf("raw N_total=%s A=%.6f (<target %.2f: %s), warned: %s",
            r$raw_N, r$A, s$target, r$A < s$target, length(r$warns) > 0))
cat("      warning text: ", paste(r$warns, collapse = " | "), "\n")

cat("== PT4  monte_carlo reproducible with the same seed ==\n")
for (id in unlist(PT$PT4_reproducible_seed$scenarios)) {
  s <- SC[[id]]
  sd <- unlist(PT$PT4_reproducible_seed$seeds)
  a1 <- run(s, method = "monte_carlo", seed = sd[1], B = 40000)
  a2 <- run(s, method = "monte_carlo", seed = sd[1], B = 40000)
  b1 <- run(s, method = "monte_carlo", seed = sd[2], B = 40000)
  same <- identical(a1$raw_N, a2$raw_N) && identical(a1$A, a2$A)
  say(paste0("PT4-", id), same,
      sprintf("seed %d twice -> N %s/%s A %.6f/%.6f ; seed %d -> N %s A %.6f",
              sd[1], a1$raw_N, a2$raw_N, a1$A, a2$A, sd[2], b1$raw_N, b1$A))
}

cat("== PT5  RNG state of the calling session is restored ==\n")
for (m in unlist(PT$PT5_rng_state_restored$methods)) {
  s <- SC[[ unlist(PT$PT5_rng_state_restored$scenarios)[1] ]]
  set.seed(4242)
  before <- .Random.seed
  expected <- { set.seed(4242); runif(3) }   # what the session would draw
  set.seed(4242)
  invisible(run(s, method = m, seed = 999, B = 20000))
  after <- .Random.seed
  got <- runif(3)
  say(paste0("PT5-", m), identical(before, after) && isTRUE(all.equal(expected, got)),
      sprintf(".Random.seed identical: %s ; next runif(3) unaffected: %s",
              identical(before, after), isTRUE(all.equal(expected, got))))
}

cat("== PT6  under method='exact', seed and B do not move the headline ==\n")
for (id in unlist(PT$PT6_seed_irrelevant_under_exact$scenarios)) {
  s <- SC[[id]]
  combos <- expand.grid(seed = unlist(PT$PT6_seed_irrelevant_under_exact$seeds),
                        B = unlist(PT$PT6_seed_irrelevant_under_exact$B_values))
  res <- Map(function(sd, B) run(s, seed = sd, B = B), combos$seed, combos$B)
  Ns <- vapply(res, `[[`, 0, "raw_N"); As <- vapply(res, `[[`, 0, "A")
  say(paste0("PT6-", id),
      length(unique(Ns)) == 1 && length(unique(As)) == 1,
      sprintf("N unique=%d  A unique=%d  (N=%s A=%.10f)",
              length(unique(Ns)), length(unique(As)), Ns[1], As[1]))
}

cat("\n---- SUMMARY ----\n")
cat(sprintf("%d/%d property checks passed\n", sum(unlist(pass)), length(pass)))
if (any(!unlist(pass)))
  cat("FAILED:", paste(names(pass)[!unlist(pass)], collapse = ", "), "\n")
write_json(pass, "conf_properties.json", auto_unbox = TRUE, pretty = TRUE)
