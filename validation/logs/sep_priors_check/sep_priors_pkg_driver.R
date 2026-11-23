
suppressMessages(library(dtasamplesize))
emit <- function(id, key, val) cat(sprintf("%s|%s|%s\n", id, key, val))
num <- function(x) if (is.na(x)) "NA" else sprintf("%.17g", x)
emit("META", "version", as.character(packageVersion("dtasamplesize")))
emit("META", "libpath", find.package("dtasamplesize"))
run_case <- function(id, args, pilot, grid) {
  if (!is.null(pilot)) {
    pp <- do.call(bam_prior_from_pilot, pilot)
    for (nm in c("prior_se","prior_sp","prior_prev","analysis_prior_se","analysis_prior_sp")) {
      v <- pp[[nm]]
      emit(id, paste0("pilot_", nm), paste(sprintf("%.17g", v), collapse = ","))
    }
    args$priors <- pp
  }
  base <- c(list(n_range = 20:60, B = 1000, N_max = 1000L), args)
  r <- suppressWarnings(do.call(bam_sample_size, base))
  emit(id, "target_reached", as.integer(r$target_reached))
  emit(id, "N_total", num(r$N_total))
  emit(id, "joint_assurance", num(r$joint_assurance))
  emit(id, "max_assurance_evaluated", num(r$max_assurance_evaluated))
  emit(id, "N_at_max_assurance", num(r$N_at_max_assurance))
  emit(id, "search_type", r$search_type)
  emit(id, "eff_analysis_prior_se", paste(sprintf("%.17g", r$analysis_prior_se), collapse = ","))
  emit(id, "eff_analysis_prior_sp", paste(sprintf("%.17g", r$analysis_prior_sp), collapse = ","))
  Ns <- grid
  if (isTRUE(r$target_reached)) Ns <- c(Ns, r$N_total + (-1:1))
  Ns <- sort(unique(Ns[Ns >= 2]))
  for (N in Ns) {
    rr <- suppressWarnings(do.call(bam_sample_size, c(base, list(N_range = N))))
    emit(id, paste0("A@", N), num(rr$max_assurance_evaluated))
  }
}

run_case("C01", list(prior_se=c(17.0, 3.0), prior_sp=c(2.0, 2.0), prior_prev=c(4.0, 16.0), delta_se=0.14, delta_sp=0.1, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C02", list(prior_se=c(17.0, 3.0), prior_sp=c(2.0, 2.0), prior_prev=c(4.0, 16.0), analysis_prior_se=c(0.5, 0.5), analysis_prior_sp=c(0.5, 0.5), delta_se=0.14, delta_sp=0.1, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C03", list(prior_se=c(17.0, 3.0), prior_sp=c(2.0, 2.0), prior_prev=c(4.0, 16.0), analysis_prior_se=c(1.0, 1.0), analysis_prior_sp=c(1.0, 1.0), delta_se=0.14, delta_sp=0.1, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C04", list(prior_se=c(17.0, 3.0), prior_sp=c(2.0, 2.0), prior_prev=c(4.0, 16.0), analysis_prior_se=c(30.0, 5.0), analysis_prior_sp=c(10.0, 2.0), delta_se=0.14, delta_sp=0.1, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C05", list(prior_se=c(17.0, 3.0), prior_sp=c(2.0, 2.0), prior_prev=c(4.0, 16.0), analysis_prior_se=c(0.5, 0.5), delta_se=0.14, delta_sp=0.1, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C06", list(prior_se=c(1.0, 1.0), prior_sp=c(1.0, 1.0), prior_prev=c(4.0, 16.0), analysis_prior_se=c(17.0, 3.0), analysis_prior_sp=c(2.0, 2.0), delta_se=0.14, delta_sp=0.1, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C07", list(prior_se=c(25.0, 7.0), prior_sp=c(30.0, 2.0), prior_prev=c(13.56, 122.06), delta_se=0.1, delta_sp=0.05, target_assurance=0.8, alpha_ci=0.95, width_type="lower_half_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C08", list(prior_se=c(25.0, 7.0), prior_sp=c(30.0, 2.0), prior_prev=c(13.56, 122.06), analysis_prior_se=c(1.0, 1.0), analysis_prior_sp=c(1.0, 1.0), delta_se=0.1, delta_sp=0.05, target_assurance=0.8, alpha_ci=0.95, width_type="lower_half_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C09", list(prior_se=c(25.9, 2.1), prior_sp=c(2.0, 2.0), prior_prev=c(29.0, 98.0), delta_se=0.16, delta_sp=0.1, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="se"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C10", list(prior_se=c(17.0, 3.0), prior_sp=c(2.0, 2.0), prior_prev=c(4.0, 16.0), analysis_prior_se=c(0.5, 0.5), delta_se=0.14, delta_sp=0.1, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="se"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C11", list(prior_se=c(17.0, 3.0), prior_sp=c(2.0, 2.0), prior_prev=c(4.0, 16.0), delta_se=0.14, delta_sp=0.1, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="sp"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C12", list(prior_se=c(25.0, 7.0), prior_sp=c(30.0, 2.0), prior_prev=c(13.56, 122.06), analysis_prior_sp=c(1.0, 1.0), delta_se=0.1, delta_sp=0.05, target_assurance=0.8, alpha_ci=0.95, width_type="lower_half_width", assurance_for="sp"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C13", list(delta_se=0.16, delta_sp=0.2, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), list(tp=16L, fn=1L, tn=20L, fp=35L, base_prior_se=c(9.9, 1.1), base_prior_sp=c(1.0, 1.0), base_prior_prev=c(12.0, 43.0), a0=1.0, design="cohort", use_in="both"), c(2L, 5L, 20L, 100L, 250L))
run_case("C14", list(delta_se=0.16, delta_sp=0.25, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), list(tp=16L, fn=1L, tn=20L, fp=35L, base_prior_se=c(9.9, 1.1), base_prior_sp=c(1.0, 1.0), base_prior_prev=c(12.0, 43.0), a0=0.5, design="cohort", use_in="design"), c(2L, 5L, 20L, 100L, 250L))
run_case("C15", list(delta_se=0.2, delta_sp=0.2, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), list(tp=30L, fn=5L, tn=40L, fp=25L, base_prior_se=c(1.0, 1.0), base_prior_sp=c(1.0, 1.0), base_prior_prev=c(12.0, 43.0), a0=1.0, design="case_control", use_in="both"), c(2L, 5L, 20L, 100L, 250L))
run_case("C16", list(prior_se=c(17.0, 3.0), prior_sp=c(2.0, 2.0), prior_prev=c(4.0, 16.0), analysis_prior_se=c(0.5, 0.5), analysis_prior_sp=c(0.5, 0.5), delta_se=0.14, delta_sp=0.1, target_assurance=0.9, alpha_ci=0.9, width_type="full_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
run_case("C17", list(prior_se=c(17.0, 3.0), prior_sp=c(2.0, 2.0), prior_prev=c(4.0, 16.0), delta_se=0.02, delta_sp=0.02, target_assurance=0.8, alpha_ci=0.95, width_type="full_width", assurance_for="joint"), NULL, c(2L, 5L, 20L, 100L, 250L))
emit("META", "done", "1")
