| Block | Label | Criterion | Objective | Analysis prior | Value | N total | A(N) | A(N - 1) | Ledger ID |
|---|---|---|---|---|---|---|---|---|---|
| A_inputs | Pilot counts (TP, FN, FP, TN) |  |  |  | 43, 11, 0, 358 (N = 412; prevalence 54/412 = 13.1%) |  |  |  | T5-01 |
| A_inputs | Base priors and power-prior weight |  |  |  | Beta(1, 1) for Se, Sp and prevalence; a0 = 0.5 |  |  |  | T5-04 |
| A_inputs | Design prior, Se |  |  |  | Beta(22.5, 6.5) |  |  |  | T5-05 |
| A_inputs | Design prior, Sp |  |  |  | Beta(180, 1) |  |  |  | T5-05 |
| A_inputs | Design prior, prevalence |  |  |  | Beta(28, 180) |  |  |  | T5-05 |
| B_main | design_onesided_joint | one-sided 95%, lower half-widths 0.1 and 0.05 | joint (Se and Sp) | flat Beta(1, 1) (use_in = "design") |  | 497 | 0.8007 | 0.7985 | T5-06 |
| B_main | both_onesided_joint | one-sided 95%, lower half-widths 0.1 and 0.05 | joint (Se and Sp) | updated pilot prior (use_in = "both") |  | 258 | 0.8012 | 0.7978 | T5-07 |
| B_main | design_twosided_joint | two-sided 95%, full widths 0.2 and 0.1 | joint (Se and Sp) | flat Beta(1, 1) (use_in = "design") |  | 612 | 0.8002 | 0.7987 | T5-08 |
| B_main | both_twosided_joint | two-sided 95%, full widths 0.2 and 0.1 | joint (Se and Sp) | updated pilot prior (use_in = "both") |  | 363 | 0.8017 | 0.7996 | T5-08 |
| B_main | design_onesided_se | one-sided 95%, lower half-width 0.1 (Se) | marginal (Se only) | flat Beta(1, 1) (use_in = "design") |  | 497 | 0.8007 | 0.7985 | T5-09 |
| B_main | design_twosided_se | two-sided 95%, full width 0.2 (Se) | marginal (Se only) | flat Beta(1, 1) (use_in = "design") |  | 612 | 0.8002 | 0.7987 | T5-09 |
| C_sensitivity | a0_0.25_design | one-sided 95%, lower half-widths 0.1 and 0.05; target 0.8 | joint (Se and Sp) | flat Beta(1, 1) (use_in = "design") |  | 520 | 0.8014 | 0.7998 | T5-10 |
| C_sensitivity | a0_1_design | one-sided 95%, lower half-widths 0.1 and 0.05; target 0.8 | joint (Se and Sp) | flat Beta(1, 1) (use_in = "design") |  | 480 | 0.8004 | 0.7977 | T5-10 |
| C_sensitivity | dse_0.075 | one-sided 95%, lower half-widths 0.075 and 0.05; target 0.8 | joint (Se and Sp) | flat Beta(1, 1) (use_in = "design") |  | 859 | 0.8009 | 0.7996 | T5-12 |
| C_sensitivity | dse_0.125 | one-sided 95%, lower half-widths 0.125 and 0.05; target 0.8 | joint (Se and Sp) | flat Beta(1, 1) (use_in = "design") |  | 325 | 0.8017 | 0.7985 | T5-12 |
| C_sensitivity | target_0.7 | one-sided 95%, lower half-widths 0.1 and 0.05; target 0.7 | joint (Se and Sp) | flat Beta(1, 1) (use_in = "design") |  | 457 | 0.7024 | 0.6997 | T5-14 |
| C_sensitivity | target_0.9 | one-sided 95%, lower half-widths 0.1 and 0.05; target 0.9 | joint (Se and Sp) | flat Beta(1, 1) (use_in = "design") |  | 558 | 0.9006 | 0.8994 | T5-14 |
| D_comparison | Buderer's formula (expected plug-in width), fixed values, one-sided 95% (Se 0.796, prevalence 0.131, Sp = Wilson lower limit 0.9894) | one-sided 95%, half-widths 0.10 and 0.05 | Se and Sp arms |  |  | 336 |  |  | T5-17 |
| D_comparison | Buderer's formula (expected plug-in width), fixed values, two-sided 95% (same inputs) | two-sided 95%, half-widths 0.10 and 0.05 | Se and Sp arms |  |  | 481 |  |  | T5-17 |
| D_comparison | Planned by Bulilete et al. (frequentist) | prevalence 15%, margin 5%, Se 90% |  |  |  | 927 |  |  | T5-03 |
| D_comparison | Enrolled by Bulilete et al. |  |  |  |  | 1369 |  |  | T5-02 |
| D_comparison | Analysed by Bulilete et al. |  |  |  |  | 1362 |  |  | T5-02 |
