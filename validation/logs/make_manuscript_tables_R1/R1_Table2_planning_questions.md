| Row | Planning question | Design assumption | Analysis interval | Computation | Function or script | N | Probability at N | MCSE |
|---|---|---|---|---|---|---|---|---|
| 1 | Expected plug-in width, each arm separately | Fixed planning values | Wald | Closed form | buderer_n() | 500 |  |  |
| 2 | Probability that both widths are attained | Fixed planning values | Wilson | Exact enumeration | validation/bridge_frequentist_bayesian.R | 578 | 0.8012 |  |
| 3 | Probability that both widths are attained | Fixed planning values | Wilson | Monte Carlo, step-10 grid (B = 20000, seed 2026) | joint_sample_size() | 580 | 0.8041 | 0.0028 |
| 4 | Probability that both widths are attained | Fixed planning values | Wilson | Monte Carlo, integer grid 100:800 (B = 20000, seed 2026; validation/comparator_calls.R) | joint_sample_size() | 579 | 0.8003 | 0.0028 |
| 5 | Credible interval with a vague prior | Fixed planning values | Jeffreys | Exact enumeration | validation/bridge_frequentist_bayesian.R | 576 | 0.8014 |  |
| 6 | Credible interval with a vague prior | Fixed planning values | Uniform | Exact enumeration | validation/bridge_frequentist_bayesian.R | 579 | 0.8025 |  |
| 7 | Informative analysis prior | Fixed planning values | Beta(17, 3), Beta(18, 2) | Exact enumeration | validation/bridge_frequentist_bayesian.R | 466 | 0.8002 |  |
| 8 | Prior uncertainty in the design | Beta design priors of the example | Wilson | Exact enumeration | validation/bridge_frequentist_bayesian.R | 831 | 0.8004 |  |
| 9 | Prior uncertainty in the design | Beta design priors of the example | Jeffreys | Exact enumeration | validation/bridge_frequentist_bayesian.R | 824 | 0.8000 |  |
| 10 | Prior uncertainty in the design | Beta design priors of the example | Uniform | Exact enumeration | validation/bridge_frequentist_bayesian.R | 831 | 0.8004 |  |
| 11 | Bayesian assurance, informative Sp prior | Beta design priors (Se Beta(17, 3), Sp Beta(18, 2), prevalence Beta(4, 16)); analysis = design | Beta(17, 3), Beta(18, 2) | Exact enumeration | bam_sample_size(method = "exact") | 672 | 0.8003 |  |
| 12 | Bayesian assurance, worked example with vague Sp prior | Beta design priors of the worked example (vague Sp prior) | Equal-tailed credible, analysis = design | Exact enumeration | bam_assurance_curve() | 678 | 0.8003 |  |
