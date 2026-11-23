| Capability | presize (0.3.11) | pROC (1.19.1) | MKpower (1.1) | epiR (3.0.1) | MKmisc (2.0) | SampleSizeDiagnostics (0.1.1) | dtasamplesize |
|---|---|---|---|---|---|---|---|
| Buderer-type Se+Sp sizing via prevalence (max of two arms) | No | No | No | Yes (epi.ssdxsesp) | No | Yes | Yes (buderer_n) |
| CI-width (precision) sizing, single arm (Se or Sp) | Yes | No | Yes (ssize.sens.ci/ssize.spec.ci) | Yes (epi.ssdxsesp) | Yes (power.diagnostic.test) | Yes | Yes |
| AUC precision/power sizing | Yes (prec_auc) | Yes (power.roc.test) | Yes (ssize.auc.ci) | No | No | No | Partly (AUC gate in joint_sample_size(); AUC modules experimental) |
| Frequentist power/assurance on target width, single arm | No | No (power.roc.test targets AUC/significance, not a Se/Sp width target) | Yes (power.diagnostic.test/ssize.*.ci, method = "exact") | No | Yes (power.diagnostic.test) | No | Partly (mc_validate_buderer() evaluates it at a given n; joint_sample_size() searches both arms jointly) |
| Joint Se+Sp+AUC precision, one simultaneous search | No | No (AUC only) | No (three separate one-arm calls) | No | No | No | Yes, Se and Sp with a deterministic AUC gate (joint_sample_size()); combined module experimental |
| Imperfect reference standard correction | No | No | No | Yes (epi.ssdxtest, Hui-Walter design, two populations and two tests) | No | No | Experimental |
| Bayesian assurance (prior on Se/Sp/prevalence) | No | No | No | No | No | No | Yes (validated) |
