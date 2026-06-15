# External cross-validation of dtasamplesize core formulas
# ---------------------------------------------------------
# Verifies the numerical machinery against INDEPENDENT references:
# Hmisc, base prop.test, pROC, numerical integration, and Monte Carlo.
# Run after installing the package:  R -f validation/cross_validation.R
# Optional reference packages: Hmisc, pROC.

library(dtasamplesize)
suppressMessages({library(Hmisc); library(pROC)})
options(width = 110, digits = 6)
ok <- function(x) if (isTRUE(x)) "PASS" else "**CHECK**"
sep <- function(t) cat("\n========== ", t, " ==========\n")

## ---- V1. Wilson CI vs Hmisc::binconf(wilson) and base prop.test ----
sep("V1. Wilson score CI")
cases <- list(c(85,100), c(45,50), c(2,20), c(19,20), c(0,30))
for (cs in cases) {
  x <- cs[1]; n <- cs[2]
  mine <- wilson_ci(x, n)
  hm   <- binconf(x, n, method = "wilson")            # Hmisc
  pt   <- prop.test(x, n, correct = FALSE)$conf.int    # base R = Wilson score
  cat(sprintf("x=%3d n=%3d | mine[%.5f,%.5f] Hmisc[%.5f,%.5f] prop.test[%.5f,%.5f]  %s / %s\n",
      x, n, mine["lower"], mine["upper"], hm[2], hm[3], pt[1], pt[2],
      ok(all(abs(c(mine["lower"]-hm[2], mine["upper"]-hm[3])) < 1e-6)),
      ok(all(abs(c(mine["lower"]-pt[1], mine["upper"]-pt[2])) < 1e-6))))
}
cat("Note: at x=19/n=20 the package matches base prop.test (canonical Wilson);\n",
    "Hmisc differs at that boundary.\n")

## ---- V2. Wald CI vs Hmisc::binconf(asymptotic) ----
sep("V2. Wald CI (interior cases; wald_ci clamps to [0,1])")
for (cs in list(c(85,100), c(45,50), c(60,90))) {
  x <- cs[1]; n <- cs[2]
  mine <- wald_ci(x, n)
  hm   <- binconf(x, n, method = "asymptotic")
  cat(sprintf("x=%3d n=%3d | mine[%.5f,%.5f] Hmisc[%.5f,%.5f]  %s\n",
      x, n, mine["lower"], mine["upper"], hm[2], hm[3],
      ok(all(abs(c(mine["lower"]-hm[2], mine["upper"]-hm[3])) < 1e-6))))
}

## ---- V3. buderer_n known published values ----
sep("V3. buderer_n")
cat("Se=0.85,d=0.07 ->", buderer_n(0.85,0.07), ok(buderer_n(0.85,0.07)==100), "\n")
cat("Se=0.90,d=0.05 ->", buderer_n(0.90,0.05), ok(buderer_n(0.90,0.05)==139), "\n")

## ---- V4. Hanley-McNeil variance: formula vs Monte Carlo empirical AUC var ----
sep("V4. Hanley-McNeil AUC variance vs Monte Carlo (binormal data)")
emp_auc_var <- function(AUC, n_case, n_ctrl, R = 4000, seed = 7) {
  set.seed(seed)
  mu <- qnorm(AUC) * sqrt(2)        # binormal: AUC = pnorm(mu/sqrt2)
  a <- numeric(R)
  for (r in seq_len(R)) {
    xc <- rnorm(n_case, mu, 1); yc <- rnorm(n_ctrl, 0, 1)
    rk <- rank(c(xc, yc))
    U <- sum(rk[seq_len(n_case)]) - n_case*(n_case+1)/2
    a[r] <- U / (n_case * n_ctrl)
  }
  var(a)
}
for (cfg in list(c(0.80,100,100), c(0.85,60,140), c(0.75,80,80))) {
  AUC <- cfg[1]; nc <- cfg[2]; nk <- cfg[3]
  f  <- hanley_mcneil_var(AUC, nc, nk)
  e  <- emp_auc_var(AUC, nc, nk)
  cat(sprintf("AUC=%.2f nC=%3d nK=%3d | formula=%.6f  MC=%.6f  ratio=%.3f  %s\n",
      AUC, nc, nk, f, e, f/e, ok(abs(f/e - 1) < 0.20)))
}
cat("Note: Hanley-McNeil (1982) assumes a negative-exponential score model;\n",
    "the small deviation from binormal MC (larger when groups are imbalanced)\n",
    "is the expected, conservative behaviour of the approximation.\n")

## ---- V5. Net benefit point estimate vs standard Vickers definition ----
sep("V5. Net benefit identity: prev-weighted form == (TP - w*FP)/N")
set.seed(1); Se<-0.85; Sp<-0.90; prev<-0.20; pt<-0.20; w<-pt/(1-pt); N<-300
n_d<-floor(N*prev); n_nd<-N-n_d
TP<-rbinom(1,n_d,Se); FP<-rbinom(1,n_nd,1-Sp)
nb_pkgform <- (n_d/N)*(TP/n_d) - (n_nd/N)*(FP/n_nd)*w     # package algebra
nb_vickers <- (TP - w*FP)/N                                # Vickers 2006 def
cat(sprintf("package=%.6f  Vickers=%.6f  %s\n", nb_pkgform, nb_vickers,
    ok(abs(nb_pkgform - nb_vickers) < 1e-9)))

## ---- V6. Net benefit variance vs Monte Carlo SD ----
sep("V6. Net benefit analytic SE vs Monte Carlo SD of NB_hat")
nb_se_test <- function(Se, Sp, prev, pt, N, R = 20000, seed = 3) {
  w <- pt/(1-pt); n_d <- floor(N*prev); n_nd <- N - n_d
  wd <- n_d/N; wnd <- n_nd/N
  set.seed(seed)
  TP <- rbinom(R, n_d, Se); FP <- rbinom(R, n_nd, 1-Sp)
  Se_h <- TP/n_d; Sp_h <- 1 - FP/n_nd
  NB <- wd*Se_h - wnd*(1-Sp_h)*w
  mc_sd <- sd(NB)
  analytic <- sqrt(wd^2*Se*(1-Se)/n_d + wnd^2*w^2*Sp*(1-Sp)/n_nd)  # plug-in truth
  c(analytic = analytic, mc_sd = mc_sd, ratio = analytic/mc_sd)
}
for (cfg in list(c(0.85,0.90,0.20,0.20,300), c(0.80,0.85,0.30,0.40,200))) {
  r <- do.call(nb_se_test, as.list(cfg))
  cat(sprintf("Se=%.2f Sp=%.2f prev=%.2f pt=%.2f N=%d | analytic=%.5f MC=%.5f ratio=%.3f %s\n",
      cfg[1],cfg[2],cfg[3],cfg[4],cfg[5], r["analytic"], r["mc_sd"], r["ratio"],
      ok(abs(r["ratio"]-1) < 0.05)))
}

## ---- V7. BAM Beta-Binomial conjugacy vs numerical Bayes posterior ----
sep("V7. Beta-Binomial posterior CI (conjugate) vs numerical integration")
a<-17; b<-3; n<-100; x<-86
cw <- qbeta(0.975, a+x, b+n-x) - qbeta(0.025, a+x, b+n-x)   # conjugate (BAM)
g <- seq(1e-5, 1-1e-5, length.out = 200001)                # numerical grid
post <- dbeta(g,a,b) * dbinom(x, n, g); post <- post/sum(post)
cdf <- cumsum(post)
lo <- g[which.min(abs(cdf-0.025))]; hi <- g[which.min(abs(cdf-0.975))]
cat(sprintf("conjugate width=%.5f  numerical width=%.5f  %s\n",
    cw, hi-lo, ok(abs(cw-(hi-lo)) < 2e-3)))

## ---- V8. Imperfect-ref: VIF (Rogan-Gladen) and apparent Se closed form ----
sep("V8. Imperfect reference: VIF and apparent-Se closed form")
Se_ref<-0.90; Sp_ref<-0.95; Se<-0.85; Sp<-0.90; prev<-0.30
vif_pkg <- ss_imperfect_ref(Se=Se,Sp=Sp,Se_ref=Se_ref,Sp_ref=Sp_ref,prev=prev,B=0)$VIF
vif_rg  <- 1/(Se_ref+Sp_ref-1)^2
cat(sprintf("VIF package=%.6f  Rogan-Gladen=%.6f  %s\n", vif_pkg, vif_rg, ok(abs(vif_pkg-vif_rg)<1e-9)))
ir <- ss_imperfect_ref(Se=Se,Sp=Sp,Se_ref=Se_ref,Sp_ref=Sp_ref,prev=prev,B=6000,seed=11,sensitivity_table=FALSE)
app_emp <- ir$mc_validation$se_apparent[ir$mc_validation$scenario=="adjusted"]
app_theo <- (prev*Se*Se_ref + (1-prev)*(1-Sp)*(1-Sp_ref)) /
            (prev*Se_ref + (1-prev)*(1-Sp_ref))
cat(sprintf("apparent Se: MC=%.4f  closed-form=%.4f  %s\n", app_emp, app_theo,
    ok(abs(app_emp-app_theo) < 0.01)))

cat("\n=== cross_validation done ===\n")
