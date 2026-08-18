# Test suite setup ------------------------------------------------------
#
# Many tests deliberately use a small B (Monte Carlo replications) to keep
# the suite fast. Silence the "B is small" advisory warning (see
# warn_small_B() in R/helpers.R) for the whole suite; the warning itself is
# exercised explicitly, with the option restored to TRUE, in the tests that
# target it.
options(dtasamplesize.warn_small_B = FALSE)
