# Version guard for tests that pin exact Monte Carlo figures ----------------
#
# Some tests pin figures (for example 0.8041, 0.573, 240, 920, 0.800112)
# captured under published R releases. A change in those figures was
# observed on a development build of R (R-devel, "4.6.0 r90605, 2026-09-30";
# NEWS 0.6.7), and was not reproduced on the published release R 4.6.2.
#
# rng_pins_apply() is TRUE on published R releases and FALSE on development
# builds: R.version$status is "" for releases and, for example,
# "Under development (unstable)" for R-devel. It deliberately does not
# compare getRversion(), so published R 4.6.x releases are tested against the
# pinned values.
rng_pins_apply <- function() !nzchar(R.version$status)
