# '.data' is the tidy-evaluation pronoun used inside ggplot2::aes() in the
# plotting functions (R/plots.R). ggplot2 re-exports it from rlang, but
# ggplot2 is a suggested package rather than an import (the package's hard
# dependencies are limited to base R), so it cannot be brought in with
# @importFrom. Declaring it here with utils::globalVariables() is the
# standard way to silence the resulting "no visible binding for global
# variable '.data'" NOTE from R CMD check without adding ggplot2/rlang to
# Imports.
utils::globalVariables(".data")
