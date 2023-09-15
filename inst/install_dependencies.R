# install_dependencies.R

# Install gGnome package from GitHub if CPLEX was not instantiated before.
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("mskilab-org/gGnome", force=TRUE)

