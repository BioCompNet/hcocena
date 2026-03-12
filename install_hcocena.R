install.packages("remotes")

remotes::install_github(
  "BioCompNet/hcocena",
  subdir = "hCoCena-r-package",
  dependencies = TRUE,
  upgrade = "never"
)
