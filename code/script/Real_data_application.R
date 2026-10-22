# Application to the Brazilian homicide data

#------ Packages and paths ----------------------------------------------------
repos <- getOption("repos")
if (is.null(repos) || isTRUE(repos["CRAN"] == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(sf, dplyr, spdep, Matrix, ggplot2, ggspatial, scales, MASS, here)
source(here("script", "functions.R"))
install_inla()

#------ Command-line arguments ------------------------------------------------
options(stringsAsFactors = FALSE, scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
data_file <- arg(1, Sys.getenv("data", "Municipios_Brasil.RData"))
impute <- arg(3, "complete")
join <- as.logical(arg(4, "FALSE"))
shapes <- as.logical(arg(5, "FALSE"))
redraw <- arg(6, "none")
out_dir <- arg(2, file.path("output", "brazil",
                            paste0(impute, if (join) "_joined" else "")))
if (!grepl("^(/|[A-Za-z]:)", out_dir)) out_dir <- here(out_dir)
threads <- max(1L, parallel::detectCores() - 1L)
alpha <- 0.05

#------ Setup -----------------------------------------------------------------
set.seed(20260725)
cfg <- list(out_dir = out_dir)
ensure_dir(out_dir)
save_csv <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, name), row.names = FALSE)
}
save_workspace <- function(res = NULL) {
  saveRDS(list(sensitivity = res, inla = inla_out, impute = impute),
          file.path(out_dir, "analysis_workspace.rds"), compress = "xz")
}

if (redraw != "none") {
  ws <- readRDS(file.path(out_dir, "analysis_workspace.rds"))
  if (redraw == "maps") {
    dat <- prepare_brazil(data_file, impute = impute, join = join)$data
    fig_app_maps(dat, ws$inla$ref_fit, out_dir)
    fig_app_effects(ws$inla$irr, out_dir)
  } else {
    if (is.null(ws$sensitivity)) stop("the run has no sensitivity results")
    write_sensitivity(ws$sensitivity, ws$inla, cfg, alpha)
  }
  quit(save = "no")
}

prep <- prepare_brazil(data_file, impute = impute, join = join)
dat <- prep$data
spec <- brazil_spectrum(prep$adj, spectrum_cache(out_dir, prep))
spec$label <- "brazil"
gs <- graph_summary(spec)
msg(gs$n, " municipalities, ", gs$components, " components (",
    gs$isolated, " isolated), ", sum(dat$homicides), " homicides")
save_csv(gs, "graph_summary.csv")
write_table(prep$flow, "table_sample_flow", cfg,
            "Municipalities kept at each step of the selection.",
            "tab:flow")
write_table(prep$missing, "table_missing", cfg,
            "Missing values in the raw covariates, before any handling.",
            "tab:missing")
write_table(descriptive_table(dat), "table_descriptive", cfg,
            paste("Median (first, third quartile) of the outcome and the",
                  "covariates, overall and by region."),
            "tab:descriptive", size = "\\scriptsize")

inla_out <- NULL
if (requireNamespace("INLA", quietly = TRUE)) {
  inla_out <- run_inla_models(dat, prep$adj, out_dir, threads)
  save_workspace()
  write_table(irr_table(inla_out$irr), "table_inla_irr", cfg,
              paste("Rate ratios per standard deviation (95% credible",
                    "intervals) from the BYM2 model with and without the",
                    "spatial adjustments."), "tab:inla")
  write_table(waic_table(inla_out$waic), "table_waic", cfg,
              paste("Widely applicable information criterion of the",
                    "unrestricted BYM2 model by likelihood."), "tab:waic")
  write_table(check_table_ppc(inla_out$checks), "table_predictive", cfg,
              paste("Observed and posterior expected number of zeros and of",
                    "high counts, and Pearson dispersion at the fitted",
                    "means, by likelihood."), "tab:predictive")
  fig_app_maps(dat, inla_out$ref_fit, out_dir)
  fig_app_effects(inla_out$irr, out_dir)
} else {
  msg("INLA unavailable: Bayesian fits skipped")
}

# The BYM2 mixing parameter fixes the shape of the class, and the fitted
# scale and dispersion define the spatial fit that is compared with the
# non-spatial one.
phi_ref <- if (is.null(inla_out)) 0.5 else inla_out$hyper[["phi"]]
theta <- if (!is.null(inla_out) && inla_out$family == "nbinomial") {
  inla_out$hyper[["theta"]]
} else Inf
sigma <- if (is.null(inla_out)) NULL else inla_out$hyper[["sigma"]]
qinv <- shape_qinv(spec, "bym2", phi_ref)
res <- run_sensitivity(dat, spec, qinv, alpha = alpha, theta = theta,
                       sigma = sigma)

save_workspace(res)
write_sensitivity(res, inla_out, cfg, alpha)
if (shapes) {
  shape_res <- run_shapes(dat, spec, alpha = alpha, theta = theta)
  save_csv(shape_res, "sensitivity_shapes.csv")
  write_table(shape_table(shape_res), "table_shapes", cfg,
              paste("Information fraction, break-even and critical radii",
                    "relative to $c_{\\mathrm{ref}}$ for other shapes of the",
                    "class."), "tab:shapes")
}

writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(out_dir, "session_info.txt"))
save_workspace(res)
msg("finished")
