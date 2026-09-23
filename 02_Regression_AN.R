##==============================================================================
## ILLUSTRATIVE EXAMPLE: 
##
## Input:  01_prepared_data_allsex_allage_allcause_iSEN.0_BRA_example.RData
##         (from the single-country patch to 01_Prepare_Data.R)
## Output: BRA_example_RegionLevel_AN.csv/.xlsx
##         BRA_example_CountryTotal_AN.csv/.xlsx
##==============================================================================

library(dplyr)
library(tidyr)
library(sf)
library(metafor)
library(matrixStats)
library(writexl)
library(conflicted)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::lag)

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), paste(..., collapse = " "), "\n")

##------------------------------------------------------------------------
## 0. Load the Brazil-only prepared data, config
##------------------------------------------------------------------------

# load("01_prepared_data_allsex_allage_allcause_iSEN.0_BRA_example.RData")

KEEP_TYPES <- c("Total Heat", "Total Cold")
cf_values  <- c(-2.5, 0, 2.5)   # overrides Part 1's 7-value set, for a
                                  # cleaner illustrative comparison
cf_labels  <- c("-2.5" = "La Nina (-2.5)", "0" = "Neutral (0)", "2.5" = "El Nino (+2.5)")

K_SUBSAMPLE_SIMS <- 50   # simulations to actually refit the spatial
                           # model on (region subsampling is NOT needed
                           # at this scale -- see header)
RHO_INIT <- 500

df_enso2 <- df_enso2 %>% filter(TYPE %in% KEEP_TYPES)

log_msg(sprintf(
  "Brazil example: %d regions, %d rows, TYPE = %s.",
  length(unique(df_enso2$REGION)), nrow(df_enso2), paste(KEEP_TYPES, collapse = ", ")
))

## Group indices, rebuilt for the (already filtered, but re-filtered
## here defensively by TYPE) dataset
df_enso2 <- df_enso2 %>% arrange(REGION, TYPE, SEASON_REG, YEAR)
grp_key <- paste(df_enso2$REGION, df_enso2$TYPE, df_enso2$SEASON_REG)
r <- rle(grp_key)
ngroups <- length(r$lengths)
group_end   <- cumsum(r$lengths)
group_start <- group_end - r$lengths + 1

n_years <- length(unique(df_enso2$YEAR))

get_predictors <- function(season) {
  if (season %in% c("JJA_before_ENSO", "SON_before_ENSO")) c("Nino34_0", "Nino34_prev")
  else if (season == "MAM") c("Nino34_0", "Nino34_next")
  else "Nino34_0"
}

##==============================================================================
## PART A: SPATIAL STRUCTURE
##==============================================================================

set.seed(20260914)
sim_subsample <- sort(sample.int(nsims, K_SUBSAMPLE_SIMS))

extract_region_coef_var <- function(df, type_sel, season_sel, response_col) {

  dat_all <- df %>% filter(TYPE == type_sel, SEASON_REG == season_sel)
  regions <- unique(dat_all$REGION)
  form <- reformulate(get_predictors(season_sel))
  out <- vector("list", length(regions))

  for (i in seq_along(regions)) {
    dat <- dat_all %>% filter(REGION == regions[i])
    X <- tryCatch(model.matrix(form, dat), error = function(e) NULL)
    if (is.null(X)) { out[[i]] <- NULL; next }

    y <- dat[[response_col]]
    npar <- ncol(X)
    df_resid <- nrow(X) - npar
    if (df_resid <= 0) { out[[i]] <- NULL; next }

    qrX <- qr(X)
    b <- qr.coef(qrX, y)
    if (any(is.na(b))) { out[[i]] <- NULL; next }

    resid <- y - X %*% b
    sigma2 <- sum(resid^2) / df_resid
    XtX_inv <- chol2inv(qrX$qr[seq_len(npar), , drop = FALSE])
    nino_row <- which(colnames(X) == "Nino34_0")

    out[[i]] <- data.frame(REGION = regions[i], beta = b[nino_row], var_beta = XtX_inv[nino_row, nino_row] * sigma2)
  }

  dplyr::bind_rows(out)

}

build_region_coef_matrix_all_sims <- function(df, type_sel, season_sel) {

  dat_all <- df %>% filter(TYPE == type_sel, SEASON_REG == season_sel)
  regions <- unique(dat_all$REGION)
  form <- reformulate(get_predictors(season_sel))
  coef_list <- vector("list", length(regions))
  names(coef_list) <- regions

  for (i in seq_along(regions)) {
    dat <- dat_all %>% filter(REGION == regions[i])
    X <- tryCatch(model.matrix(form, dat), error = function(e) NULL)
    if (is.null(X)) next
    Y <- as.matrix(dat[, sim_cols])
    if (nrow(X) - ncol(X) <= 0) next

    qrX <- qr(X)
    coef_mat <- qr.coef(qrX, Y)
    if (any(is.na(coef_mat))) next

    nino_row <- which(rownames(coef_mat) == "Nino34_0")
    if (length(nino_row) == 0) next
    coef_list[[i]] <- coef_mat[nino_row, ]
  }

  coef_list <- coef_list[!vapply(coef_list, is.null, logical(1))]
  if (length(coef_list) == 0) return(NULL)

  coef_mat_all <- do.call(cbind, coef_list)
  colnames(coef_mat_all) <- names(coef_list)
  coef_mat_all

}

get_gadm_centroids_moll <- function(polygon_sf) {
  s2_was_on <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_was_on), add = TRUE)

  n_invalid <- sum(!sf::st_is_valid(polygon_sf))
  if (n_invalid > 0) polygon_sf <- sf::st_make_valid(polygon_sf)

  polygon_moll <- sf::st_transform(polygon_sf, crs = "+proj=moll")
  cent_moll <- sf::st_centroid(polygon_moll)
  coords <- sf::st_coordinates(cent_moll)

  data.frame(REGION = cent_moll$GADM_ID, x_km = coords[, 1] / 1000, y_km = coords[, 2] / 1000)
}

fit_spatial_once <- function(coef_var_df, coord_df, dmat) {

  dat <- coef_var_df %>%
    inner_join(coord_df, by = "REGION") %>%
    filter(!is.na(beta), !is.na(var_beta), var_beta > 0) %>%
    mutate(const = "ALL")

  if (nrow(dat) < 10) return(NULL)   # lower threshold than production (30), given Brazil's small region count

  fit <- tryCatch(
    metafor::rma.mv(
      beta, var_beta, random = ~ REGION | const, struct = "SPGAU",
      data = dat, dist = list(dmat), control = list(rho.init = RHO_INIT)
    ),
    error = function(e) { cat(sprintf("    fit failed: %s\n", conditionMessage(e))); NULL }
  )
  if (is.null(fit)) return(NULL)
  
  ## Sanity cap on rho.
  ## When tau2 is very small, rho becomes poorly identified and REML can
  ## converge to an absurd, effectively-infinite range. Cap at a large-
  ## but-physically-plausible multiple of THIS fitting sample's own
  ## spatial extent, rather than trusting an unbounded estimate.
  max_extent <- max(dmat, na.rm = TRUE)
  rho_cap <- 10 * max_extent
  rho_fitted <- fit$rho
  
  if (!is.finite(rho_fitted) || rho_fitted > rho_cap) {
    cat(sprintf(
      "    rho (%.1f km) implausibly large relative to sample extent (%.1f km); capping at %.1f km.\n",
      rho_fitted, max_extent, rho_cap
    ))
    rho_fitted <- rho_cap
  }

  list(tau2 = fit$tau2, rho = rho_fitted, mu_hat = as.numeric(fit$b))

}

gadm_coords <- get_gadm_centroids_moll(polygon)

combos <- df_enso2 %>% distinct(TYPE, SEASON_REG)
spatial_lookup <- list()

log_msg("=== Part A: spatial structure ===")

for (ci in seq_len(nrow(combos))) {

  type_i <- combos$TYPE[ci]
  season_i <- combos$SEASON_REG[ci]

  cat(sprintf("[%d/%d] %s / %s\n", ci, nrow(combos), type_i, season_i))

  region_coef_matrix <- tryCatch(
    build_region_coef_matrix_all_sims(df_enso2, type_i, season_i),
    error = function(e) NULL
  )
  if (is.null(region_coef_matrix) || ncol(region_coef_matrix) < 10) {
    cat("  -> too few regions with valid coefficients, skipping.\n")
    next
  }

  all_regions <- colnames(region_coef_matrix)   # ALL Brazil regions, no subsampling
  coords_all <- gadm_coords %>% filter(REGION %in% all_regions)
  D_all <- as.matrix(dist(coords_all[, c("x_km", "y_km")]))
  rownames(D_all) <- colnames(D_all) <- coords_all$REGION

  param_draws <- vector("list", length(sim_subsample))

  for (j in seq_along(sim_subsample)) {
    sim_col <- paste0("sim_", sim_subsample[j])
    coef_var_j <- tryCatch(extract_region_coef_var(df_enso2, type_i, season_i, sim_col), error = function(e) NULL)
    if (is.null(coef_var_j) || nrow(coef_var_j) < 10) { param_draws[[j]] <- NULL; next }
    regs_j <- coef_var_j$REGION
    dmat_j <- D_all[regs_j, regs_j, drop = FALSE]
    param_draws[[j]] <- fit_spatial_once(coef_var_j, coords_all, dmat_j)
  }

  param_draws <- param_draws[!vapply(param_draws, is.null, logical(1))]

  if (length(param_draws) < 5) {
    cat(sprintf("  -> too few successful fits (%d), skipping.\n", length(param_draws)))
    next
  }

  tau2_vec <- vapply(param_draws, function(x) x$tau2, numeric(1))
  rho_vec  <- vapply(param_draws, function(x) x$rho, numeric(1))
  cat(sprintf(
    "  k=%d fits. tau2: mean=%.4g (CV=%.2f) | rho: mean=%.1f km (CV=%.2f)\n",
    length(param_draws), mean(tau2_vec), sd(tau2_vec) / mean(tau2_vec),
    mean(rho_vec), sd(rho_vec) / mean(rho_vec)
  ))

  ## Apply at full resolution (== all Brazil regions, same set used for fitting)
  n_reg <- length(all_regions)
  L_list <- vector("list", length(param_draws))

  for (k in seq_along(param_draws)) {
    p <- param_draws[[k]]
    R_spatial <- exp(-(D_all^2) / (p$rho^2))
    Sigma_k <- p$tau2 * R_spatial
    diag(Sigma_k) <- p$tau2
    L_list[[k]] <- tryCatch(t(chol(Sigma_k + diag(1e-8, n_reg))), error = function(e) NULL)
  }

  valid_k <- which(!vapply(L_list, is.null, logical(1)))
  if (length(valid_k) == 0) { cat("  -> no valid Cholesky factors, skipping.\n"); next }

  assign_k <- sample(valid_k, nsims, replace = TRUE)
  beta_sim <- matrix(NA_real_, nrow = nsims, ncol = n_reg)
  colnames(beta_sim) <- all_regions

  for (s in seq_len(nsims)) {
    k <- assign_k[s]
    Z <- rnorm(n_reg)
    beta_sim[s, ] <- region_coef_matrix[s, all_regions] + as.numeric(L_list[[k]] %*% Z)
  }

  spatial_lookup[[paste(type_i, season_i)]] <- beta_sim

}

log_msg(sprintf("Part A complete: spatial correction available for %d of %d strata.", length(spatial_lookup), nrow(combos)))

# You may see a metafor warning about extreme sampling-variance ratios during Part A. 
# This reflects heterogeneity in data quality across regions (some Brazilian states have shorter or less variable mortality records than others) 
# and was directly tested in the full-scale analysis. 
# Flooring the affected variances did not meaningfully change results, 
# so this warning can be treated as expected rather than indicating a specific fixable problem.

brazil_heat_jja <- df_enso2 %>%
  filter(TYPE == "Total Heat", SEASON_REG == "JJA_before_ENSO") %>%
  group_by(REGION) %>%
  summarise(
    sd_att_val = sd(att_val, na.rm = TRUE),
    n_distinct = length(unique(round(att_val, 8)))
  ) %>%
  arrange(sd_att_val)

print(brazil_heat_jja, n = 30)

#Near-zero heat-attributable mortality variance in JJA is concentrated in two distinct groups of Brazilian states, 
# driven by different climatic mechanisms: (1) southern subtropical states (Rio Grande do Sul, Santa Catarina, Paraná), 
# where winter temperatures are genuinely too cool for meaningful heat-attributable mortality, 
# and (2) coastal states with oceanic thermal moderation (e.g., Alagoas, Rio de Janeiro), 
# where proximity to the ocean narrows the temperature distribution overall. 
# Conversely, the highest variance is found in interior Amazonian and central-west states (Rondônia, Amazonas, Acre, Tocantins, Pará), 
# consistent with episodic cold-air incursion events 
# punctuating otherwise consistently hot conditions and producing substantial year-to-year variability in heat-related mortality.

##==============================================================================
## PART B: predict at 3 ENSO scenarios, aggregate to region and
## country (Brazil) total level
##==============================================================================

log_msg("=== Part B: prediction and aggregation ===")

region_key <- expand.grid(
  REGION = unique(df_enso2$REGION), TYPE = KEEP_TYPES, SEASON_REG = unique(df_enso2$SEASON_REG),
  YEAR = unique(df_enso2$YEAR), cf = cf_values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
region_key$.key <- paste(region_key$REGION, region_key$TYPE, region_key$SEASON_REG, region_key$YEAR, region_key$cf)
region_lookup <- setNames(seq_len(nrow(region_key)), region_key$.key)

country_key <- expand.grid(
  TYPE = KEEP_TYPES, SEASON_REG = unique(df_enso2$SEASON_REG),
  YEAR = unique(df_enso2$YEAR), cf = cf_values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
country_key$.key <- paste(country_key$TYPE, country_key$SEASON_REG, country_key$YEAR, country_key$cf)
country_lookup <- setNames(seq_len(nrow(country_key)), country_key$.key)

region_sum <- numeric(nrow(region_key)); region_sim <- matrix(0, nrow(region_key), nsims)
country_sum <- numeric(nrow(country_key)); country_sim <- matrix(0, nrow(country_key), nsims)

n_spatial_hits <- 0L; n_spatial_misses <- 0L

for (g in seq_len(ngroups)) {

  idx <- group_start[g]:group_end[g]
  dat <- df_enso2[idx, , drop = FALSE]

  form <- reformulate(get_predictors(dat$SEASON_REG[1]))
  X <- model.matrix(form, dat)
  Y <- as.matrix(dat[, c("att_val", sim_cols), drop = FALSE])

  qrX <- qr(X)
  coef_mat <- qr.coef(qrX, Y)
  coef_mat[is.na(coef_mat)] <- 0

  strat_key <- paste(dat$TYPE[1], dat$SEASON_REG[1])
  region_id <- dat$REGION[1]
  nino_idx <- which(rownames(coef_mat) == "Nino34_0")

  coef_full <- coef_mat[, -1, drop = FALSE]
  colnames(coef_full) <- sim_cols

  if (!is.null(spatial_lookup[[strat_key]]) && region_id %in% colnames(spatial_lookup[[strat_key]])) {
    n_spatial_hits <- n_spatial_hits + 1L
    coef_full[nino_idx, ] <- spatial_lookup[[strat_key]][, region_id]
  } else {
    n_spatial_misses <- n_spatial_misses + 1L
  }

  X_cf <- lapply(cf_values, function(cf) { nd <- dat; nd$Nino34_0 <- cf; model.matrix(form, nd) })

  for (j in seq_along(cf_values)) {

    pred_point <- X_cf[[j]] %*% coef_mat[, 1]
    pred_sim   <- X_cf[[j]] %*% coef_full

    att_num     <- pred_point[, 1] / 100 * dat$adjusted_factor
    att_num_sim <- sweep(pred_sim, 1, dat$adjusted_factor / 100, "*")

    region_rows <- region_lookup[paste(region_id, dat$TYPE, dat$SEASON_REG, dat$YEAR, cf_values[j])]
    region_sum[region_rows] <- region_sum[region_rows] + att_num
    region_sim[region_rows, ] <- region_sim[region_rows, ] + att_num_sim

    country_rows <- country_lookup[paste(dat$TYPE, dat$SEASON_REG, dat$YEAR, cf_values[j])]
    country_sum[country_rows] <- country_sum[country_rows] + att_num
    country_sim[country_rows, ] <- country_sim[country_rows, ] + att_num_sim

  }

}

log_msg(sprintf(
  "Spatial coverage: %d hits, %d misses (%.1f%%).",
  n_spatial_hits, n_spatial_misses, 100 * n_spatial_hits / (n_spatial_hits + n_spatial_misses)
))

##==============================================================================
## PART C: FORMAT AND SAVE
##==============================================================================

format_table <- function(key_df, sum_vec, sim_mat, group_vars) {
  key_df %>%
    mutate(
      att_num     = sum_vec / n_years,
      att_num_low = matrixStats::rowQuantiles(sim_mat, probs = 0.025) / n_years,
      att_num_upp = matrixStats::rowQuantiles(sim_mat, probs = 0.975) / n_years,
      ENSO = cf_labels[as.character(cf)]
    ) %>%
    group_by(across(all_of(c(group_vars, "TYPE", "ENSO")))) %>%
    summarise(att_num = sum(att_num), att_num_low = sum(att_num_low), att_num_upp = sum(att_num_upp), .groups = "drop") %>%
    mutate(AN_formatted = sprintf("%s (%s; %s)",
                                    format(round(att_num), big.mark = ","),
                                    format(round(att_num_low), big.mark = ","),
                                    format(round(att_num_upp), big.mark = ",")))
}

region_table <- format_table(region_key, region_sum, region_sim, "REGION")
country_table <- format_table(country_key, country_sum, country_sim, character(0))

cat("\n=== Brazil country-total AN, by TYPE and ENSO scenario ===\n")
print(as.data.frame(country_table), digits = 4)

writexl::write_xlsx(region_table, "BRA_example_RegionLevel_AN.xlsx")
writexl::write_xlsx(country_table, "BRA_example_CountryTotal_AN.xlsx")

log_msg("Illustrative Brazil pipeline complete.")
