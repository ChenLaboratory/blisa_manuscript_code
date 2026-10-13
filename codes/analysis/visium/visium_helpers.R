# Visium preparation and CSV export helpers extracted from the executed manuscript analysis.
# Source after codes/helpers/source_data_helpers.R. Sourcing performs no analysis.


visium_position_table <- function(spe) {
  if (!inherits(spe, "SpatialExperiment")) stop("Expected a SpatialExperiment.")
  cd <- SummarizedExperiment::colData(spe)
  xy <- SpatialExperiment::spatialCoords(spe)
  fields <- c("in_tissue", "array_row", "array_col")
  coords <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  if (!all(fields %in% names(cd)) || !all(coords %in% colnames(xy))) {
    stop("Expected the vignette SPE's tissue, array and full-resolution pixel coordinates.")
  }
  ids <- colnames(spe)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Unique, non-empty spot barcodes are required.")
  }
  positions <- data.frame(barcode = ids, as.data.frame(cd[, fields]),
                         xy[, coords, drop = FALSE], row.names = NULL)
  # The example SPE's pixel-named columns have already been converted to
  # micrometres. Respect the recorded unit when recovering image pixels.
  unit <- S4Vectors::metadata(spe)$coord_unit
  if (identical(unit, "micron")) {
    um_per_pixel <- S4Vectors::metadata(spe)$um_per_pixel
    if (length(um_per_pixel) != 1L || !is.numeric(um_per_pixel) ||
        !is.finite(um_per_pixel) || um_per_pixel <= 0) {
      stop("Micrometre coordinates require a positive metadata$um_per_pixel.")
    }
    positions[, coords] <- positions[, coords] / um_per_pixel
  } else if (!identical(unit, "pixel")) {
    stop("Expected SPE metadata$coord_unit to be 'micron' or 'pixel'.")
  }
  if (any(!is.finite(as.matrix(positions[, -1])))) stop("Invalid spot coordinates.")
  positions
}


visium_cell_type_palette <- function() {
  # Category order and colours shared with the Xenium manuscript report.
  c(Stromal = "#4D4D4D", Invasive_Tumor = "#377EB8", DCIS_1 = "#4DAF4A",
    DCIS_2 = "#006D2C", Macrophages = "#984EA3", Myoepi = "#FF7F00",
    Endothelial = "#FFFF33", `CD4+_T_Cells` = "#A65628", `CD8+_T_Cells` = "#66C2A5",
    Others = "#FC8D62", B_Cells = "#8DA0CB", `Perivascular-Like` = "#E78AC3",
    DCs = "#A6D854", Mast_Cells = "#FFD92F")
}


visium_cell_type_mapping <- function() {
  c(B_Cells = "B_Cells", `CD4+_T_Cells` = "CD4+_T_Cells", `CD8+_T_Cells` = "CD8+_T_Cells",
    DCIS_1 = "DCIS_1", DCIS_2 = "DCIS_2", Endothelial = "Endothelial",
    Invasive_Tumor = "Invasive_Tumor", `IRF7+_DCs` = "DCs", `LAMP3+_DCs` = "DCs",
    Macrophages_1 = "Macrophages", Macrophages_2 = "Macrophages", Mast_Cells = "Mast_Cells",
    `Myoepi_ACTA2+` = "Myoepi", `Myoepi_KRT15+` = "Myoepi", `Perivascular-Like` = "Perivascular-Like",
    Prolif_Invasive_Tumor = "Invasive_Tumor", Stromal = "Stromal",
    `Stromal_&_T_Cell_Hybrid` = "Others", `T_Cell_&_Tumor_Hybrid` = "Others")
}


harmonize_visium_cell_types <- function(estimates) {
  if (!is.matrix(estimates) || !is.numeric(estimates) || !ncol(estimates)) {
    stop("Cell-type estimates must be a numeric matrix with named columns.")
  }
  source_types <- colnames(estimates)
  if (is.null(source_types) || anyNA(source_types) || any(!nzchar(source_types)) ||
      anyDuplicated(source_types)) stop("Cell-type estimate columns must have unique non-empty names.")
  mapping <- visium_cell_type_mapping()
  unmapped <- setdiff(source_types, names(mapping))
  if (length(unmapped)) stop("Unmapped Visium cell types: ", paste(unmapped, collapse = ", "))
  # Entirely missing rows represent unmatched barcodes. Partially missing rows
  # cannot be interpreted as a complete deconvolution and are rejected.
  missing <- rowSums(is.na(estimates)) == ncol(estimates)
  observed <- estimates[!missing, , drop = FALSE]
  if (any(!is.finite(observed)) || any(observed < 0)) stop("Invalid cell-type estimates.")
  type_order <- names(visium_cell_type_palette())
  grouped <- matrix(0, nrow = nrow(estimates), ncol = length(type_order),
                    dimnames = list(rownames(estimates), type_order))
  for (type in type_order) {
    source <- source_types[unname(mapping[source_types]) == type]
    if (length(source)) grouped[, type] <- rowSums(estimates[, source, drop = FALSE])
  }
  grouped[missing, ] <- NA_real_
  total <- rowSums(estimates)
  known <- !missing & is.finite(total) & total > 0
  normalize <- function(x) {
    x[known, ] <- x[known, , drop = FALSE] / total[known]
    x[!known, ] <- NA_real_
    x
  }
  proportions <- normalize(grouped)
  proportions_original <- normalize(estimates)
  dominant <- dominant_original <- rep("Unknown", nrow(estimates))
  # Ties in the harmonized categories follow the shared Xenium category order;
  # the original labels retain the published table's column-order tie rule.
  dominant[known] <- type_order[max.col(proportions[known, , drop = FALSE], ties.method = "first")]
  dominant_original[known] <- source_types[max.col(proportions_original[known, , drop = FALSE],
                                                 ties.method = "first")]
  levels <- c(type_order, if (any(!known)) "Unknown")
  list(estimates = grouped, proportions = proportions,
       estimates_original = estimates, proportions_original = proportions_original,
       dominant = factor(dominant, levels = levels), dominant_original = dominant_original,
       total = total, known = known,
       category_mapping = data.frame(original_cell_type = source_types,
                                     xenium_cell_type = unname(mapping[source_types])))
}


lr_panel_coverage <- function(database, gene_panel) {
  required <- c("ligand.symbol", "receptor.symbol", "annotation", "pathway_name")
  if (!is.data.frame(database) || !all(required %in% names(database))) {
    stop("Supply the native BLISA interaction database with expanded partner symbols.")
  }
  ids <- rownames(database)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Database interactions must have unique, nonempty row identifiers.")
  }
  gene_panel <- unique(as.character(gene_panel))
  if (!length(gene_panel) || anyNA(gene_panel) || any(!nzchar(gene_panel))) {
    stop("The data gene panel must contain nonempty gene symbols.")
  }
  # Native parsing handles comma- or underscore-separated subunits and whitespace.
  parse_units <- getFromNamespace("parse_units", "blisa")
  ligand_units <- lapply(as.character(database$ligand.symbol), parse_units)
  receptor_units <- lapply(as.character(database$receptor.symbol), parse_units)
  valid_units <- function(x) length(x) > 0L && !anyNA(x) && all(nzchar(x))
  if (!all(vapply(c(ligand_units, receptor_units), valid_units, logical(1)))) {
    stop("Database contains missing ligand/receptor symbols; cannot determine panel coverage.")
  }
  missing_units <- function(units) {
    vapply(units, function(x) paste(setdiff(x, gene_panel), collapse = ", "), character(1))
  }
  out <- data.frame(
    database_id = ids,
    database[, required, drop = FALSE],
    missing_ligand_genes = missing_units(ligand_units),
    missing_receptor_genes = missing_units(receptor_units),
    stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL
  )
  out$ligand_covered <- !nzchar(out$missing_ligand_genes)
  out$receptor_covered <- !nzchar(out$missing_receptor_genes)
  # Use BLISA's presence-only selection, without the filterLRpairs count screen.
  included <- getFromNamespace("getLRpairs", "blisa")(
    gene_panel = gene_panel, LR_df = database
  )
  out$panel_covered <- ids %in% rownames(included)
  stopifnot(identical(out$panel_covered, out$ligand_covered & out$receptor_covered))
  out
}


lr_panel_coverage_summary <- function(coverage) {
  annotation <- as.character(coverage$annotation)
  annotation[is.na(annotation) | !nzchar(annotation)] <- "Unannotated"
  categories <- unique(annotation)
  groups <- c(list(seq_len(nrow(coverage))),
              lapply(categories, function(x) which(annotation == x)))
  n_database <- vapply(groups, length, integer(1))
  n_covered <- vapply(groups, function(i) sum(coverage$panel_covered[i]), integer(1))
  data.frame(
    annotation = c("All interactions", categories),
    database_pairs = n_database,
    covered_pairs = n_covered,
    coverage_percent = ifelse(n_database > 0L, 100 * n_covered / n_database, NA_real_),
    stringsAsFactors = FALSE
  )
}


native_report_tables <- function(result, spe) {
  if (is.null(result) || is.null(result$LR_results) || is.null(result$bins) || is.null(result$spatial_weights)) {
    stop("Expected a native BLISA result with LR_results, bins and spatial_weights.")
  }
  output <- as.data.frame(result$LR_results, check.names = FALSE)
  ids <- rownames(output)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Invalid native BLISA LR rownames.")
  cd <- as.data.frame(SummarizedExperiment::colData(spe), optional = TRUE)
  before_count_filter <- seq_len(ncol(spe))
  if ("in_tissue" %in% names(cd)) {
    if (anyNA(cd$in_tissue)) stop("Cannot map native Visium bins with missing in_tissue values.")
    before_count_filter <- before_count_filter[cd$in_tissue == 1]
  }
  bin_index <- result$bins$bin_id
  if (!is.numeric(bin_index) || anyNA(bin_index) || any(bin_index != as.integer(bin_index)) ||
      anyDuplicated(bin_index) || any(!bin_index %in% seq_along(before_count_filter))) {
    stop("Cannot map native Visium bin IDs back to prepared spots.")
  }
  # visiumSpotBins numbers bins after its tissue filter, before its native count filter.
  mapped <- before_count_filter[as.integer(bin_index)]
  expected <- before_count_filter[Matrix::colSums(SummarizedExperiment::assay(spe, "counts"))[
    before_count_filter] >= eval(formals(blisa::visiumSpotBins)$min_total_counts)]
  if (!identical(mapped, expected)) stop("Native Visium result bins do not match its default spot filters.")
  copy <- result
  copy$bins$spot_id <- colnames(spe)[mapped]
  output$lr_id <- ids
  copy$LR_results <- output
  scalar <- vapply(output, function(x) is.atomic(x) && is.null(dim(x)), logical(1))
  manifest <- output[, scalar, drop = FALSE]
  for (field in c("ligand.symbol", "receptor.symbol", "annotation", "pathway_name")) {
    if (!field %in% names(manifest)) manifest[[field]] <- rep(NA_character_, nrow(manifest))
  }
  manifest$is_complex <- grepl(",", manifest$ligand.symbol, fixed = TRUE) |
    grepl(",", manifest$receptor.symbol, fixed = TRUE)
  pair_status <- data.frame(lr_id = ids, status = rep("tested", length(ids)),
    reason = rep("Returned by the native BLISA Visium pipeline", length(ids)), stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(output)), function(i) {
    p <- output$all_pval[[i]]; score <- output$all_lisa[[i]]; q <- output$all_quadrant[[i]]
    n <- nrow(copy$bins)
    if (!is.numeric(p) || !is.numeric(score) || length(p) != n || length(score) != n || length(q) != n) {
      stop("Native BLISA per-spot output dimensions are inconsistent.")
    }
    mode <- as.character(output$ccc_mode[i])
    if (length(mode) != 1L || is.na(mode) || !mode %in% c("nearby", "diffuse")) stop("Unknown native BLISA spatial mode.")
    keep <- if (mode == "nearby") result$spatial_weights$keep_idx_queen else result$spatial_weights$keep_idx_dist
    if (is.null(keep)) stop("Missing native BLISA tested-spot mask.")
    if (anyNA(keep) || any(keep != as.integer(keep)) || any(!keep %in% seq_len(n))) stop("Invalid native tested-spot indices.")
    tested <- rep(FALSE, n); tested[keep] <- TRUE
    tested <- tested & is.finite(p)
    sig_index <- output$sig_index[[i]]
    if (anyNA(sig_index) || any(sig_index != as.integer(sig_index)) || any(!sig_index %in% seq_len(n))) {
      stop("Invalid native BLISA hotspot indices.")
    }
    sig <- rep(FALSE, n); sig[sig_index] <- TRUE
    if (any(sig & !tested)) stop("A native BLISA hotspot must be a tested spot.")
    p[!tested] <- NA_real_; score[!tested] <- NA_real_
    data.frame(lr_id = ids[i], spot_id = copy$bins$spot_id, p_value = p,
      adjusted_p_value = rep(NA_real_, n), score = score, tested = tested,
      raw_significant = sig, native_significant = sig, quadrant = as.character(q), stringsAsFactors = FALSE)
  })
  lr_spot <- do.call(rbind, rows)
  list(result = copy, manifest = manifest, pair_status = pair_status, lr_spot = lr_spot)
}


lri_result_summary <- function(result) {
  if (is.null(result$LR_results) || is.null(result$bins)) {
    stop("Expected a BLISA result with LR_results and bins.")
  }
  x <- as.data.frame(result$LR_results, check.names = FALSE)
  ids <- rownames(x)
  if (length(ids) != nrow(x) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("BLISA LR row names must be nonempty, unique interaction identifiers.")
  }
  if (anyDuplicated(names(x)) || !all(c("sig_numbers", "sig_index", "all_pval") %in% names(x))) {
    stop("LR_results must have unique column names and native hotspot counts, indices and P values.")
  }
  if ("lr_id" %in% names(x) && !identical(as.character(x$lr_id), ids)) {
    stop("Existing lr_id values differ from the native BLISA row names.")
  }
  if ("rank" %in% names(x) && !identical(as.numeric(x$rank), as.numeric(seq_len(nrow(x))))) {
    stop("Existing rank values differ from native BLISA row order.")
  }
  n_bins <- nrow(result$bins)
  if (length(n_bins) != 1L || !is.finite(n_bins) || n_bins < 0L) {
    stop("BLISA bins must have a valid row count.")
  }
  if (!is.numeric(x$sig_numbers) || anyNA(x$sig_numbers) ||
      any(!is.finite(x$sig_numbers)) || any(x$sig_numbers < 0) ||
      any(x$sig_numbers != as.integer(x$sig_numbers)) ||
      !is.list(x$sig_index) || !is.list(x$all_pval)) {
    stop("Invalid native BLISA hotspot counts or per-bin lists.")
  }
  valid_indices <- function(index) {
    is.numeric(index) && !anyNA(index) && all(is.finite(index)) &&
      all(index == as.integer(index)) && !anyDuplicated(index) &&
      all(index %in% seq_len(n_bins))
  }
  tested <- integer(nrow(x))
  for (i in seq_len(nrow(x))) {
    significant <- x$sig_index[[i]]
    if (!valid_indices(significant) || x$sig_numbers[i] != length(significant)) {
      stop("Inconsistent sig_numbers or sig_index for ", ids[i], ".")
    }
    p <- x$all_pval[[i]]
    if (!is.numeric(p) || length(p) != n_bins ||
        any(is.finite(p) & (p < 0 | p > 1))) {
      stop("Invalid per-bin P values for ", ids[i], ".")
    }
    available <- is.finite(p)
    # Native BLISA can leave finite placeholders in isolated bins. Intersect
    # with the spatial-mode mask so these do not inflate the tested count.
    if (!is.null(result$spatial_weights)) {
      mode <- if ("ccc_mode" %in% names(x)) as.character(x$ccc_mode[i]) else NA_character_
      if (is.na(mode) || !mode %in% c("nearby", "diffuse")) {
        stop("Missing or invalid native BLISA spatial mode for ", ids[i], ".")
      }
      key <- if (mode == "nearby") "keep_idx_queen" else "keep_idx_dist"
      keep <- result$spatial_weights[[key]]
      if (is.null(keep) || !valid_indices(keep)) {
        stop("Missing or invalid tested-bin mask ", key, " for ", ids[i], ".")
      }
      available <- available & seq_len(n_bins) %in% keep
    }
    if (any(!available[significant])) {
      stop("A significant hotspot is outside the tested bins for ", ids[i], ".")
    }
    tested[i] <- sum(available)
  }
  scalar <- vapply(x, function(column) is.atomic(column) && is.null(dim(column)), logical(1))
  output <- x[, scalar, drop = FALSE]
  for (key in c("n_tested", "n_total_bins")) {
    expected <- if (key == "n_tested") tested else rep(n_bins, nrow(x))
    if (key %in% names(output) && !identical(as.numeric(output[[key]]), as.numeric(expected))) {
      stop("Existing ", key, " conflicts with the per-bin BLISA results.")
    }
  }
  output$rank <- seq_len(nrow(x))
  output$lr_id <- ids
  output$n_tested <- tested
  output$n_total_bins <- rep(n_bins, nrow(x))
  output <- output[, c("rank", "lr_id", setdiff(names(output), c("rank", "lr_id"))), drop = FALSE]
  rownames(output) <- NULL
  output
}


measure_blisa_call <- function(expr, sample_interval_seconds = 0.1) {
  stopifnot(length(sample_interval_seconds) == 1L,
            is.finite(sample_interval_seconds), sample_interval_seconds > 0)
  require_packages(c("callr", "ps"))
  parent <- ps::ps_handle()
  rss_mib <- function() unname(ps::ps_memory_info(parent)["rss"]) / 1024^2
  directory <- tempfile("blisa-memory-")
  dir.create(directory)
  ready_file <- file.path(directory, "ready")
  stop_file <- file.path(directory, "stop")
  monitor <- NULL
  on.exit({
    if (!is.null(monitor) && monitor$is_alive()) monitor$kill()
    unlink(directory, recursive = TRUE)
  }, add = TRUE)

  # A separate R process can sample while BLISA is busy in native threaded code.
  # Only the calling R process is measured; the sampler's memory is excluded.
  monitor <- callr::r_bg(function(pid, interval, ready_file, stop_file) {
    parent <- ps::ps_handle(pid)
    peak <- unname(ps::ps_memory_info(parent)["rss"]) / 1024^2
    samples <- 1L
    writeLines("ready", ready_file)
    repeat {
      if (file.exists(stop_file)) break
      peak <- max(peak, unname(ps::ps_memory_info(parent)["rss"]) / 1024^2)
      samples <- samples + 1L
      Sys.sleep(interval)
    }
    list(peak_rss_mib = peak, rss_samples = samples)
  }, args = list(pid = Sys.getpid(), interval = sample_interval_seconds,
                 ready_file = ready_file, stop_file = stop_file),
  system_profile = FALSE, user_profile = FALSE)

  deadline <- proc.time()["elapsed"] + 30
  while (!file.exists(ready_file)) {
    if (!monitor$is_alive() || proc.time()["elapsed"] > deadline) {
      if (monitor$is_alive()) monitor$kill()
      stop("BLISA memory sampler failed to start: ", monitor$read_all_error())
    }
    Sys.sleep(0.01)
  }

  before_rss <- rss_mib()
  started_at <- Sys.time()
  started_cpu <- proc.time()
  value <- force(expr)
  elapsed <- proc.time() - started_cpu
  finished_at <- Sys.time()
  after_rss <- rss_mib()
  file.create(stop_file)
  monitor$wait(timeout = 5000)
  sampled <- tryCatch({
    if (monitor$is_alive()) stop("memory sampler did not stop")
    monitor$get_result()
  }, error = function(e) {
    warning("BLISA completed, but peak memory is unavailable: ", conditionMessage(e))
    list(peak_rss_mib = NA_real_, rss_samples = NA_integer_)
  })
  peak_rss <- max(sampled$peak_rss_mib, before_rss, after_rss)

  list(value = value, performance = data.frame(
    started_at_utc = format(started_at, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
    finished_at_utc = format(finished_at, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
    elapsed_seconds = unname(elapsed["elapsed"]),
    user_cpu_seconds = unname(elapsed["user.self"]),
    system_cpu_seconds = unname(elapsed["sys.self"]),
    baseline_rss_mib = before_rss,
    sampled_peak_rss_mib = peak_rss,
    peak_above_baseline_mib = peak_rss - before_rss,
    final_rss_mib = after_rss,
    sample_interval_seconds = sample_interval_seconds,
    rss_samples = sampled$rss_samples,
    timing_scope = "blisa() call only; excludes sampler startup, input preparation, exports and plotting",
    memory_scope = paste("Sampled resident memory of the calling R process during blisa();",
                         "includes already-loaded data and all threads, excludes separate processes;",
                         "peaks shorter than the sampling interval may be missed"),
    database_state = "Default database already loaded by the report's panel-coverage step"
  ))
}


project_xy <- function(xy, transform) {
  z <- cbind(as.matrix(xy), 1) %*% t(transform)
  if (any(!is.finite(z)) || any(abs(z[, 3]) < 1e-8))
    stop("Invalid homogeneous coordinates in published transform.")
  z[, 1:2, drop = FALSE] / z[, 3]
}


paperROIBorderPaths <- function(rois) {
  do.call(rbind, lapply(seq_len(nrow(rois)), function(i) {
    xy <- sf::st_coordinates(sf::st_geometry(rois)[i])
    rings <- xy[, startsWith(colnames(xy), "L"), drop = FALSE]
    data.frame(
      x = xy[, "X"], y = xy[, "Y"], region = as.character(rois$region[i]),
      # Keep disconnected polygons and interior rings as separate paths.
      piece = paste(i, apply(rings, 1, paste, collapse = "_"), sep = "_")
    )
  }))
}


annotate_visium <- function(spe, decon, workbook) {
  if (anyDuplicated(names(decon))) stop("Deconvolution table column names must be unique.")
  if (!"barcode" %in% names(decon) || anyNA(decon$barcode) ||
      any(!nzchar(decon$barcode)) || anyDuplicated(decon$barcode)) {
    stop("Deconvolution table needs a unique, non-empty barcode column.")
  }
  type_names <- setdiff(names(decon), "barcode")
  if (!length(type_names) || !all(vapply(decon[type_names], is.numeric, logical(1)))) {
    stop("Deconvolution columns must be numeric cell-type estimates.")
  }
  values <- as.matrix(decon[, type_names, drop = FALSE])
  if (any(!is.finite(values)) || any(values < 0)) stop("Invalid cell-type estimates.")
  index <- match(colnames(spe), decon$barcode)
  estimates <- values[index, , drop = FALSE]
  rownames(estimates) <- colnames(spe)
  harmonized <- harmonize_visium_cell_types(estimates)
  md <- S4Vectors::metadata(spe)
  cd <- SummarizedExperiment::colData(spe)
  existing_region <- if ("region" %in% names(cd)) as.character(cd$region) else NULL
  previous_dominant <- if ("dominant_ct" %in% names(cd)) as.character(cd$dominant_ct) else NULL
  # Reannotating a cached report must not replace its pathology assignments.
  # Only the legacy initial dominant-cell-type fallback is recomputed.
  preserve_region <- !is.null(existing_region) &&
    (!is.null(md$pathology) || !identical(existing_region, previous_dominant))
  stale <- grepl("^cell_(estimate|proportion)(_original)?::", names(cd))
  cd <- cd[, !stale, drop = FALSE]
  cd$cell_type <- cd$dominant_ct <- harmonized$dominant
  cd$cell_type_original <- cd$dominant_ct_original <- harmonized$dominant_original
  if (!preserve_region) cd$region <- harmonized$dominant
  cd$deconvolution_available <- harmonized$known
  cd$deconvolution_estimate_sum <- harmonized$total
  for (type in colnames(harmonized$estimates)) {
    cd[[paste0("cell_estimate::", type)]] <- harmonized$estimates[, type]
    cd[[paste0("cell_proportion::", type)]] <- harmonized$proportions[, type]
  }
  for (type in type_names) {
    cd[[paste0("cell_estimate_original::", type)]] <- harmonized$estimates_original[, type]
    cd[[paste0("cell_proportion_original::", type)]] <- harmonized$proportions_original[, type]
  }
  SummarizedExperiment::colData(spe) <- cd
  md$cell_type_estimates <- harmonized$estimates
  md$cell_type_proportions <- harmonized$proportions
  md$cell_type_estimates_original <- harmonized$estimates_original
  md$cell_type_proportions_original <- harmonized$proportions_original
  # Retain publication URL, source sheet, and any other existing provenance.
  deconvolution <- md$deconvolution
  if (is.null(deconvolution)) deconvolution <- list()
  updated <- list(
    source = "Data/downloaded/janesick_2023_springer/janesick_2023_source_data.xlsx", sha256 = digest::digest(file = workbook, algo = "sha256"),
    interpretation = paste("Published nonnegative estimates are summed into Xenium categories before row normalization and dominant-category assignment.",
                           "Original estimates, proportions and dominant labels are retained separately; estimates are not measured cell counts."),
    category_scheme = "Xenium manuscript cell types",
    cell_type_order = names(visium_cell_type_palette()),
    original_cell_type_order = type_names,
    category_mapping = harmonized$category_mapping,
    aggregation = "Sum original estimates within each Xenium category, then normalize by each spot's total estimate.",
    dominant_tie_rule = "First category in cell_type_order; original labels use original_cell_type_order.",
    unknown_definition = "Unmatched barcode or zero total deconvolution estimate.",
    region_definition = if (preserve_region) {
      if (!is.null(deconvolution$region_definition)) deconvolution$region_definition else
        "Preserved existing spatial region annotation"
    } else "Dominant inferred cell type, not pathology ROI",
    unmatched_spots = colnames(spe)[!harmonized$known]
  )
  deconvolution[names(updated)] <- updated
  md$deconvolution <- deconvolution
  S4Vectors::metadata(spe) <- md
  spe
}


load_visium_cellchat <- function() {
  path <- repo_path("Data", "downloaded", "cellchat_github", "CellChatDB.human.rda")
  hash <- digest::digest(file = path, algo = "sha256")
  stopifnot(identical(hash, "582e99db2bb1b1e05c0364553893d804952bd007fef870720a9cc659a7bd3ad8"))
  database <- new.env()
  load(path, envir = database)
  cache <- getFromNamespace(".blisa_cache", "blisa")
  cache[["CellChatDB.human"]] <- database$CellChatDB.human$interaction
  database <- getFromNamespace(".load_cellchat_db", "blisa")("human")
  stopifnot(nrow(database) == 3233L)
  database
}

visium_qc_dimensions <- function(spe, database) {
  detected <- Matrix::rowSums(SummarizedExperiment::assay(spe, "counts")) > 0
  coverage <- lr_panel_coverage(database, rownames(spe))
  pathways <- unique(as.character(coverage$pathway_name[coverage$panel_covered]))
  pathways <- pathways[!is.na(pathways) & nzchar(trimws(pathways))]
  c("Gene feature rows" = nrow(spe),
    "Unique gene symbols" = length(unique(rownames(spe))),
    "Detected gene symbols" = length(unique(rownames(spe)[detected])),
    "Spots" = ncol(spe), "Panel-covered LR pairs" = sum(coverage$panel_covered),
    "Panel-covered CellChatDB pathways" = length(pathways))
}

# Keep the complete spot results in one compressed CSV; large files remain local.
write_visium_spot_csv <- function(x, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file.path(directory, "lri_by_spot.csv.gz"),
    compress = "gzip", na = "NA", quote = "auto")
}
