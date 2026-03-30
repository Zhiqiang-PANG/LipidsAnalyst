#' Characterize Lipid Molecules
#'
#' Main entry function for lipidomics annotation. Identifies and characterizes
#' lipid species by comparing query spectrum against a reference lipid database.
#'
#' @param query_mz Numeric. Query m/z value from mass spectrum
#' @param query_rt Numeric. Query retention time
#' @param ms_scan_prop List. MS scan properties containing fragmentation data
#' @param ref_molecules Data frame. Reference lipid database with columns:
#'   - LipidName, Mz, LipidClass, Adduct, TotalCarbonCount,
#'     TotalDoubleBondCount, TotalOxidizedCount
#' @param ion_mode Character. Ionization mode ("positive" or "negative")
#' @param ms1_tol Numeric. MS1 mass tolerance (Da)
#' @param ms2_tol Numeric. MS2 mass tolerance (Da)
#'
#' @return List with identified lipid properties, or NULL if no match found.
#'         Best match ranked by Score field.
#'
#' @examples
#' \dontrun{
#' result <- characterize_lipid(
#'   query_mz = 760.5,
#'   query_rt = 5.2,
#'   ms_scan_prop = scan_data,
#'   ref_molecules = lipid_db,
#'   ion_mode = "positive",
#'   ms1_tol = 0.01,
#'   ms2_tol = 0.05
#' )
#' }
#'
#' @export
characterize_lipid <- function(query_mz, query_rt, ms_scan_prop,
                               ref_molecules, ion_mode,
                               ms1_tol, ms2_tol) {

  # Find starting index for m/z search
  start_id <- get_database_start_index(query_mz, ms1_tol, ref_molecules)

  candidate_molecules <- list()

  # Define acyl chain constraints for glycerolipids
  sn1_max_carbon <- 36
  sn1_max_db_bond <- 12
  sn1_min_carbon <- 2
  sn1_min_db_bond <- 0
  sn1_max_oxidized <- 0

  sn2_max_carbon <- 36
  sn2_max_db_bond <- 12
  sn2_min_carbon <- 2
  sn2_min_db_bond <- 0
  sn2_max_oxidized <- 6

  sn3_max_carbon <- 36
  sn3_max_db_bond <- 12
  sn3_min_carbon <- 2
  sn3_min_db_bond <- 0

  # Iterate through candidate molecules within m/z tolerance
  for (i in start_id:nrow(ref_molecules)) {
    molecule <- ref_molecules[i, ]
    ref_mz <- molecule$Mz
    ref_class <- molecule$LipidClass
    adduct <- molecule$Adduct

    # Check m/z window boundaries
    if (ref_mz < query_mz - ms1_tol) next
    if (ref_mz > query_mz + ms1_tol) break

    result <- NULL

    # Extract structural parameters
    total_carbon <- molecule$TotalCarbonCount
    total_db_bond <- molecule$TotalDoubleBondCount
    total_oxidized <- molecule$TotalOxidizedCount

    # Handle oxidized modification marker
    if (grepl("\\+O", molecule$LipidName)) {
      total_oxidized <- 1
    }

    # Dispatch to lipid-class-specific characterization function
    result <- switch(ref_class,
                     "PC" = judge_phosphatidylcholine(ms_scan_prop, ms2_tol, ref_mz,
                                                      total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                                      sn1_min_db_bond, sn1_max_db_bond, adduct),

                     "PE" = judge_phosphatidylethanolamine(ms_scan_prop, ms2_tol, ref_mz,
                                                           total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                                           sn1_min_db_bond, sn1_max_db_bond, adduct),

                     "PS" = judge_phosphatidylserine(ms_scan_prop, ms2_tol, ref_mz,
                                                     total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                                     sn1_min_db_bond, sn1_max_db_bond, adduct),

                     "PG" = judge_phosphatidylglycerol(ms_scan_prop, ms2_tol, ref_mz,
                                                       total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                                       sn1_min_db_bond, sn1_max_db_bond, adduct),

                     "PI" = judge_phosphatidylinositol(ms_scan_prop, ms2_tol, ref_mz,
                                                       total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                                       sn1_min_db_bond, sn1_max_db_bond, adduct),

                     "TG" = judge_triacylglycerol(ms_scan_prop, ms2_tol, ref_mz,
                                                  total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                                  sn1_min_db_bond, sn1_max_db_bond, sn2_min_carbon, sn2_max_carbon,
                                                  sn2_min_db_bond, sn2_max_db_bond, adduct),

                     "SM" = judge_sphingomyelin(ms_scan_prop, ms2_tol, ref_mz,
                                                total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                                sn1_min_db_bond, sn1_max_db_bond, adduct),

                     "LPC" = judge_lysopc(ms_scan_prop, ms2_tol, ref_mz,
                                          total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                          sn1_min_db_bond, sn1_max_db_bond, adduct),

                     "LPE" = judge_lysope(ms_scan_prop, ms2_tol, ref_mz,
                                          total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                          sn1_min_db_bond, sn1_max_db_bond, adduct),

                     "DG" = judge_dag(ms_scan_prop, ms2_tol, ref_mz,
                                      total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                      sn1_min_db_bond, sn1_max_db_bond, adduct),

                     "CE" = judge_cholesteryl_ester(ms_scan_prop, ms2_tol, ref_mz,
                                                    total_carbon, total_db_bond, adduct),

                     "CAR" = judge_acylcarnitine(ms_scan_prop, ms2_tol, ref_mz,
                                                 total_carbon, total_db_bond, adduct),

                     "FA" = judge_fatty_acid(ms_scan_prop, ms2_tol, ref_mz,
                                             total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
                                             sn1_min_db_bond, sn1_max_db_bond, adduct),

                     # Default: no match
                     NULL
    )

    # Collect valid results
    if (!is.null(result)) {
      candidate_molecules[[length(candidate_molecules) + 1]] <- result
    }
  }

  # Return best match by score, or NULL if no candidates
  if (length(candidate_molecules) > 0) {
    # Sort by score (descending)
    scores <- sapply(candidate_molecules, function(x) x$Score %||% -1)
    best_idx <- which.max(scores)
    return(candidate_molecules[[best_idx]])
  } else {
    return(NULL)
  }
}


#' Get database start index for m/z search
#'
#' Performs binary search to find starting position for m/z window filtering
#'
#' @param mz Numeric. Query m/z value
#' @param tolerance Numeric. Mass tolerance
#' @param molecules Data frame. Sorted reference molecules
#'
#' @return Integer. Starting index for m/z window search
#'
#' @keywords internal
get_database_start_index <- function(mz, tolerance, molecules) {
  lower <- mz - tolerance
  idx <- which(molecules$Mz >= lower)[1]
  return(ifelse(is.na(idx), nrow(molecules), idx))
}


