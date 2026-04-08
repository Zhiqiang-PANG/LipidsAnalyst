
# characterize_lipid <- function(query_mz, query_rt, ms_scan_prop,
#                                ref_molecules, ion_mode,
#                                ms1_tol, ms2_tol) {
#
#   # Find starting index for m/z search
#   start_id <- get_database_start_index(query_mz, ms1_tol, ref_molecules)
#
#   candidate_molecules <- list()
#
#   # Define acyl chain constraints for glycerolipids
#   sn1_max_carbon <- 36
#   sn1_max_db_bond <- 12
#   sn1_min_carbon <- 2
#   sn1_min_db_bond <- 0
#   sn1_max_oxidized <- 0
#
#   sn2_max_carbon <- 36
#   sn2_max_db_bond <- 12
#   sn2_min_carbon <- 2
#   sn2_min_db_bond <- 0
#   sn2_max_oxidized <- 6
#
#   sn3_max_carbon <- 36
#   sn3_max_db_bond <- 12
#   sn3_min_carbon <- 2
#   sn3_min_db_bond <- 0
#
#   # Iterate through candidate molecules within m/z tolerance
#   for (i in start_id:nrow(ref_molecules)) {
#     molecule <- ref_molecules[i, ]
#     ref_mz <- molecule$Mz
#     ref_class <- molecule$LipidClass
#     adduct <- molecule$Adduct
#
#     # Check m/z window boundaries
#     if (ref_mz < query_mz - ms1_tol) next
#     if (ref_mz > query_mz + ms1_tol) break
#
#     result <- NULL
#
#     # Extract structural parameters
#     total_carbon <- molecule$TotalCarbonCount
#     total_db_bond <- molecule$TotalDoubleBondCount
#     total_oxidized <- molecule$TotalOxidizedCount
#
#     # Handle oxidized modification marker
#     if (grepl("\\+O", molecule$LipidName)) {
#       total_oxidized <- 1
#     }
#
#     # Dispatch to lipid-class-specific characterization function
#     result <- switch(ref_class,
#                      "PC" = judge_phosphatidylcholine(ms_scan_prop, ms2_tol, ref_mz,
#                                                       total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                                       sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      "PE" = judge_phosphatidylethanolamine(ms_scan_prop, ms2_tol, ref_mz,
#                                                            total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                                            sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      "PS" = judge_phosphatidylserine(ms_scan_prop, ms2_tol, ref_mz,
#                                                      total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                                      sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      "PG" = judge_phosphatidylglycerol(ms_scan_prop, ms2_tol, ref_mz,
#                                                        total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                                        sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      "PI" = judge_phosphatidylinositol(ms_scan_prop, ms2_tol, ref_mz,
#                                                        total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                                        sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      "TG" = judge_triacylglycerol(ms_scan_prop, ms2_tol, ref_mz,
#                                                   total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                                   sn1_min_db_bond, sn1_max_db_bond, sn2_min_carbon, sn2_max_carbon,
#                                                   sn2_min_db_bond, sn2_max_db_bond, adduct),
#
#                      "SM" = judge_sphingomyelin(ms_scan_prop, ms2_tol, ref_mz,
#                                                 total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                                 sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      "LPC" = judge_lysopc(ms_scan_prop, ms2_tol, ref_mz,
#                                           total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                           sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      "LPE" = judge_lysope(ms_scan_prop, ms2_tol, ref_mz,
#                                           total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                           sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      "DG" = judge_dag(ms_scan_prop, ms2_tol, ref_mz,
#                                       total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                       sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      "CE" = judge_cholesteryl_ester(ms_scan_prop, ms2_tol, ref_mz,
#                                                     total_carbon, total_db_bond, adduct),
#
#                      "CAR" = judge_acylcarnitine(ms_scan_prop, ms2_tol, ref_mz,
#                                                  total_carbon, total_db_bond, adduct),
#
#                      "FA" = judge_fatty_acid(ms_scan_prop, ms2_tol, ref_mz,
#                                              total_carbon, total_db_bond, sn1_min_carbon, sn1_max_carbon,
#                                              sn1_min_db_bond, sn1_max_db_bond, adduct),
#
#                      # Default: no match
#                      NULL
#     )
#
#     # Collect valid results
#     if (!is.null(result)) {
#       candidate_molecules[[length(candidate_molecules) + 1]] <- result
#     }
#   }
#
#   # Return best match by score, or NULL if no candidates
#   if (length(candidate_molecules) > 0) {
#     # Sort by score (descending)
#     scores <- sapply(candidate_molecules, function(x) x$Score %||% -1)
#     best_idx <- which.max(scores)
#     return(candidate_molecules[[best_idx]])
#   } else {
#     return(NULL)
#   }
# }
#


camel_to_snake <- function(x) {
  gsub("_+", "_", tolower(gsub("([a-z0-9])([A-Z])", "\\1_\\2", x, perl = TRUE)))
}

resolve_judge_function <- function(name) {
  candidates <- unique(c(
    name,
    paste0(tolower(substring(name, 1, 1)), substring(name, 2)),
    camel_to_snake(name),
    camel_to_snake(sub("^JudgeIf", "judge_if_", name))
  ))

  for (nm in candidates) {
    if (exists(nm, mode = "function", inherits = TRUE)) {
      return(get(nm, mode = "function", inherits = TRUE))
    }
  }
  NULL
}

call_check_lipid <- function(fn_name, args) {
  fn <- resolve_judge_function(fn_name)
  if (is.null(fn)) return(NULL)
  do.call(fn, args)
}

resolve_csharp_arg <- function(arg, scope) {
  arg <- trimws(arg)
  if (arg == "molecule.LipidName") return(scope$molecule$LipidName)
  if (arg == "molecule.LipidClass") return(scope$molecule$LipidClass)
  if (arg == "molecule.TotalCarbonCount") return(scope$molecule$TotalCarbonCount)
  if (arg == "molecule.TotalDoubleBondCount") return(scope$molecule$TotalDoubleBondCount)
  if (exists(arg, envir = scope, inherits = FALSE)) return(get(arg, envir = scope, inherits = FALSE))
  NULL
}

invoke_case <- function(case_entry, scope) {
  args <- lapply(case_entry$args, resolve_csharp_arg, scope = scope)
  value <- call_check_lipid(case_entry$fn, args)
  list(mode = case_entry$mode, value = value)
}

read_test_spectrum <- function(input) {
  if (exists("raw_data_file_reader", mode = "function", inherits = TRUE)) {
    return(raw_data_file_reader(input))
  }
  stop("raw_data_file_reader() is not available in the current R environment.")
}

.characterize_dispatch <- list(
  "ADGGA" = list(mode = "result", fn = "JudgeIfAcylglcadg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "AHexBRS" = list(mode = "result", fn = "JudgeIfAhexbrseSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "AHexCAS" = list(mode = "result", fn = "JudgeIfAhexcaseSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "AHexCS" = list(mode = "result", fn = "JudgeIfAhexceSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "AHexCer" = list(mode = "result", fn = "JudgeIfAcylhexcer", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "sn3MinCarbon", "sn3MaxCarbon", "sn3MinDbBond", "sn3MaxDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "AHexSIS" = list(mode = "result", fn = "JudgeIfAhexsiseSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "AHexSTS" = list(mode = "result", fn = "JudgeIfAhexstseSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "ASHexCer" = list(mode = "result", fn = "JudgeIfAshexcer", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "sn3MinCarbon", "sn3MaxCarbon", "sn3MinDbBond", "sn3MaxDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "ASM" = list(mode = "result", fn = "JudgeIfAcylsm", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "BAHex" = list(mode = "result", fn = "JudgeIfSterolHexoside", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "BASulfate" = list(mode = "result", fn = "JudgeIfSterolSulfate", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "BMP" = list(mode = "result", fn = "JudgeIfBismonoacylglycerophosphate", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "BRSE" = list(mode = "result", fn = "JudgeIfBrseSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "BRSLPHex" = list(mode = "result", fn = "JudgeIfSteroidWithLpa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "BRSPHex" = list(mode = "result", fn = "JudgeIfSteroidWithPa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "BisMeLPA" = list(mode = "result", fn = "JudgeIfBismelpa", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "CAR" = list(mode = "result", fn = "JudgeIfAcylcarnitine", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "CASE" = list(mode = "result", fn = "JudgeIfCaseSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "CASLPHex" = list(mode = "result", fn = "JudgeIfSteroidWithLpa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "CASPHex" = list(mode = "result", fn = "JudgeIfSteroidWithPa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "CE" = list(mode = "result", fn = "JudgeIfCholesterylEster", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "CE_d7" = list(mode = "result", fn = "JudgeIfCholesterylEsterD7", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "CL" = list(mode = "result", fn = "JudgeIfCardiolipin", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "sn3MinCarbon", "sn3MaxCarbon", "sn3MinDbBond", "sn3MaxDbBond", "adduct")),
  "CSLPHex" = list(mode = "result", fn = "JudgeIfSteroidWithLpa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "CSPHex" = list(mode = "result", fn = "JudgeIfSteroidWithPa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "CerP" = list(mode = "result", fn = "JudgeIfCeramidePhosphate", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_ABP" = list(mode = "result", fn = "JudgeIfCeramideabp", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_ADS" = list(mode = "result", fn = "JudgeIfCeramideads", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_AH" = list(mode = "result", fn = "JudgeIfCeramideAh", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_AH_d9" = list(mode = "result", fn = "JudgeIfCeramideAhD9", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_AP" = list(mode = "result", fn = "JudgeIfCeramideap", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_AS" = list(mode = "result", fn = "JudgeIfCeramideas", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_BDS" = list(mode = "result", fn = "JudgeIfCeramidebds", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_BS" = list(mode = "result", fn = "JudgeIfCeramidebs", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_EBDS" = list(mode = "result", fn = "JudgeIfAcylcerbds", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "Cer_EODS" = list(mode = "result", fn = "JudgeIfCeramideeods", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "Cer_EOS" = list(mode = "result", fn = "JudgeIfCeramideeos", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "Cer_HDS" = list(mode = "result", fn = "JudgeIfCeramideo", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_HS" = list(mode = "result", fn = "JudgeIfCeramideo", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_NDOS" = list(mode = "result", fn = "JudgeIfCeramidedos", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_NDS" = list(mode = "result", fn = "JudgeIfCeramidends", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_NH" = list(mode = "result", fn = "JudgeIfCeramideNh", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_NH_d9" = list(mode = "result", fn = "JudgeIfCeramideNhD9", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_NP" = list(mode = "result", fn = "JudgeIfCeramidenp", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_NS" = list(mode = "result", fn = "JudgeIfCeramidens", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_NS_d7" = list(mode = "result", fn = "JudgeIfCeramidensD7", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Cer_OS" = list(mode = "result", fn = "JudgeIfCeramideos", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "CoQ" = list(mode = "result", fn = "JudgeIfCoenzymeq", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "DCAE" = list(mode = "result", fn = "JudgeIfDcae", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct", "totalOxidized")),
  "DEGSE" = list(mode = "result", fn = "JudgeIfDehydroErgoSESpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "DG" = list(mode = "result", fn = "JudgeIfDag", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DGCC" = list(mode = "result", fn = "JudgeIfDgcc", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DGDG" = list(mode = "result", fn = "JudgeIfDgdg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DGGA" = list(mode = "result", fn = "JudgeIfGlcadg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DGMG" = list(mode = "result", fn = "JudgeIfDgmg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "DGTS" = list(mode = "result", fn = "JudgeIfDgts", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DG_d5" = list(mode = "result", fn = "JudgeIfDagD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DHSph" = list(mode = "result", fn = "JudgeIfSphinganine", args = c("msScanProp", "ms2tol", "refMz", "molecule.TotalCarbonCount", "molecule.TotalDoubleBondCount", "adduct")),
  "DLCL" = list(mode = "result", fn = "JudgeIfDilysocardiolipin", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DMEDFA" = list(mode = "result", fn = "JudgeIfDmedFattyacid", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DMEDFAHFA" = list(mode = "result", fn = "JudgeIfFahfaDMED", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DMEDOxFA" = list(mode = "result", fn = "JudgeIfDmedOxfattyacid", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized")),
  "DMPE" = list(mode = "result", fn = "JudgeIfDiMethylPE", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "DSMSE" = list(mode = "result", fn = "JudgeIfDesmosterolSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "EGSE" = list(mode = "result", fn = "JudgeIfErgoSESpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "EtherDGDG" = list(mode = "result", fn = "JudgeIfEtherdgdg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherLPC" = list(mode = "result", fn = "JudgeIfEtherlysopc", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherLPE" = list(mode = "result", fn = "JudgeIfEtherlysope", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherLPG" = list(mode = "result", fn = "JudgeIfEtherlysopg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherMGDG" = list(mode = "result", fn = "JudgeIfEthermgdg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherOxPC" = list(mode = "result", fn = "JudgeIfEtheroxpc", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized", "sn1MaxOxidized", "sn2MaxOxidized")),
  "EtherOxPE" = list(mode = "result", fn = "JudgeIfEtheroxpe", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized", "sn1MaxOxidized", "sn2MaxOxidized")),
  "EtherPC" = list(mode = "result", fn = "JudgeIfEtherpc", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherPE" = list(mode = "result", fn = "JudgeIfEtherpe", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherPG" = list(mode = "result", fn = "JudgeIfEtherpg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherPI" = list(mode = "result", fn = "JudgeIfEtherpi", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherPS" = list(mode = "result", fn = "JudgeIfEtherps", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherSMGDG" = list(mode = "result", fn = "JudgeIfEtherSmgdg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "EtherTG" = list(mode = "result", fn = "JudgeIfEthertag", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "FA" = list(mode = "result", fn = "JudgeIfFattyacid", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "FAHFA" = list(mode = "result", fn = "JudgeIfFahfa", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "GD1a" = list(mode = "result", fn = "JudgeIfGD1a", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "GD1b" = list(mode = "result", fn = "JudgeIfGD1b", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "GD2" = list(mode = "result", fn = "JudgeIfGD2", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "GD3" = list(mode = "result", fn = "JudgeIfGD3", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "GDCAE" = list(mode = "result", fn = "JudgeIfGdcae", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct", "totalOxidized")),
  "GLCAE" = list(mode = "result", fn = "JudgeIfGlcae", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct", "totalOxidized")),
  "GM1" = list(mode = "result", fn = "JudgeIfGM1", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "GM3" = list(mode = "result", fn = "JudgeIfGm3", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "GPNAE" = list(mode = "result", fn = "JudgeIfGpnae", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "GQ1b" = list(mode = "result", fn = "JudgeIfGQ1b", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "GT1b" = list(mode = "result", fn = "JudgeIfGT1b", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "HBMP" = list(mode = "result", fn = "JudgeIfHemiismonoacylglycerophosphate", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "Hex2Cer" = list(mode = "result", fn = "JudgeIfHexhexceramidens", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Hex3Cer" = list(mode = "result", fn = "JudgeIfHexhexhexceramidens", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "HexCer_AP" = list(mode = "result", fn = "JudgeIfHexceramideap", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "HexCer_EOS" = list(mode = "result", fn = "JudgeIfHexceramideeos", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "HexCer_HDS" = list(mode = "result", fn = "JudgeIfHexceramideo", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "HexCer_HS" = list(mode = "result", fn = "JudgeIfHexceramideo", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "HexCer_NDS" = list(mode = "result", fn = "JudgeIfHexceramidends", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "HexCer_NS" = list(mode = "result", fn = "JudgeIfHexceramidens", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "KDCAE" = list(mode = "result", fn = "JudgeIfKdcae", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct", "totalOxidized")),
  "KLCAE" = list(mode = "result", fn = "JudgeIfKlcae", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct", "totalOxidized")),
  "LCAE" = list(mode = "result", fn = "JudgeIfLcae", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct", "totalOxidized")),
  "LDGCC" = list(mode = "result", fn = "JudgeIfLdgcc", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LDGTS" = list(mode = "result", fn = "JudgeIfLdgts", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LNAPE" = list(mode = "result", fn = "JudgeIfNacylphosphatidylethanolamine", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LNAPS" = list(mode = "result", fn = "JudgeIfNacylphosphatidylserine", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPA" = list(mode = "result", fn = "JudgeIfLysopa", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPC" = list(mode = "result", fn = "JudgeIfLysopc", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPC_d5" = list(mode = "result", fn = "JudgeIfLysopcD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPE" = list(mode = "result", fn = "JudgeIfLysope", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPE_d5" = list(mode = "result", fn = "JudgeIfLysopeD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPG" = list(mode = "result", fn = "JudgeIfLysopg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPG_d5" = list(mode = "result", fn = "JudgeIfLysopgD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPI" = list(mode = "result", fn = "JudgeIfLysopi", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPI_d5" = list(mode = "result", fn = "JudgeIfLysopiD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPS" = list(mode = "result", fn = "JudgeIfLysops", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "LPS_d5" = list(mode = "result", fn = "JudgeIfLysopsD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "MG" = list(mode = "result", fn = "JudgeIfMag", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "MGDG" = list(mode = "result", fn = "JudgeIfMgdg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "MGMG" = list(mode = "result", fn = "JudgeIfMgmg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "MIPC" = list(mode = "result", fn = "JudgeIfMipc", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "MLCL" = list(mode = "result", fn = "JudgeIfLysocardiolipin", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "MMPE" = list(mode = "result", fn = "JudgeIfMonoMethylPE", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "NA5HT" = list(mode = "result", fn = "JudgeIfNAcyl5HT", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NAAla" = list(mode = "result", fn = "JudgeIfNAcylAla", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NAAnt" = list(mode = "result", fn = "JudgeIfNAcylAnthranilicacid", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NAE" = list(mode = "result", fn = "JudgeIfAnandamide", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "NAGABA" = list(mode = "return", fn = "JudgeIfNAcylGaba", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NAGln" = list(mode = "result", fn = "JudgeIfNAcylGln", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NALeu" = list(mode = "result", fn = "JudgeIfNAcylLeu", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NAPhe" = list(mode = "result", fn = "JudgeIfNAcylPheFa", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NASer" = list(mode = "result", fn = "JudgeIfNAcylSer", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NATau" = list(mode = "result", fn = "JudgeIfNAcylTauFa", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NAVal" = list(mode = "result", fn = "JudgeIfNAcylVal", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "NGcGM3" = list(mode = "result", fn = "JudgeIfNGcGM3", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "OxFA" = list(mode = "result", fn = "JudgeIfOxfattyacid", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized")),
  "OxPC" = list(mode = "result", fn = "JudgeIfOxpc", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized", "sn1MaxOxidized", "sn2MaxOxidized")),
  "OxPE" = list(mode = "result", fn = "JudgeIfOxpe", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized", "sn1MaxOxidized", "sn2MaxOxidized")),
  "OxPG" = list(mode = "result", fn = "JudgeIfOxpg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized", "sn1MaxOxidized", "sn2MaxOxidized")),
  "OxPI" = list(mode = "result", fn = "JudgeIfOxpi", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized", "sn1MaxOxidized", "sn2MaxOxidized")),
  "OxPS" = list(mode = "result", fn = "JudgeIfOxps", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized", "sn1MaxOxidized", "sn2MaxOxidized")),
  "OxTG" = list(mode = "result", fn = "JudgeIfOxTriacylglycerol", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "totalOxidized", "adduct")),
  "PA" = list(mode = "result", fn = "JudgeIfPhosphatidicacid", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PBtOH" = list(mode = "result", fn = "JudgeIfPbtoh", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PC" = list(mode = "result", fn = "JudgeIfPhosphatidylcholine", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PC_d5" = list(mode = "result", fn = "JudgeIfPhosphatidylcholineD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PE" = list(mode = "result", fn = "JudgeIfPhosphatidylethanolamine", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PE_Cer" = list(mode = "result", fn = "JudgeIfPecermide", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized")),
  "PE_d5" = list(mode = "result", fn = "JudgeIfPhosphatidylethanolamineD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PEtOH" = list(mode = "result", fn = "JudgeIfPetoh", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PG" = list(mode = "result", fn = "JudgeIfPhosphatidylglycerol", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PG_d5" = list(mode = "result", fn = "JudgeIfPhosphatidylglycerolD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PI" = list(mode = "result", fn = "JudgeIfPhosphatidylinositol", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PI_Cer" = list(mode = "result", fn = "JudgeIfPicermide", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized")),
  "PI_d5" = list(mode = "result", fn = "JudgeIfPhosphatidylinositolD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PMeOH" = list(mode = "result", fn = "JudgeIfPmeoh", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PS" = list(mode = "result", fn = "JudgeIfPhosphatidylserine", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PS_d5" = list(mode = "result", fn = "JudgeIfPhosphatidylserineD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PT" = list(mode = "result", fn = "JudgeIfPhosphatidylThreonine", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "PhytoSph" = list(mode = "result", fn = "JudgeIfPhytosphingosine", args = c("msScanProp", "ms2tol", "refMz", "molecule.TotalCarbonCount", "molecule.TotalDoubleBondCount", "adduct")),
  "SHex" = list(mode = "result", fn = "JudgeIfSterolHexoside", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "SHexCer" = list(mode = "result", fn = "JudgeIfShexcer", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized")),
  "SISE" = list(mode = "result", fn = "JudgeIfSiseSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "SISLPHex" = list(mode = "result", fn = "JudgeIfSteroidWithLpa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "SISPHex" = list(mode = "result", fn = "JudgeIfSteroidWithPa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "SL" = list(mode = "result", fn = "JudgeIfSulfonolipid", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct", "totalOxidized")),
  "SM" = list(mode = "result", fn = "JudgeIfSphingomyelin", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "SMGDG" = list(mode = "result", fn = "JudgeIfSmgdg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "SM_d9" = list(mode = "result", fn = "JudgeIfSphingomyelinD9", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "SPE" = list(mode = "result", fn = "JudgeIfSpeSpecies", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "SPEHex" = list(mode = "result", fn = "JudgeIfSpehex", args = c("molecule.LipidName", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "SPGHex" = list(mode = "result", fn = "JudgeIfSpghex", args = c("molecule.LipidName", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "adduct")),
  "SQDG" = list(mode = "result", fn = "JudgeIfSqdg", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "SSulfate" = list(mode = "result", fn = "JudgeIfSterolSulfate", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "ST" = list(mode = "result", fn = "JudgeIfnoChainSterol", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "STSE" = list(mode = "result", fn = "JudgeIfStseSpecies", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "STSLPHex" = list(mode = "result", fn = "JudgeIfSteroidWithLpa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "STSPHex" = list(mode = "result", fn = "JudgeIfSteroidWithPa", args = c("molecule.LipidName", "molecule.LipidClass", "msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "Sph" = list(mode = "result", fn = "JudgeIfSphingosine", args = c("msScanProp", "ms2tol", "refMz", "molecule.TotalCarbonCount", "molecule.TotalDoubleBondCount", "adduct")),
  "TDCAE" = list(mode = "result", fn = "JudgeIfTdcae", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct", "totalOxidized")),
  "TG" = list(mode = "result", fn = "JudgeIfTriacylglycerol", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "TG_EST" = list(mode = "result", fn = "JudgeIfFahfaTriacylglycerol", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "sn3MinCarbon", "sn3MaxCarbon", "sn3MinDbBond", "sn3MaxDbBond", "adduct")),
  "TG_d5" = list(mode = "result", fn = "JudgeIfTriacylglycerolD5", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "sn2MinCarbon", "sn2MaxCarbon", "sn2MinDbBond", "sn2MaxDbBond", "adduct")),
  "TLCAE" = list(mode = "result", fn = "JudgeIfTlcae", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct", "totalOxidized")),
  "VAE" = list(mode = "result", fn = "JudgeIfVitaminaestermolecules", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "Vitamin_D" = list(mode = "result", fn = "JudgeIfVitaminDmolecules", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "Vitamin_E" = list(mode = "result", fn = "JudgeIfVitaminEmolecules", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "adduct")),
  "WE" = list(mode = "result", fn = "JudgeIfWaxEster", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "totalOxidized", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct")),
  "bmPC" = list(mode = "return", fn = "JudgeIfBetaMethylPhosphatidylcholine", args = c("msScanProp", "ms2tol", "refMz", "totalCarbon", "totalDbBond", "sn1MinCarbon", "sn1MaxCarbon", "sn1MinDbBond", "sn1MaxDbBond", "adduct"))
)


#' Characterize Lipid Molecules
#'
#' Main entry function for lipidomics annotation. Identifies and characterizes
#' lipid species by comparing query spectrum against a reference lipid database.
#'
#' @param query_mz Numeric. Query m/z value from mass spectrum
#' @param query_rt Numeric. Query retention time
#' @param msScanProp List. MS scan properties containing fragmentation data
#' @param RefMolecules Data frame. Reference lipid database with columns:
#'   - LipidName, Mz, LipidClass, Adduct, TotalCarbonCount,
#'     TotalDoubleBondCount, TotalOxidizedCount
#' @param ionMode Character. Ionization mode ("positive" or "negative")
#' @param ms1tol Numeric. MS1 mass tolerance (Da)
#' @param ms2tol Numeric. MS2 mass tolerance (Da)
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
characterize <- function(queryMz, queryRt, msScanProp, RefMolecules, ionMode, ms1tol, ms2tol) {
  startID <- get_database_start_index(queryMz, ms1tol, RefMolecules)
  molecules <- list()

  for (i in seq.int(startID, length(RefMolecules))) {
    molecule <- RefMolecules[[i]]
    refMz <- molecule$Mz
    lipidclass <- molecule$LipidClass
    adduct <- molecule$Adduct

    if (refMz < queryMz - ms1tol) next
    if (refMz > queryMz + ms1tol) break

    result <- NULL
    totalCarbon <- molecule$TotalCarbonCount
    totalDbBond <- molecule$TotalDoubleBondCount
    totalOxidized <- molecule$TotalOxidizedCount

    if (!is.null(molecule$LipidName) && grepl("+O", molecule$LipidName, fixed = TRUE)) {
      totalOxidized <- 1
    }

    sn1MaxCarbon <- 36
    sn1MaxDbBond <- 12
    sn1MinCarbon <- 2
    sn1MinDbBond <- 0
    sn1MaxOxidized <- 0

    sn2MaxCarbon <- 36
    sn2MaxDbBond <- 12
    sn2MinCarbon <- 2
    sn2MinDbBond <- 0
    sn2MaxOxidized <- 6

    sn3MaxCarbon <- 36
    sn3MaxDbBond <- 12
    sn3MinCarbon <- 2
    sn3MinDbBond <- 0

    scope <- list2env(list(
      msScanProp = msScanProp,
      ms2tol = ms2tol,
      refMz = refMz,
      totalCarbon = totalCarbon,
      totalDbBond = totalDbBond,
      totalOxidized = totalOxidized,
      sn1MinCarbon = sn1MinCarbon,
      sn1MaxCarbon = sn1MaxCarbon,
      sn1MinDbBond = sn1MinDbBond,
      sn1MaxDbBond = sn1MaxDbBond,
      sn2MinCarbon = sn2MinCarbon,
      sn2MaxCarbon = sn2MaxCarbon,
      sn2MinDbBond = sn2MinDbBond,
      sn2MaxDbBond = sn2MaxDbBond,
      sn3MinCarbon = sn3MinCarbon,
      sn3MaxCarbon = sn3MaxCarbon,
      sn3MinDbBond = sn3MinDbBond,
      sn3MaxDbBond = sn3MaxDbBond,
      sn1MaxOxidized = sn1MaxOxidized,
      sn2MaxOxidized = sn2MaxOxidized,
      adduct = adduct,
      molecule = molecule
    ), parent = emptyenv())

    if (identical(lipidclass, "NAGly")) {
      if (totalCarbon == sn1MinCarbon && totalCarbon == sn1MaxCarbon) {
        result <- call_check_lipid("JudgeIfNAcylGlyOxFa", list(msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, totalOxidized, adduct))
      } else {
        result <- call_check_lipid("JudgeIfFahfamidegly", list(msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, sn1MinCarbon, sn1MaxCarbon, sn1MinDbBond, sn1MaxDbBond, adduct))
      }
    } else if (identical(lipidclass, "NAGlySer")) {
      if (totalCarbon == sn1MinCarbon && totalCarbon == sn1MaxCarbon) {
        result <- call_check_lipid("JudgeIfNAcylGlySerOxFa", list(msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, totalOxidized, adduct))
      } else {
        result <- call_check_lipid("JudgeIfFahfamideglyser", list(msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, sn1MinCarbon, sn1MaxCarbon, sn1MinDbBond, sn1MaxDbBond, adduct))
      }
    } else if (identical(lipidclass, "NAOrn")) {
      if (totalCarbon == sn1MinCarbon && totalCarbon == sn1MaxCarbon) {
        result <- call_check_lipid("JudgeIfNAcylOrnOxFa", list(msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, totalOxidized, adduct))
      } else {
        result <- call_check_lipid("JudgeIfFahfamideorn", list(msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, sn1MinCarbon, sn1MaxCarbon, sn1MinDbBond, sn1MaxDbBond, adduct))
      }
    } else if (identical(lipidclass, "NATryA")) {
      if (totalCarbon < 29) {
        return(call_check_lipid("JudgeIfNAcylTryA", list(msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, totalOxidized, sn1MinCarbon, sn1MaxCarbon, sn1MinDbBond, sn1MaxDbBond, adduct)))
      } else {
        return(call_check_lipid("JudgeIfFahfamideTrya", list(msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, sn1MinCarbon, sn1MaxCarbon, sn1MinDbBond, sn1MaxDbBond, adduct)))
      }
    } else {
      case_entry <- .characterize_dispatch[[lipidclass]]
      if (is.null(case_entry)) return(NULL)
      out <- invoke_case(case_entry, scope)
      if (identical(out$mode, "return")) return(out$value)
      result <- out$value
    }

    if (!is.null(result)) {
      molecules[[length(molecules) + 1]] <- result
    }
  }

  if (length(molecules) > 0) {
    scores <- vapply(molecules, function(m) m$Score %||% -Inf, numeric(1))
    return(molecules[[which.max(scores)]])
  }

  NULL
}

convert_to_required_spectrum_format <- function(peaks) {
  if (is.null(peaks) || length(peaks) == 0) return(list())

  intensities <- vapply(peaks, function(p) p$Intensity %||% 0, numeric(1))
  maxintensity <- max(intensities)
  if (maxintensity == 0) maxintensity <- 1

  spectrum <- lapply(peaks, function(peak) {
    list(Mass = peak$Mass, Intensity = (peak$Intensity / maxintensity) * 100.0)
  })

  spectrum[order(vapply(spectrum, function(x) x$Intensity, numeric(1)), decreasing = TRUE)]
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
  targetMass <- mz - tolerance
  startIndex <- 1
  endIndex <- length(molecules)

  if (endIndex < 1) return(1)
  if (targetMass > molecules[[endIndex]]$Mz) return(endIndex)

  counter <- 0
  while (counter < 10) {
    mid <- floor((startIndex + endIndex) / 2)
    if (molecules[[startIndex]]$Mz <= targetMass && targetMass < molecules[[mid]]$Mz) {
      endIndex <- mid
    } else if (molecules[[mid]]$Mz <= targetMass && targetMass < molecules[[endIndex]]$Mz) {
      startIndex <- mid
    }
    counter <- counter + 1
  }
  startIndex
}

# [0] Name [1] m/z [2] adduct
read_library <- function(file) {
  if (!file.exists(file)) stop("Library file not found: ", file)

  lines <- readLines(file, warn = FALSE, encoding = "ASCII")
  if (length(lines) <= 1) return(list())

  molecules <- list()
  for (line in lines[-1]) {
    lineArray <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(lineArray) < 3) next

    nameString <- lineArray[[1]]
    if (!grepl(":", nameString, fixed = TRUE)) {
      nameString <- paste0(nameString, " :")
    }

    parts <- strsplit(nameString, " ", fixed = TRUE)[[1]]
    if (length(parts) < 2) next

    lipidClass <- parts[[1]]
    lipidClassEnum <- get_lipid_class_enum(lipidClass)
    if (identical(lipidClassEnum, "Undefined")) next

    mzValue <- suppressWarnings(as.numeric(lineArray[[2]]))
    if (is.na(mzValue)) next

    adductString <- lineArray[[3]]
    adduct <- if (exists("get_adduct_ion", mode = "function", inherits = TRUE)) {
      get_adduct_ion(adductString)
    } else {
      list(AdductIonName = adductString, FormatCheck = TRUE)
    }
    if (is.null(adduct$FormatCheck) || !isTRUE(adduct$FormatCheck)) next

    chainString <- parts[[2]]
    totalCarbonString <- strsplit(chainString, ":", fixed = TRUE)[[1]][1]
    totalCarbonString <- gsub("[dtm]", "", totalCarbonString)

    bondString <- strsplit(nameString, ":", fixed = TRUE)[[1]][2]
    bondParts <- strsplit(bondString, "+", fixed = TRUE)[[1]]
    totalDoubleBondString <- gsub("e", "", bondParts[[1]])
    totalOxidizedString <- if (length(bondParts) > 1) gsub("O", "", bondParts[[2]]) else "0"

    totalCarbon <- suppressWarnings(as.integer(totalCarbonString))
    totalDoubleBond <- suppressWarnings(as.integer(totalDoubleBondString))
    totalOxidized <- suppressWarnings(as.integer(totalOxidizedString))

    molecules[[length(molecules) + 1]] <- list(
      LipidName = lineArray[[1]],
      SublevelLipidName = lineArray[[1]],
      LipidClass = lipidClassEnum,
      Adduct = adduct,
      Mz = mzValue,
      TotalChainString = chainString,
      TotalCarbonCount = totalCarbon,
      TotalDoubleBondCount = totalDoubleBond,
      TotalOxidizedCount = totalOxidized
    )
  }

  if (length(molecules) == 0) return(list())
  molecules[order(vapply(molecules, function(m) m$Mz, numeric(1)))]
}

get_lipid_class_enum <- function(lipidClass) {
  known <- unique(c(names(.characterize_dispatch), "NAGly", "NAGlySer", "NAOrn", "NATryA"))
  if (lipidClass %in% known) lipidClass else "Undefined"
}



{

  `%||%` <- function(x, y) if (is.null(x)) y else x

  camel_to_snake <- function(x) {
    gsub("_+", "_", tolower(gsub("([a-z0-9])([A-Z])", "\\1_\\2", x, perl = TRUE)))
  }

  resolve_judge_function <- function(name) {
    candidates <- unique(c(
      name,
      paste0(tolower(substring(name, 1, 1)), substring(name, 2)),
      camel_to_snake(name),
      camel_to_snake(sub("^JudgeIf", "judge_if_", name))
    ))

    for (nm in candidates) {
      if (exists(nm, mode = "function", inherits = TRUE)) {
        return(get(nm, mode = "function", inherits = TRUE))
      }
    }
    NULL
  }

  call_judge <- function(fn_name, args) {
    fn <- resolve_judge_function(fn_name)
    if (is.null(fn)) return(NULL)
    do.call(fn, args)
  }

  resolve_csharp_arg <- function(arg, scope) {
    arg <- trimws(arg)
    if (arg == "molecule.LipidName") return(scope$molecule$LipidName)
    if (arg == "molecule.LipidClass") return(scope$molecule$LipidClass)
    if (arg == "molecule.TotalCarbonCount") return(scope$molecule$TotalCarbonCount)
    if (arg == "molecule.TotalDoubleBondCount") return(scope$molecule$TotalDoubleBondCount)
    if (exists(arg, envir = scope, inherits = FALSE)) return(get(arg, envir = scope, inherits = FALSE))
    NULL
  }

  invoke_case <- function(case_entry, scope) {
    args <- lapply(case_entry$args, resolve_csharp_arg, scope = scope)
    call_judge(case_entry$fn, args)
  }

  get_lipid_annotation_dispatch <- function() {
    # if (exists(".characterize_dispatch", inherits = TRUE)) {
    #   return(get(".characterize_dispatch", inherits = TRUE))
    # }
    #
    # script_path <- "R/lipidomics.R"
    # if (!file.exists(script_path)) {
    #   stop("Cannot find switch-case dispatch table: ", script_path)
    # }
    #
    # env <- new.env(parent = globalenv())
    # sys.source(script_path, envir = env)
    # if (!exists(".characterize_dispatch", envir = env, inherits = FALSE)) {
    #   stop(".characterize_dispatch was not found after sourcing R/lipidomics.R")
    # }

    get(".characterize_dispatch", envir = env, inherits = FALSE)
  }

  GetLipidMoleculeAnnotationResult <- function(msScanProp, molecule, ms2tol) {
    lipidclass <- molecule$LipidClass
    refMz <- molecule$Mz
    adduct <- molecule$Adduct

    totalCarbon <- molecule$TotalCarbonCount
    totalDbBond <- molecule$TotalDoubleBondCount
    totalOxidized <- molecule$TotalOxidizedCount

    sn1Carbon <- molecule$Sn1CarbonCount %||% 0
    sn1DbBond <- molecule$Sn1DoubleBondCount %||% 0
    sn1Oxidized <- molecule$Sn1Oxidizedount %||% 0
    sn2Oxidized <- molecule$Sn2Oxidizedount %||% 0

    sn2Carbon <- molecule$Sn2CarbonCount %||% 0
    sn2DbBond <- molecule$Sn2DoubleBondCount %||% 0
    sn3Carbon <- molecule$Sn3CarbonCount %||% 0
    sn3DbBond <- molecule$Sn3DoubleBondCount %||% 0

    if (exists("get_lipid_header_string", mode = "function", inherits = TRUE)) {
      lipidheader <- get_lipid_header_string(molecule$LipidName)
    } else {
      lipidheader <- NULL
    }

    dispatch <- get_lipid_annotation_dispatch()

    scope <- list2env(list(
      msScanProp = msScanProp,
      ms2tol = ms2tol,
      refMz = refMz,
      totalCarbon = totalCarbon,
      totalDbBond = totalDbBond,
      totalOxidized = totalOxidized,
      sn1MinCarbon = sn1Carbon,
      sn1MaxCarbon = sn1Carbon,
      sn1MinDbBond = sn1DbBond,
      sn1MaxDbBond = sn1DbBond,
      sn2MinCarbon = sn2Carbon,
      sn2MaxCarbon = sn2Carbon,
      sn2MinDbBond = sn2DbBond,
      sn2MaxDbBond = sn2DbBond,
      sn3MinCarbon = sn3Carbon,
      sn3MaxCarbon = sn3Carbon,
      sn3MinDbBond = sn3DbBond,
      sn3MaxDbBond = sn3DbBond,
      sn1MaxOxidized = sn1Oxidized,
      sn2MaxOxidized = sn2Oxidized,
      adduct = adduct,
      molecule = molecule,
      lipidheader = lipidheader
    ), parent = emptyenv())

    if (identical(lipidclass, "SM")) {
      total_chain <- molecule$TotalChainString %||% ""
      if (grepl("O3", total_chain, fixed = TRUE)) {
        return(call_judge("JudgeIfSphingomyelinPhyto", list(
          msScanProp, ms2tol, refMz,
          totalCarbon, totalDbBond, sn1Carbon, sn1Carbon, sn1DbBond, sn1DbBond, adduct
        )))
      }
      return(call_judge("JudgeIfSphingomyelin", list(
        msScanProp, ms2tol, refMz,
        totalCarbon, totalDbBond, sn1Carbon, sn1Carbon, sn1DbBond, sn1DbBond, adduct
      )))
    }

    if (identical(lipidclass, "CL") && sn3Carbon < 1) {
      return(call_judge("JudgeIfCardiolipin", list(
        msScanProp, ms2tol, refMz,
        totalCarbon, totalDbBond, sn1Carbon, sn1Carbon, sn1DbBond, sn1DbBond, adduct
      )))
    }

    if (identical(lipidclass, "NAGly")) {
      if (totalCarbon == sn1Carbon) {
        return(call_judge("JudgeIfNAcylGlyOxFa", list(
          msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, totalOxidized, adduct
        )))
      }
      return(call_judge("JudgeIfFahfamidegly", list(
        msScanProp, ms2tol, refMz,
        totalCarbon, totalDbBond, sn1Carbon, sn1Carbon, sn1DbBond, sn1DbBond, adduct
      )))
    }

    if (identical(lipidclass, "NAGlySer")) {
      if (totalCarbon == sn1Carbon) {
        return(call_judge("JudgeIfNAcylGlySerOxFa", list(
          msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, totalOxidized, adduct
        )))
      }
      return(call_judge("JudgeIfFahfamideglyser", list(
        msScanProp, ms2tol, refMz,
        totalCarbon, totalDbBond, sn1Carbon, sn1Carbon, sn1DbBond, sn1DbBond, adduct
      )))
    }

    if (identical(lipidclass, "NAOrn")) {
      if (totalCarbon == sn1Carbon) {
        return(call_judge("JudgeIfNAcylOrnOxFa", list(
          msScanProp, ms2tol, refMz, totalCarbon, totalDbBond, totalOxidized, adduct
        )))
      }
      return(call_judge("JudgeIfFahfamideorn", list(
        msScanProp, ms2tol, refMz,
        totalCarbon, totalDbBond, sn1Carbon, sn1Carbon, sn1DbBond, sn1DbBond, adduct
      )))
    }

    if (identical(lipidclass, "NATryA")) {
      if (totalCarbon < 29) {
        return(call_judge("JudgeIfNAcylTryA", list(
          msScanProp, ms2tol, refMz,
          totalCarbon, totalDbBond, totalOxidized, sn1Carbon, sn1Carbon, sn1DbBond, sn1DbBond, adduct
        )))
      }
      return(call_judge("JudgeIfFahfamideTrya", list(
        msScanProp, ms2tol, refMz,
        totalCarbon, totalDbBond, sn1Carbon, sn1Carbon, sn1DbBond, sn1DbBond, adduct
      )))
    }

    case_entry <- dispatch[[lipidclass]]
    if (is.null(case_entry)) {
      return(NULL)
    }

    invoke_case(case_entry, scope)
  }



}

{

  get_field <- function(obj, name, default = NULL) {
    if (is.null(obj)) return(default)
    if (is.list(obj) && !is.null(obj[[name]])) return(obj[[name]])
    if (is.environment(obj) && exists(name, envir = obj, inherits = FALSE)) return(get(name, envir = obj, inherits = FALSE))
    default
  }

  resolve_function <- function(name) {
    if (exists(name, mode = "function", inherits = TRUE)) {
      return(get(name, mode = "function", inherits = TRUE))
    }
    NULL
  }

  is_compared_available_scan_ref <- function(msScanProp, molMsRef) {
    s1 <- get_field(msScanProp, "Spectrum")
    s2 <- get_field(molMsRef, "Spectrum")
    !is.null(s1) && !is.null(s2) && length(s1) > 0 && length(s2) > 0
  }

  convert_reference_to_lipid <- function(molMsRef) {
    fn <- resolve_function("convert_msdial_lipidname_to_lipid_molecule_object_vs2")
    if (!is.null(fn)) return(fn(molMsRef))

    fn <- resolve_function("ConvertMsdialLipidnameToLipidMoleculeObjectVS2")
    if (!is.null(fn)) return(fn(molMsRef))

    NULL
  }

  ensure_get_lipid_molecule_annotation_result <- function() {
    if (exists("GetLipidMoleculeAnnotationResult", mode = "function", inherits = TRUE)) {
      return(get("GetLipidMoleculeAnnotationResult", mode = "function", inherits = TRUE))
    }

    script_path <- "GetLipidMoleculeAnnotationResult.R"
    if (!file.exists(script_path)) {
      stop("GetLipidMoleculeAnnotationResult.R is required but not found.")
    }

    env <- new.env(parent = globalenv())
    sys.source(script_path, envir = env)
    if (!exists("GetLipidMoleculeAnnotationResult", envir = env, inherits = FALSE)) {
      stop("GetLipidMoleculeAnnotationResult() not found in GetLipidMoleculeAnnotationResult.R")
    }

    get("GetLipidMoleculeAnnotationResult", envir = env, inherits = FALSE)
  }

  ensure_get_matched_peaks_scores <- function() {
    fn <- resolve_function("GetMatchedPeaksScores")
    if (!is.null(fn)) return(fn)

    fn <- resolve_function("get_matched_peaks_scores")
    if (!is.null(fn)) return(fn)

    stop("GetMatchedPeaksScores()/get_matched_peaks_scores() is required in the R environment.")
  }

  # Translated from C# MsScanMatching.GetLipidomicsMatchedPeaksScores
  GetLipidomicsMatchedPeaksScores <- function(msScanProp, molMsRef, bin, massBegin, massEnd) {
    if (!is_compared_available_scan_ref(msScanProp, molMsRef)) {
      return(c(-1, -1))
    }

    get_matched <- ensure_get_matched_peaks_scores()
    resultArray <- get_matched(msScanProp, molMsRef, bin, massBegin, massEnd)

    compClass <- get_field(molMsRef, "CompoundClass", "")
    comment <- get_field(molMsRef, "Comment", "")

    if (comment != "SPLASH" && compClass != "Unknown" && compClass != "Others") {
      molecule <- convert_reference_to_lipid(molMsRef)
      if (is.null(molecule) || is.null(get_field(molecule, "Adduct"))) {
        return(resultArray)
      }

      isEtherPE <- identical(get_field(molecule, "LipidClass"), "EtherPE")
      refSpectrum <- get_field(molMsRef, "Spectrum", list())
      ionMode <- get_field(msScanProp, "IonMode")
      if (isTRUE(isEtherPE) && length(refSpectrum) == 3 && identical(ionMode, "Positive")) {
        return(resultArray)
      }

      get_result <- ensure_get_lipid_molecule_annotation_result()
      result <- get_result(msScanProp, molecule, bin)

      if (!is.null(result)) {
        annotationLevel <- get_field(result, "AnnotationLevel", NA)

        if (isTRUE(annotationLevel == 1)) {
          lipidName <- get_field(molecule, "LipidName", "")
          if (identical(compClass, "SM") && (grepl("3O", lipidName, fixed = TRUE) || grepl("O3", lipidName, fixed = TRUE))) {
            resultArray[1] <- 2.0
            return(resultArray)
          }

          resultArray[1] <- 1.0
          return(resultArray)
        }

        if (isTRUE(annotationLevel >= 2)) {
          resultArray[1] <- 2.0
          return(resultArray)
        }

        return(resultArray)
      }

      return(resultArray)
    }

    resultArray
  }

  source_type_generated_lipid_flag <- 32

  is_source_generated_lipid <- function(source) {
    if (is.null(source)) return(FALSE)
    if (is.character(source)) return(identical(source, "GeneratedLipid"))
    if (is.numeric(source)) return(bitwAnd(as.integer(source), source_type_generated_lipid_flag) != 0L)
    FALSE
  }

  is_collision_in <- function(collisionType, values) {
    is.character(collisionType) && length(collisionType) > 0 && collisionType %in% values
  }

  read_chrom_rt_value <- function(obj) {
    chrom <- get_field(obj, "ChromXs")
    if (is.null(chrom)) return(NULL)
    rt <- get_field(chrom, "RT")
    if (is.null(rt)) return(NULL)
    get_field(rt, "Value")
  }

  compute_ms1_tolerance <- function(parameter, precursorMz) {
    fn <- resolve_function("FixMassTolerance")
    if (!is.null(fn)) {
      return(fn(parameter$Ms1Tolerance, precursorMz))
    }

    fn <- resolve_function("fix_mass_tolerance")
    if (!is.null(fn)) {
      return(fn(parameter$Ms1Tolerance, precursorMz))
    }

    parameter$Ms1Tolerance
  }

  ensure_validate_helpers <- function() {
    if (!exists("ValidateOnLipidomics", mode = "function", inherits = TRUE) && file.exists("ValidateOnLipidomics.R")) {
      sys.source("ValidateOnLipidomics.R", envir = .GlobalEnv)
    }

    if (!exists("ValidateOnLipidomics", mode = "function", inherits = TRUE)) {
      stop("ValidateOnLipidomics() is required but not found.")
    }
  }

  ensure_scoring_fn <- function(name_candidates) {
    for (nm in name_candidates) {
      fn <- resolve_function(nm)
      if (!is.null(fn)) return(fn)
    }
    NULL
  }

  safe_total_average <- function(values) {
    if (length(values) == 0) return(0)
    mean(values)
  }

  get_dot_product <- function(result, squared_name, direct_name) {
    direct <- get_field(result, direct_name)
    if (!is.null(direct)) return(direct)
    squared <- get_field(result, squared_name, -1)
    if (is.null(squared) || squared < 0) return(-1)
    sqrt(max(as.numeric(squared), 0))
  }

  # Translated from MsReferenceScorer.Score
  Score <- function(query, reference, scorerContext) {
    CalculateScore(
      property = query$Property,
      scan = query$NormalizedScan,
      scanIsotopes = query$Isotopes,
      reference = reference,
      referenceIsotopes = reference$IsotopicPeaks,
      parameter = query$Parameter,
      scorerContext = scorerContext
    )
  }

  # Translated from MsReferenceScorer.CalculateScore
  CalculateScore <- function(property, scan, scanIsotopes, reference, referenceIsotopes, parameter, scorerContext) {
    get_weighted <- ensure_scoring_fn(c("GetWeightedDotProduct", "get_weighted_dot_product"))
    get_simple <- ensure_scoring_fn(c("GetSimpleDotProduct", "get_simple_dot_product"))
    get_reverse <- ensure_scoring_fn(c("GetReverseDotProduct", "get_reverse_dot_product"))
    get_gaussian <- ensure_scoring_fn(c("GetGaussianSimilarity", "get_gaussian_similarity"))
    get_isotope <- ensure_scoring_fn(c("GetIsotopeRatioSimilarity", "get_isotope_ratio_similarity"))

    if (is.null(get_weighted) || is.null(get_simple) || is.null(get_reverse) || is.null(get_gaussian) || is.null(get_isotope)) {
      stop("Required scoring functions are missing in the R environment.")
    }

    sqweightedDotProduct <- get_weighted(scan, reference, parameter$Ms2Tolerance, parameter$MassRangeBegin, parameter$MassRangeEnd)
    sqsimpleDotProduct <- get_simple(scan, reference, parameter$Ms2Tolerance, parameter$MassRangeBegin, parameter$MassRangeEnd)
    sqreverseDotProduct <- get_reverse(scan, reference, parameter$Ms2Tolerance, parameter$MassRangeBegin, parameter$MassRangeEnd)

    refSpectrum <- get_field(reference, "Spectrum")
    spectrumPenalty <- !is.null(refSpectrum) && length(refSpectrum) == 1

    matchedPeaksScores <- NULL
    omics <- scorerContext$omics
    source <- scorerContext$source
    collisionType <- scorerContext$collisionType

    if (identical(omics, "Lipidomics")) {
      if (is_source_generated_lipid(source)) {
        if (identical(collisionType, "EIEIO")) {
          fn <- ensure_scoring_fn(c("GetEieioBasedLipidomicsMatchedPeaksScores", "get_eieio_based_lipidomics_matched_peaks_scores"))
        } else if (identical(collisionType, "EID")) {
          fn <- ensure_scoring_fn(c("GetEidBasedLipidomicsMatchedPeaksScores", "get_eid_based_lipidomics_matched_peaks_scores"))
        } else if (identical(collisionType, "OAD")) {
          fn <- ensure_scoring_fn(c("GetOadBasedLipidomicsMatchedPeaksScores", "get_oad_based_lipidomics_matched_peaks_scores"))
        } else {
          fn <- ensure_scoring_fn(c("GetLipidomicsMatchedPeaksScores", "get_lipidomics_matched_peaks_scores"))
        }
      } else {
        if (identical(collisionType, "EIEIO")) {
          fn <- ensure_scoring_fn(c("GetLipidomicsMoleculerSpeciesLevelAnnotationPeaksScoresForEIEIO", "get_lipidomics_moleculer_species_level_annotation_peaks_scores_for_eieio"))
        } else if (identical(collisionType, "EID")) {
          fn <- ensure_scoring_fn(c("GetLipidomicsMoleculerSpeciesLevelAnnotationPeaksScoresForEIEIO", "get_lipidomics_moleculer_species_level_annotation_peaks_scores_for_eieio"))
        } else if (identical(collisionType, "OAD")) {
          fn <- ensure_scoring_fn(c("GetLipidomicsMoleculerSpeciesLevelAnnotationPeaksScoresForOAD", "get_lipidomics_moleculer_species_level_annotation_peaks_scores_for_oad"))
        } else {
          fn <- ensure_scoring_fn(c("GetLipidomicsMatchedPeaksScores", "get_lipidomics_matched_peaks_scores"))
        }
      }
    } else {
      fn <- ensure_scoring_fn(c("GetMatchedPeaksScores", "get_matched_peaks_scores"))
    }

    if (is.null(fn)) {
      stop("Matched-peaks scoring function is missing for the current context.")
    }
    matchedPeaksScores <- fn(scan, reference, parameter$Ms2Tolerance, parameter$MassRangeBegin, parameter$MassRangeEnd)

    ms1Tol <- compute_ms1_tolerance(parameter, property$PrecursorMz)
    ms1Similarity <- get_gaussian(property$PrecursorMz, reference$PrecursorMz, ms1Tol)
    isotopeSimilarity <- get_isotope(scanIsotopes, referenceIsotopes, property$PrecursorMz, ms1Tol)

    result <- list(
      Name = reference$Name,
      LibraryID = reference$ScanID,
      InChIKey = reference$InChIKey,
      SquaredWeightedDotProduct = as.numeric(sqweightedDotProduct),
      SquaredSimpleDotProduct = as.numeric(sqsimpleDotProduct),
      SquaredReverseDotProduct = as.numeric(sqreverseDotProduct),
      MatchedPeaksPercentage = as.numeric(matchedPeaksScores[1]),
      MatchedPeaksCount = as.numeric(matchedPeaksScores[2]),
      AcurateMassSimilarity = as.numeric(ms1Similarity),
      IsotopeSimilarity = as.numeric(isotopeSimilarity),
      Source = source,
      AnnotatorID = scorerContext$id,
      Priority = scorerContext$priority,
      IsSpectrumMatch = FALSE,
      IsPrecursorMzMatch = FALSE,
      IsRtMatch = FALSE,
      IsCcsMatch = FALSE,
      IsLipidClassMatch = FALSE,
      IsLipidChainsMatch = FALSE,
      IsLipidPositionMatch = FALSE,
      IsLipidDoubleBondPositionMatch = FALSE,
      IsOtherLipidMatch = FALSE,
      IsReferenceMatched = FALSE,
      IsAnnotationSuggested = FALSE
    )

    if (isTRUE(parameter$IsUseTimeForAnnotationScoring)) {
      rtProperty <- read_chrom_rt_value(property)
      rtReference <- read_chrom_rt_value(reference)
      if (!is.null(rtProperty) && !is.null(rtReference)) {
        result$RtSimilarity <- as.numeric(get_gaussian(rtProperty, rtReference, parameter$RtTolerance))
      }
    }
    if (isTRUE(parameter$IsUseCcsForAnnotationScoring)) {
      result$CcsSimilarity <- as.numeric(get_gaussian(property$CollisionCrossSection, reference$CollisionCrossSection, parameter$CcsTolerance))
    }

    dotProductFactor <- 3.0
    revesrseDotProdFactor <- 2.0
    presensePercentageFactor <- 1.0
    msmsFactor <- 3.0
    rtFactor <- 1.0
    ccsFactor <- 1.0
    massFactor <- 1.0
    isotopeFactor <- 0.0

    if (identical(omics, "Lipidomics")) {
      if (is_collision_in(collisionType, c("CID", "HCD"))) {
        dotProductFactor <- 1.0
        revesrseDotProdFactor <- 2.0
        presensePercentageFactor <- 3.0
        rtFactor <- 0.5
        ccsFactor <- 0.5
      } else {
        dotProductFactor <- 1.0
        revesrseDotProdFactor <- 2.0
        presensePercentageFactor <- 2.0
        rtFactor <- 0.5
        ccsFactor <- 0.5
      }
    }
    if (identical(omics, "Metabolomics") && isTRUE(spectrumPenalty)) {
      dotProductFactor <- 1.5
      revesrseDotProdFactor <- 1.0
      presensePercentageFactor <- 0.5
    }

    result$WeightedDotProduct <- get_dot_product(result, "SquaredWeightedDotProduct", "WeightedDotProduct")
    result$SimpleDotProduct <- get_dot_product(result, "SquaredSimpleDotProduct", "SimpleDotProduct")
    result$ReverseDotProduct <- get_dot_product(result, "SquaredReverseDotProduct", "ReverseDotProduct")

    msmsScore <- (
      result$WeightedDotProduct * dotProductFactor +
        result$SimpleDotProduct * dotProductFactor +
        result$ReverseDotProduct * revesrseDotProdFactor +
        result$MatchedPeaksPercentage * presensePercentageFactor
    ) / (dotProductFactor + dotProductFactor + revesrseDotProdFactor + presensePercentageFactor)

    scores <- c()
    if (!is.null(result$AcurateMassSimilarity) && result$AcurateMassSimilarity >= 0 && massFactor > 0)
      scores <- c(scores, result$AcurateMassSimilarity * massFactor)
    if (result$WeightedDotProduct >= 0 && result$SimpleDotProduct >= 0 && result$ReverseDotProduct >= 0)
      scores <- c(scores, msmsScore * msmsFactor)
    if (isTRUE(parameter$IsUseTimeForAnnotationScoring) && !is.null(result$RtSimilarity) && result$RtSimilarity >= 0 && rtFactor > 0)
      scores <- c(scores, result$RtSimilarity * rtFactor)
    if (isTRUE(parameter$IsUseCcsForAnnotationScoring) && !is.null(result$CcsSimilarity) && result$CcsSimilarity >= 0 && ccsFactor > 0)
      scores <- c(scores, result$CcsSimilarity * ccsFactor)
    if (!is.null(result$IsotopeSimilarity) && result$IsotopeSimilarity >= 0 && isotopeFactor > 0)
      scores <- c(scores, result$IsotopeSimilarity * isotopeFactor)

    result$TotalScore <- as.numeric(safe_total_average(scores))
    if (is.null(result$InChIKey) || (is.character(result$InChIKey) && !nzchar(result$InChIKey))) {
      result$TotalScore <- result$TotalScore * 0.9
    }

    result <- Validate(result, property, scan, reference, parameter, scorerContext)
    result
  }

  validate_base <- function(result, property, reference, parameter, omics) {
    if (identical(omics, "Lipidomics")) {
      result$IsSpectrumMatch <- result$SquaredWeightedDotProduct >= parameter$SquaredWeightedDotProductCutOff ||
        result$SquaredSimpleDotProduct >= parameter$SquaredSimpleDotProductCutOff ||
        result$SquaredReverseDotProduct >= parameter$SquaredReverseDotProductCutOff

      compoundClass <- get_field(reference, "CompoundClass", "")
      if ((identical(compoundClass, "EtherTG") || identical(compoundClass, "EtherDG")) &&
          result$SquaredSimpleDotProduct < parameter$SquaredSimpleDotProductCutOff) {
        result$IsSpectrumMatch <- FALSE
      }
    } else {
      result$IsSpectrumMatch <- result$SquaredWeightedDotProduct >= parameter$SquaredWeightedDotProductCutOff &&
        result$SquaredSimpleDotProduct >= parameter$SquaredSimpleDotProductCutOff &&
        result$SquaredReverseDotProduct >= parameter$SquaredReverseDotProductCutOff &&
        result$MatchedPeaksPercentage >= parameter$MatchedPeaksPercentageCutOff &&
        result$MatchedPeaksCount >= parameter$MinimumSpectrumMatch
    }

    ms1Tol <- compute_ms1_tolerance(parameter, property$PrecursorMz)
    result$IsPrecursorMzMatch <- abs(property$PrecursorMz - reference$PrecursorMz) <= ms1Tol

    if (isTRUE(parameter$IsUseTimeForAnnotationScoring)) {
      rtProperty <- read_chrom_rt_value(property)
      rtReference <- read_chrom_rt_value(reference)
      if (!is.null(rtProperty) && !is.null(rtReference)) {
        result$IsRtMatch <- abs(rtProperty - rtReference) <= parameter$RtTolerance
      }
    }

    if (isTRUE(parameter$IsUseCcsForAnnotationScoring)) {
      result$IsCcsMatch <- abs(property$CollisionCrossSection - reference$CollisionCrossSection) <= parameter$CcsTolerance
    }

    result
  }

  # Translated from MsReferenceScorer.Validate
  Validate <- function(result, property, scan, reference, parameter, scorerContext) {
    omics <- scorerContext$omics
    source <- scorerContext$source
    collisionType <- scorerContext$collisionType
    useMs2 <- isTRUE(scorerContext$useMs2)

    result <- validate_base(result, property, reference, parameter, omics)

    if (identical(omics, "Lipidomics")) {
      if (is_source_generated_lipid(source) && is_collision_in(collisionType, c("EIEIO", "OAD", "EID"))) {
        if (exists("ValidateOnEadLipidomics", mode = "function", inherits = TRUE)) {
          result <- ValidateOnEadLipidomics(result, scan, reference, parameter, collisionType)
        } else {
          stop("ValidateOnEadLipidomics() is required for GeneratedLipid with EIEIO/OAD/EID.")
        }
      } else {
        ensure_validate_helpers()
        result <- ValidateOnLipidomics(result, scan, reference, parameter, collisionType)
      }
    }

    result$IsReferenceMatched <- isTRUE(result$IsPrecursorMzMatch) &&
      (!isTRUE(parameter$IsUseTimeForAnnotationScoring) || isTRUE(result$IsRtMatch)) &&
      (!isTRUE(parameter$IsUseCcsForAnnotationScoring) || isTRUE(result$IsCcsMatch)) &&
      (!useMs2 || isTRUE(result$IsSpectrumMatch))

    result$IsAnnotationSuggested <- isTRUE(result$IsPrecursorMzMatch) &&
      (!isTRUE(parameter$IsUseTimeForAnnotationScoring) || isTRUE(result$IsRtMatch)) &&
      (!isTRUE(parameter$IsUseCcsForAnnotationScoring) || isTRUE(result$IsCcsMatch)) &&
      !isTRUE(result$IsReferenceMatched)

    result
  }

  resolve_parameter_or_default <- function(query, annotatorContext) {
    parameter <- get_field(query, "Parameter")
    if (!is.null(parameter)) return(parameter)
    get_field(annotatorContext, "Parameter")
  }

  ensure_annotator_search_core <- function(annotatorContext) {
    fn <- get_field(annotatorContext, "SearchCore")
    if (is.function(fn)) return(fn)

    fn <- resolve_function("SearchCore")
    if (!is.null(fn)) return(fn)

    fn <- resolve_function("search_core")
    if (!is.null(fn)) return(fn)

    stop("SearchCore()/search_core() or annotatorContext$SearchCore is required.")
  }

  ensure_annotator_calculate_score <- function(annotatorContext) {
    fn <- get_field(annotatorContext, "CalculateScore")
    if (is.function(fn)) return(fn)

    fn <- resolve_function("LcmsMspAnnotatorCalculateScore")
    if (!is.null(fn)) return(fn)

    fn <- resolve_function("lcms_msp_annotator_calculate_score")
    if (!is.null(fn)) return(fn)

    stop("LcmsMspAnnotatorCalculateScore()/lcms_msp_annotator_calculate_score() or annotatorContext$CalculateScore is required.")
  }

  order_results_by_total_score <- function(results) {
    if (length(results) == 0) return(results)

    scores <- vapply(results, function(result) {
      total <- get_field(result, "TotalScore", NA_real_)
      if (is.null(total) || is.na(total)) -Inf else as.numeric(total)
    }, numeric(1))

    results[order(scores, decreasing = TRUE)]
  }

  # Translated from LcmsMspAnnotator.FindCandidatesCore
  FindCandidatesCore <- function(query, parameter, annotatorContext) {
    search_core_fn <- ensure_annotator_search_core(annotatorContext)
    calculate_score_fn <- ensure_annotator_calculate_score(annotatorContext)

    property <- get_field(query, "Property")
    candidates <- search_core_fn(property, parameter, annotatorContext)
    if (is.null(candidates) || length(candidates) == 0) {
      return(list())
    }

    results <- lapply(candidates, function(candidate) {
      calculate_score_fn(query, candidate, annotatorContext)
    })

    order_results_by_total_score(results)
  }

  # Translated from LcmsMspAnnotator.FindCandidates
  FindCandidates <- function(query, annotatorContext) {
    parameter <- resolve_parameter_or_default(query, annotatorContext)
    FindCandidatesCore(query, parameter, annotatorContext)
  }

  # Translated from LcmsMspAnnotator.Annotate
  Annotate <- function(query, annotatorContext) {
    parameter <- resolve_parameter_or_default(query, annotatorContext)
    results <- FindCandidatesCore(query, parameter, annotatorContext)

    if (length(results) == 0) {
      return(NULL)
    }

    results[[1]]
  }

  empty_ms_scan_match_results <- function() {
    list()
  }

  ensure_annotation_query_normalizer <- function() {
    fn <- resolve_function("GetNormalizedMSScanProperty")
    if (!is.null(fn)) return(fn)

    fn <- resolve_function("get_normalized_ms_scan_property")
    if (!is.null(fn)) return(fn)

    fn <- resolve_function("DataAccess.GetNormalizedMSScanProperty")
    if (!is.null(fn)) return(fn)

    stop("GetNormalizedMSScanProperty()/get_normalized_ms_scan_property() is required for AnnotationQuery.NormalizedScan.")
  }

  ensure_annotation_query_finder <- function(query) {
    finder <- get_field(query, ".finder")
    if (is.function(finder)) return(finder)

    if (is.environment(finder) && exists("FindCandidates", envir = finder, inherits = FALSE)) {
      fn <- get("FindCandidates", envir = finder, inherits = FALSE)
      if (is.function(fn)) return(fn)
    }

    finder <- get_field(query, "Finder")
    if (is.function(finder)) return(finder)

    stop("AnnotationQuery finder is required and must expose FindCandidates behavior.")
  }

  # Translated from AnnotationQuery.NormalizedScan getter
  AnnotationQueryNormalizedScan <- function(query) {
    cached <- get_field(query, ".normalizedScan")
    if (!is.null(cached)) {
      return(cached)
    }

    normalize_fn <- ensure_annotation_query_normalizer()
    normalized <- normalize_fn(get_field(query, "Scan"), get_field(query, "Parameter"))

    if (is.environment(query)) {
      assign(".normalizedScan", normalized, envir = query)
    } else if (is.list(query)) {
      query[[".normalizedScan"]] <- normalized
    }

    normalized
  }

  # Translated from AnnotationQuery.FindCandidates
  AnnotationQueryFindCandidates <- function(query, forceFind = FALSE) {
    finder <- get_field(query, ".finder")
    ignoreIsotopicPeak <- isTRUE(get_field(query, ".ignoreIsotopicPeak", TRUE))
    ionFeature <- get_field(query, "IonFeature")
    isMonoIsotopicIon <- isTRUE(get_field(ionFeature, "IsMonoIsotopicIon", FALSE))

    if (is.null(finder) || (!isTRUE(forceFind) && ignoreIsotopicPeak && !isMonoIsotopicIon)) {
      return(empty_ms_scan_match_results())
    }

    finder_fn <- ensure_annotation_query_finder(query)
    finder_fn(query)
  }

  # Translated from AnnotationQuery constructor
  AnnotationQuery <- function(property, scan, isotopes, ionFeature, parameter, finder, ignoreIsotopicPeak = TRUE) {
    if (is.null(property)) {
      stop("property cannot be NULL")
    }
    if (is.null(scan)) {
      stop("scan cannot be NULL")
    }
    if (is.null(parameter)) {
      stop("parameter cannot be NULL")
    }

    query <- new.env(parent = emptyenv())
    query$Property <- property
    query$Scan <- scan
    query$Isotopes <- isotopes
    query$Parameter <- parameter
    query$IonFeature <- ionFeature
    query$.finder <- finder
    query$.ignoreIsotopicPeak <- ignoreIsotopicPeak
    query$.normalizedScan <- NULL

    makeActiveBinding("NormalizedScan", function() {
      AnnotationQueryNormalizedScan(query)
    }, env = query)

    query$FindCandidates <- function(forceFind = FALSE) {
      AnnotationQueryFindCandidates(query, forceFind)
    }

    class(query) <- c("AnnotationQuery", class(query))
    query
  }

  STANDARD_ANNOTATION_PROCESS_NUMBER_OF_RESULTS <- 3L

  throw_if_cancellation_requested <- function(token = NULL) {
    if (is.null(token)) return(invisible(NULL))

    if (is.function(token)) {
      token()
      return(invisible(NULL))
    }

    throw_fn <- get_field(token, "ThrowIfCancellationRequested")
    if (is.function(throw_fn)) {
      throw_fn()
      return(invisible(NULL))
    }

    cancelled <- get_field(token, "IsCancellationRequested", FALSE)
    if (isTRUE(cancelled)) {
      stop("Operation cancelled.")
    }

    invisible(NULL)
  }

  create_reporter_from_range <- function(reportAction = NULL, start = 0, end = 1) {
    report_fn <- function(current, total) {
      if (!is.function(reportAction)) return(invisible(NULL))
      ratio <- if (is.null(total) || total <= 0) 1 else current / total
      scaled <- start + (end - start) * ratio
      reportAction(scaled)
      invisible(NULL)
    }

    list(Report = report_fn)
  }

  ensure_list_of_factories <- function(queryFactories) {
    if (is.null(queryFactories)) return(list())
    if (is.list(queryFactories) && !is.function(queryFactories)) return(queryFactories)
    list(queryFactories)
  }

  invoke_factory_create <- function(factory, chromPeakFeature, msdecResult, ms1Spectrum, peakCharacter, parameter) {
    create_fn <- get_field(factory, "Create")
    if (!is.function(create_fn)) {
      stop("Each query factory must provide a Create method/function.")
    }
    create_fn(chromPeakFeature, msdecResult, ms1Spectrum, peakCharacter, parameter)
  }

  invoke_factory_prepare_parameter <- function(factory) {
    prepare_fn <- get_field(factory, "PrepareParameter")
    if (!is.function(prepare_fn)) {
      stop("Each query factory must provide a PrepareParameter method/function.")
    }
    prepare_fn()
  }

  invoke_provider_load_ms_spectrum <- function(provider, index) {
    load_fn <- get_field(provider, "LoadMsSpectrumFromIndex")
    if (!is.function(load_fn)) {
      stop("provider must provide LoadMsSpectrumFromIndex().")
    }
    load_fn(index)
  }

  invoke_query_find_candidates <- function(query) {
    find_fn <- get_field(query, "FindCandidates")
    if (is.function(find_fn)) {
      return(find_fn())
    }
    AnnotationQueryFindCandidates(query)
  }

  invoke_evaluator_filter_by_threshold <- function(evaluator, candidates) {
    fn <- get_field(evaluator, "FilterByThreshold")
    if (!is.function(fn)) {
      stop("evaluator must provide FilterByThreshold().")
    }
    fn(candidates)
  }

  invoke_evaluator_select_top_n <- function(evaluator, results, n) {
    fn <- get_field(evaluator, "SelectTopN")
    if (!is.function(fn)) {
      stop("evaluator must provide SelectTopN().")
    }
    fn(results, n)
  }

  invoke_evaluator_is_reference_matched <- function(evaluator, representative) {
    fn <- get_field(evaluator, "IsReferenceMatched")
    if (!is.function(fn)) {
      stop("evaluator must provide IsReferenceMatched().")
    }
    fn(representative)
  }

  invoke_evaluator_is_annotation_suggested <- function(evaluator, representative) {
    fn <- get_field(evaluator, "IsAnnotationSuggested")
    if (!is.function(fn)) {
      stop("evaluator must provide IsAnnotationSuggested().")
    }
    fn(representative)
  }

  invoke_refer_refer <- function(refer, representative) {
    fn <- get_field(refer, "Refer")
    if (!is.function(fn)) {
      stop("refer must provide Refer().")
    }
    fn(representative)
  }

  invoke_match_results_add_results <- function(matchResults, results) {
    fn <- get_field(matchResults, "AddResults")
    if (!is.function(fn)) {
      stop("MatchResults must provide AddResults().")
    }
    fn(results)
  }

  set_molecule_ms_property_dispatch <- function(chromPeakFeature, reference, representative, suggested = FALSE) {
    fn <- NULL
    if (isTRUE(suggested)) {
      fn <- ensure_scoring_fn(c("SetMoleculeMsPropertyAsSuggested", "set_molecule_ms_property_as_suggested"))
    } else {
      fn <- ensure_scoring_fn(c("SetMoleculeMsProperty", "set_molecule_ms_property"))
    }
    if (is.null(fn)) {
      stop("Required DataAccess setter function is missing in the R environment.")
    }
    fn(chromPeakFeature, reference, representative)
  }

  get_msdec_results_at <- function(msdecResults, index) {
    results <- get_field(msdecResults, "MSDecResults")
    if (is.null(results)) {
      stop("msdecResults must expose MSDecResults.")
    }
    results[[index]]
  }

  # Translated from StandardAnnotationProcess constructors
  StandardAnnotationProcess <- function(queryFactories, evaluator, refer) {
    if (missing(queryFactories) || is.null(queryFactories)) {
      stop("queryFactories cannot be NULL")
    }
    if (is.null(evaluator)) {
      stop("evaluator cannot be NULL")
    }
    if (is.null(refer)) {
      stop("refer cannot be NULL")
    }

    process <- new.env(parent = emptyenv())
    process$.queryFactories <- ensure_list_of_factories(queryFactories)
    process$.evaluator <- evaluator
    process$.refer <- refer

    process$RunAnnotationAsync <- function(chromPeakFeatures, msdecResults, provider, numThreads = 1, reportAction = NULL, token = NULL) {
      RunAnnotationAsync(process, chromPeakFeatures, msdecResults, provider, numThreads, reportAction, token)
    }
    process$RunBySingleThreadAsync <- function(chromPeakFeatures, msdecResults, provider, token = NULL, reporter = create_reporter_from_range()) {
      RunBySingleThreadAsync(process, chromPeakFeatures, msdecResults, provider, token, reporter)
    }
    process$RunByMultiThreadAsync <- function(chromPeakFeatures, msdecResults, provider, numThreads, token = NULL, reporter = create_reporter_from_range()) {
      RunByMultiThreadAsync(process, chromPeakFeatures, msdecResults, provider, numThreads, token, reporter)
    }
    process$RunAnnotationCoreAsync <- function(chromPeakFeature, msdecResult, provider, collisionEnergy, token = NULL) {
      RunAnnotationCoreAsync(process, chromPeakFeature, msdecResult, provider, collisionEnergy, token)
    }
    process$SetAnnotationResult <- function(chromPeakFeature, query, msdecResult, collisionEnergy) {
      SetAnnotationResult(process, chromPeakFeature, query, msdecResult, collisionEnergy)
    }
    process$SetRepresentativeProperty <- function(chromPeakFeature) {
      SetRepresentativeProperty(process, chromPeakFeature)
    }

    class(process) <- c("StandardAnnotationProcess", class(process))
    process
  }

  # Translated from StandardAnnotationProcess.RunAnnotationAsync
  RunAnnotationAsync <- function(chromPeakFeatures, msdecResults, params, numThreads = 1) {

    if (numThreads <= 1) {
      RunBySingleThreadAsync(chromPeakFeatures, msdecResults, params)
    } else {
      RunByMultiThreadAsync(chromPeakFeatures, msdecResults, params, numThreads)
    }

  }

  # Translated from StandardAnnotationProcess.RunBySingleThreadAsync
  RunBySingleThreadAsync <- function(chromPeakFeatures, msdecResults, params) {
    total <- length(chromPeakFeatures)
    res <- rep(list(NA), length(mSDecResultCollections))

    for (i in seq_len(total)) {
      chromPeakFeature <- chromPeakFeatures[[i]]
      msdecResult <- msdecResults[[i]]
      res[[i]] <- RunAnnotationCoreAsync(chromPeakFeature, msdecResult, get_field(msdecResults, "CollisionEnergy"), params)
    }
  }

  # Translated from StandardAnnotationProcess.RunByMultiThreadAsync
  RunByMultiThreadAsync <- function(chromPeakFeatures, msdecResults, params, numThreads) {
    total <- length(chromPeakFeatures)
    res <- rep(list(NA), length(mSDecResultCollections))
    # TODO: make this as multi-cores
    for (i in seq_len(total)) {
      chromPeakFeature <- chromPeakFeatures[[i]]
      msdecResult <- msdecResults[[i]]
      res[[i]] <- RunAnnotationCoreAsync(chromPeakFeature, msdecResult, get_field(msdecResults, "CollisionEnergy"), params)
    }
  }

  # Translated from StandardAnnotationProcess.RunAnnotationCoreAsync
  RunAnnotationCoreAsync <- function(chromPeakFeature, msdecResult, collisionEnergy, params) {
    spectrum <- get_field(msdecResult, "Spectrum")
    if (!is.null(spectrum) && length(spectrum) > 0) {
      chromPeakFeature$MSDecResultIdUsed <- get_field(chromPeakFeature, "MasterPeakID")
    }

    #factories <- get_field(".queryFactories", list())
    #for (factory in factories) {
      # loaded <- invoke_provider_load_ms_spectrum(provider, get_field(chromPeakFeature, "MS1RawSpectrumIdTop"))
      # ms1Spectrum <- get_field(loaded, "Spectrum")
      # query <- invoke_factory_create(
      #   factory,
      #   chromPeakFeature,
      #   msdecResult,
      #   ms1Spectrum,
      #   get_field(chromPeakFeature, "PeakCharacter"),
      #   invoke_factory_prepare_parameter(factory)
      # )
      #throw_if_cancellation_requested(token)

      # TODO: neet to prepare some information for this steps
      query <- list(
        chromPeakFeature,
        msdecResult,
        ms1Spectrum,
        get_field(chromPeakFeature, "PeakCharacter"),
        params
      )

      SetAnnotationResult(query, collisionEnergy)
    #}

    #throw_if_cancellation_requested(token)
    #SetRepresentativeProperty(chromPeakFeature)
    #invisible(NULL)
  }

  # Translated from StandardAnnotationProcess.SetAnnotationResult
  SetAnnotationResult <- function(query, collisionEnergy) {


    #candidates <- invoke_query_find_candidates(query)

    candidates <- AnnotationQueryFindCandidates(query)



    # results <- invoke_evaluator_filter_by_threshold(get_field(process, ".evaluator"), candidates)
    # topResults <- invoke_evaluator_select_top_n(
    #   get_field(process, ".evaluator"),
    #   results,
    #   STANDARD_ANNOTATION_PROCESS_NUMBER_OF_RESULTS
    # )
    #
    # topResults <- lapply(topResults, function(result) {
    #   result$CollisionEnergy <- collisionEnergy
    #   result$SpectrumID <- get_field(msdecResult, "RawSpectrumID")
    #   result
    # })

    # invoke_match_results_add_results(get_field(chromPeakFeature, "MatchResults"), topResults)
    # invisible(NULL)
  }

  # Translated from StandardAnnotationProcess.SetRepresentativeProperty
  SetRepresentativeProperty <- function(process, chromPeakFeature) {
    representative <- get_field(get_field(chromPeakFeature, "MatchResults"), "Representative")
    evaluator <- get_field(process, ".evaluator")
    if (is.null(evaluator)) {
      return(invisible(NULL))
    }

    refer <- get_field(process, ".refer")
    if (invoke_evaluator_is_reference_matched(evaluator, representative)) {
      set_molecule_ms_property_dispatch(chromPeakFeature, invoke_refer_refer(refer, representative), representative, suggested = FALSE)
    } else if (invoke_evaluator_is_annotation_suggested(evaluator, representative)) {
      set_molecule_ms_property_dispatch(chromPeakFeature, invoke_refer_refer(refer, representative), representative, suggested = TRUE)
    }

    invisible(NULL)
  }



}

{

  get_field <- function(x, name, default = NULL) {
    if (is.null(x)) return(default)
    if (is.environment(x)) {
      if (exists(name, envir = x, inherits = FALSE)) return(get(name, envir = x))
      return(default)
    }
    if (is.list(x) && !is.null(x[[name]])) return(x[[name]])
    default
  }

  set_field <- function(x, name, value) {
    if (is.environment(x)) {
      assign(name, value, envir = x)
      return(x)
    }
    x[[name]] <- value
    x
  }

  call_method <- function(obj, candidates, ..., default = NULL) {
    for (nm in candidates) {
      fn <- get_field(obj, nm, NULL)
      if (is.function(fn)) {
        return(fn(...))
      }
    }
    default
  }

  as_rt_value <- function(chrom_xs) {
    rt <- get_field(chrom_xs, "RT", NULL)
    if (is.list(rt) || is.environment(rt)) {
      return(get_field(rt, "Value", NA_real_))
    }
    if (!is.null(rt)) return(rt)
    NA_real_
  }

  mass_reference_searcher <- function(db) {
    structure(
      list(
        Search = function(query) {
          lapply(db, identity)
        }
      ),
      class = "MassReferenceSearcher"
    )
  }

  mass_rt_reference_searcher <- function(db) {
    structure(
      list(
        Search = function(query) {
          lapply(db, identity)
        }
      ),
      class = "MassRtReferenceSearcher"
    )
  }

  fix_mass_tolerance <- function(ms1_tolerance, precursor_mz) {
    fn <- get0("MolecularFormulaUtility_FixMassTolerance", ifnotfound = NULL, inherits = TRUE)
    if (is.function(fn)) {
      return(fn(ms1_tolerance, precursor_mz))
    }
    ms1_tolerance
  }

  create_mass_rt_query <- function(precursor_mz, tol, rt, rt_tol) {
    list(type = "MassRtQuery", precursor_mz = precursor_mz, tolerance = tol, rt = rt, rt_tolerance = rt_tol)
  }

  create_mass_query <- function(precursor_mz, tol) {
    list(type = "MassQuery", precursor_mz = precursor_mz, tolerance = tol)
  }

  searcher <- function(self) {
    existing <- get_field(self, "searcher", NULL)
    if (!is.null(existing)) return(existing)
    db <- get_field(self, "db", list())
    created <- mass_reference_searcher(db)
    self <- set_field(self, "searcher", created)
    get_field(self, "searcher", created)
  }

  searcher_with_rt <- function(self) {
    existing <- get_field(self, "searcherWithRt", NULL)
    if (!is.null(existing)) return(existing)
    db <- get_field(self, "db", list())
    created <- mass_rt_reference_searcher(db)
    self <- set_field(self, "searcherWithRt", created)
    get_field(self, "searcherWithRt", created)
  }

  calculate_score <- function(self, query, reference) {
    scorer_obj <- get_field(self, "scorer", NULL)
    result <- call_method(scorer_obj, c("Score", "score"), query, reference)
    # this step will call one more calculate score located in another class [MsReferenceScorer]

    parameter <- get_field(query, "Parameter", get_field(self, "Parameter", NULL))
    use_time <- isTRUE(get_field(parameter, "IsUseTimeForAnnotationScoring", FALSE))

    is_precursor_mz_match <- isTRUE(get_field(result, "IsPrecursorMzMatch", FALSE))
    is_rt_match <- isTRUE(get_field(result, "IsRtMatch", FALSE))
    is_spectrum_match <- isTRUE(get_field(result, "IsSpectrumMatch", FALSE))

    is_reference_matched <- is_precursor_mz_match && (!use_time || is_rt_match) && is_spectrum_match
    result <- set_field(result, "IsReferenceMatched", is_reference_matched)

    is_annotation_suggested <- is_precursor_mz_match && (!use_time || is_rt_match) && !is_reference_matched
    result <- set_field(result, "IsAnnotationSuggested", is_annotation_suggested)

    result
  }

  search_core <- function(self, property, parameter) {
    precursor_mz <- get_field(property, "PrecursorMz", NA_real_)
    ms1_tolerance <- get_field(parameter, "Ms1Tolerance", NA_real_)
    fixed_tol <- fix_mass_tolerance(ms1_tolerance, precursor_mz)

    if (isTRUE(get_field(parameter, "IsUseTimeForAnnotationFiltering", FALSE))) {
      rt <- as_rt_value(get_field(property, "ChromXs", NULL))
      rt_tolerance <- get_field(parameter, "RtTolerance", NA_real_)
      query <- create_mass_rt_query(precursor_mz, fixed_tol, rt, rt_tolerance)
      return(call_method(searcher_with_rt(self), c("Search", "search"), query, default = list()))
    }

    query <- create_mass_query(precursor_mz, fixed_tol)
    call_method(searcher(self), c("Search", "search"), query, default = list())
  }

  find_candidates_core <- function(self, query, parameter) {
    property <- get_field(query, "Property", NULL)
    candidates <- search_core(self, property, parameter)

    scored <- lapply(candidates, function(candidate) {
      calculate_score(self, query, candidate)
    })

    if (length(scored) <= 1L) return(scored)

    scores <- vapply(scored, function(result) {
      as.numeric(get_field(result, "TotalScore", -Inf))
    }, numeric(1))

    scored[order(scores, decreasing = TRUE)]
  }



}

{
  get_field <- function(x, name, default = NULL) {
    if (is.null(x)) return(default)
    if (is.environment(x)) {
      if (exists(name, envir = x, inherits = FALSE)) return(get(name, envir = x))
      return(default)
    }
    if (is.list(x) && !is.null(x[[name]])) return(x[[name]])
    default
  }

  set_field <- function(x, name, value) {
    if (is.environment(x)) {
      assign(name, value, envir = x)
      return(x)
    }
    x[[name]] <- value
    x
  }

  call_method <- function(obj, candidates, ..., default = NULL) {
    for (nm in candidates) {
      fn <- get_field(obj, nm, NULL)
      if (is.function(fn)) {
        return(fn(...))
      }
    }
    default
  }

  as_rt_value <- function(chrom_xs) {
    rt <- get_field(chrom_xs, "RT", NULL)
    if (is.list(rt) || is.environment(rt)) {
      return(get_field(rt, "Value", NA_real_))
    }
    if (!is.null(rt)) return(rt)
    NA_real_
  }

  mass_reference_searcher <- function(db) {
    structure(
      list(
        Search = function(query) {
          lapply(db, identity)
        }
      ),
      class = "MassReferenceSearcher"
    )
  }

  mass_rt_reference_searcher <- function(db) {
    structure(
      list(
        Search = function(query) {
          lapply(db, identity)
        }
      ),
      class = "MassRtReferenceSearcher"
    )
  }

  fix_mass_tolerance <- function(ms1_tolerance, precursor_mz) {
    fn <- get0("MolecularFormulaUtility_FixMassTolerance", ifnotfound = NULL, inherits = TRUE)
    if (is.function(fn)) {
      return(fn(ms1_tolerance, precursor_mz))
    }
    ms1_tolerance
  }

  create_mass_rt_query <- function(precursor_mz, tol, rt, rt_tol) {
    list(type = "MassRtQuery", precursor_mz = precursor_mz, tolerance = tol, rt = rt, rt_tolerance = rt_tol)
  }

  create_mass_query <- function(precursor_mz, tol) {
    list(type = "MassQuery", precursor_mz = precursor_mz, tolerance = tol)
  }

  searcher <- function(self) {
    existing <- get_field(self, "searcher", NULL)
    if (!is.null(existing)) return(existing)
    db <- get_field(self, "db", list())
    created <- mass_reference_searcher(db)
    self <- set_field(self, "searcher", created)
    get_field(self, "searcher", created)
  }

  searcher_with_rt <- function(self) {
    existing <- get_field(self, "searcherWithRt", NULL)
    if (!is.null(existing)) return(existing)
    db <- get_field(self, "db", list())
    created <- mass_rt_reference_searcher(db)
    self <- set_field(self, "searcherWithRt", created)
    get_field(self, "searcherWithRt", created)
  }

  calculate_score <- function(self, query, reference) {
    scorer_obj <- get_field(self, "scorer", NULL)
    result <- call_method(scorer_obj, c("Score", "score"), query, reference)
    # this step will call one more calculate score located in another class [MsReferenceScorer]

    parameter <- get_field(query, "Parameter", get_field(self, "Parameter", NULL))
    use_time <- isTRUE(get_field(parameter, "IsUseTimeForAnnotationScoring", FALSE))

    is_precursor_mz_match <- isTRUE(get_field(result, "IsPrecursorMzMatch", FALSE))
    is_rt_match <- isTRUE(get_field(result, "IsRtMatch", FALSE))
    is_spectrum_match <- isTRUE(get_field(result, "IsSpectrumMatch", FALSE))

    is_reference_matched <- is_precursor_mz_match && (!use_time || is_rt_match) && is_spectrum_match
    result <- set_field(result, "IsReferenceMatched", is_reference_matched)

    is_annotation_suggested <- is_precursor_mz_match && (!use_time || is_rt_match) && !is_reference_matched
    result <- set_field(result, "IsAnnotationSuggested", is_annotation_suggested)

    result
  }

  search_core <- function(self, property, parameter) {
    precursor_mz <- get_field(property, "PrecursorMz", NA_real_)
    ms1_tolerance <- get_field(parameter, "Ms1Tolerance", NA_real_)
    fixed_tol <- fix_mass_tolerance(ms1_tolerance, precursor_mz)

    if (isTRUE(get_field(parameter, "IsUseTimeForAnnotationFiltering", FALSE))) {
      rt <- as_rt_value(get_field(property, "ChromXs", NULL))
      rt_tolerance <- get_field(parameter, "RtTolerance", NA_real_)
      query <- create_mass_rt_query(precursor_mz, fixed_tol, rt, rt_tolerance)
      return(call_method(searcher_with_rt(self), c("Search", "search"), query, default = list()))
    }

    query <- create_mass_query(precursor_mz, fixed_tol)
    call_method(searcher(self), c("Search", "search"), query, default = list())
  }

  find_candidates_core <- function(self, query, parameter) {
    property <- get_field(query, "Property", NULL)
    candidates <- search_core(self, property, parameter)

    scored <- lapply(candidates, function(candidate) {
      calculate_score(self, query, candidate)
    })

    if (length(scored) <= 1L) return(scored)

    scores <- vapply(scored, function(result) {
      as.numeric(get_field(result, "TotalScore", -Inf))
    }, numeric(1))

    scored[order(scores, decreasing = TRUE)]
  }

}

{
  get_field <- function(x, name, default = NULL) {
    if (is.null(x)) return(default)
    if (is.environment(x)) {
      if (exists(name, envir = x, inherits = FALSE)) return(get(name, envir = x))
      return(default)
    }
    if (is.list(x) && !is.null(x[[name]])) return(x[[name]])
    default
  }

  set_field <- function(x, name, value) {
    if (is.environment(x)) {
      assign(name, value, envir = x)
      return(x)
    }
    x[[name]] <- value
    x
  }

  call_method <- function(obj, candidates, ..., default = NULL) {
    for (nm in candidates) {
      fn <- get_field(obj, nm)
      if (is.function(fn)) {
        return(fn(...))
      }
    }
    default
  }

  has_flag <- function(option, flag) {
    if (is.null(option)) return(FALSE)
    if (is.character(option)) return(flag %in% option)
    if (is.list(option)) {
      flags <- get_field(option, "flags", NULL)
      if (is.character(flags)) return(flag %in% flags)
    }
    isTRUE(FALSE)
  }

  token_throw_if_cancellation_requested <- function(token) {
    if (is.null(token)) return(invisible(NULL))
    if (isTRUE(get_field(token, "isCancelled", FALSE))) {
      stop("Operation canceled", call. = FALSE)
    }
    throw_fn <- get_field(token, "ThrowIfCancellationRequested", NULL)
    if (is.function(throw_fn)) throw_fn()
    invisible(NULL)
  }

  report_progress <- function(progress, value) {
    if (is.null(progress)) return(invisible(NULL))
    if (is.function(progress)) {
      progress(value)
      return(invisible(NULL))
    }
    report_fn <- get_field(progress, "Report", NULL)
    if (is.function(report_fn)) report_fn(value)
    invisible(NULL)
  }

  make_reporter_from_length <- function(progress, initial, span) {
    list(
      Report = function(v) {
        report_progress(progress, initial + span * (v / 100))
      }
    )
  }

  # Translation of PeakPickProcess.Pick()
  pick <- function(self, file, provider, progress = NULL, token = NULL) {
    storage <- get_field(self, "storage")
    parameter <- get_field(storage, "Parameter")
    iupac_db <- get_field(storage, "IupacDatabase")

    peak_spotting <- get("PeakSpotting", envir = parent.frame())
    isotope_estimator <- get("IsotopeEstimator_Process", envir = parent.frame())
    chrom_collection_ctor <- get("ChromatogramPeakFeatureCollection", envir = parent.frame())

    chrom_peak_features <- peak_spotting(file, 0, 30)$Run(provider, parameter, progress, token)
    isotope_estimator(chrom_peak_features, parameter, iupac_db)
    chrom_collection_ctor(chrom_peak_features)
  }

  # Translation of SpectrumDeconvolutionProcess.Deconvolute()
  deconvolute <- function(self, provider, chromPeakFeatures, analysisFile, summaryDto, progress = NULL, token = NULL) {
    storage <- get_field(self, "storage")
    parameter <- get_field(storage, "Parameter")
    iupac_db <- get_field(storage, "IupacDatabase")

    msdec_result_collections <- list()
    initial_msdec <- 30.0
    max_msdec <- 30.0

    ce_list <- call_method(provider, c("LoadCollisionEnergyTargets", "load_collision_energy_targets"), default = list())
    summary <- get("ChromatogramPeaksDataSummary_ConvertFromDto", envir = parent.frame())(summaryDto)

    acquisition_type <- get_field(analysisFile, "AcquisitionType")
    if (identical(acquisition_type, "AIF")) {
      for (i in seq_along(ce_list)) {
        target_ce <- round(ce_list[[i]], 2)
        if (target_ce <= 0) {
          message("No correct CE information in AIF-MSDEC")
          next
        }
        max_msdec_aif <- max_msdec / length(ce_list)
        initial_msdec_aif <- initial_msdec + max_msdec_aif * (i - 1)

        ms2dec_ctor <- get("Ms2Dec", envir = parent.frame())
        results <- ms2dec_ctor(initial_msdec_aif, max_msdec_aif)$GetMS2DecResults(
          analysisFile, provider, chromPeakFeatures, parameter, summary, iupac_db, progress, token, target_ce
        )

        msdec_coll_ctor <- get("MSDecResultCollection", envir = parent.frame())
        msdec_result_collections[[length(msdec_result_collections) + 1]] <- msdec_coll_ctor(results, target_ce)
      }
    } else {
      target_ce <- if (length(ce_list) == 0) -1 else ce_list[[1]]
      ms2dec_ctor <- get("Ms2Dec", envir = parent.frame())
      results <- ms2dec_ctor(initial_msdec, max_msdec)$GetMS2DecResults(
        analysisFile, provider, chromPeakFeatures, parameter, summary, iupac_db, progress, token
      )

      msdec_coll_ctor <- get("MSDecResultCollection", envir = parent.frame())
      msdec_result_collections[[1]] <- msdec_coll_ctor(results, target_ce)
    }

    msdec_result_collections
  }

  # Translation of PeakAnnotationProcess.AnnotateAsync()
  annotate_async <- function(mSDecResultCollections, chromPeakFeatures, params, ncores = 1) {

    initial_annotation <- 60.0
    max_annotation <- 30.0
    max_annotation_local <- max_annotation / length(mSDecResultCollections)
    res <- rep(list(NA), length(mSDecResultCollections))

    for (index in seq_along(mSDecResultCollections)) {
      ce2msdecs <- mSDecResultCollections[[index]]
      initial_annotation_local <- initial_annotation + max_annotation_local * (index - 1)
      res[[i]] <- RunAnnotationAsync(chromPeakFeatures, ce2msdecs, params, ncores)
    }

    peak_character_estimator <- get("PeakCharacterEstimator", envir = parent.frame())
    argmin_fn <- get("argmin_collision_energy", envir = parent.frame())

    msdec_for_char <- NULL
    if (length(mSDecResultCollections) > 0) {
      min_idx <- argmin_fn(mSDecResultCollections)
      msdec_for_char <- get_field(mSDecResultCollections[[min_idx]], "MSDecResults")
    }

    peak_character_estimator(90, 10)$Process(file, provider, chromPeakFeatures, msdec_for_char, evaluator, parameter, progress)

    invisible(NULL)
  }

  # Translation of FileProcess.LoadPeakAndScans()
  load_peak_and_scans <- function(self, analysisFile, token = NULL) {
    load_peak_fn <- call_method(analysisFile, c("LoadChromatogramPeakFeatureCollectionAsync", "load_chromatogram_peak_feature_collection_async"), token)
    chrom_peak_features <- load_peak_fn

    clear_fn <- get_field(chrom_peak_features, "ClearMatchResultProperties")
    if (is.function(clear_fn)) clear_fn()

    deserialize_fn <- get("MSDecResultCollection_DeserializeAsync", envir = parent.frame())
    msdec_result_collections <- deserialize_fn(analysisFile, token)

    list(chromPeakFeatures = chrom_peak_features, mSDecResultCollections = msdec_result_collections)
  }

  # Translation of FileProcess.FindPeakAndScans()
  find_peak_and_scans <- function(self, analysisFile, provider, progress = NULL, token = NULL) {
    token_throw_if_cancellation_requested(token)
    message("Peak picking started")

    peak_pick_process <- get_field(self, "peakPickProcess")
    chrom_peak_features <- call_method(peak_pick_process, c("Pick", "pick"), analysisFile, provider, progress, token)

    summarizer <- get("ChromFeatureSummarizer_GetChromFeaturesSummary", envir = parent.frame())
    summary_dto <- summarizer(provider, get_field(chrom_peak_features, "Items"))
    analysisFile <- set_field(analysisFile, "ChromPeakFeaturesSummary", summary_dto)

    token_throw_if_cancellation_requested(token)
    message("Deconvolution started")

    spectrum_deconv_process <- get_field(self, "spectrumDeconvolutionProcess")
    msdec_result_collections <- call_method(
      spectrum_deconv_process,
      c("Deconvolute", "deconvolute"),
      provider,
      get_field(chrom_peak_features, "Items"),
      analysisFile,
      summary_dto,
      progress,
      token
    )

    list(chromPeakFeatures = chrom_peak_features, mSDecResultCollections = msdec_result_collections)
  }

  # Translation of FileProcess.SaveToFileAsync()
  save_to_file_async <- function(file, chromPeakFeatures, mSDecResultCollections) {
    call_method(chromPeakFeatures, c("SerializeAsync", "serialize_async"), file)

    if (length(mSDecResultCollections) == 1) {
      call_method(mSDecResultCollections[[1]], c("SerializeAsync", "serialize_async"), file)
    } else {
      deconv_path_list <- get_field(file, "DeconvolutionFilePathList")
      if (is.environment(deconv_path_list)) {
        clear_fn <- get_field(deconv_path_list, "Clear")
        if (is.function(clear_fn)) clear_fn()
      } else {
        file <- set_field(file, "DeconvolutionFilePathList", list())
      }

      for (mSDecResultCollection in mSDecResultCollections) {
        call_method(mSDecResultCollection, c("SerializeWithCEAsync", "serialize_with_ce_async"), file)
      }
    }

    invisible(NULL)
  }

  # Translation of FileProcess.RunAsync()
  run_async <- function(self, analysisFile, option, progress = NULL, token = NULL) {
    has_peak_spotting <- has_flag(option, "PeakSpotting")
    has_identification <- has_flag(option, "Identification")

    if (!has_peak_spotting && !has_identification) {
      return(invisible(NULL))
    }

    factory <- get_field(self, "factory")
    provider <- call_method(factory, c("Create", "create"), analysisFile)

    if (has_peak_spotting) {
      loaded <- find_peak_and_scans(self, analysisFile, provider, progress, token)
    } else {
      loaded <- load_peak_and_scans(self, analysisFile, token)
    }

    chromPeakFeatures <- loaded$chromPeakFeatures
    mSDecResultCollections <- loaded$mSDecResultCollections

    if (has_identification) {
      token_throw_if_cancellation_requested(token)
      message("Annotation started")

      peak_annotation_process <- get_field(self, "peakAnnotationProcess")
      call_method(
        peak_annotation_process,
        c("AnnotateAsync", "annotate_async"),
        analysisFile,
        mSDecResultCollections,
        get_field(chromPeakFeatures, "Items"),
        provider,
        progress,
        token
      )
    }

    token_throw_if_cancellation_requested(token)
    save_to_file_async(analysisFile, chromPeakFeatures, mSDecResultCollections)
    report_progress(progress, 100)

    invisible(NULL)
  }

  # Utility: index of minimum CollisionEnergy from a list of MSDecResultCollection-like objects.
  argmin_collision_energy <- function(collections) {
    ces <- vapply(collections, function(x) as.numeric(get_field(x, "CollisionEnergy", Inf)), numeric(1))
    which.min(ces)
  }

  # Example R representation of MSDIAL MsRefSearchParameterBase defaults.

  create_msref_search_parameter <- function(
    MassRangeBegin = 0,
    MassRangeEnd = 2000,
    RtTolerance = 100.0,
    RiTolerance = 100.0,
    CcsTolerance = 20.0,
    Ms1Tolerance = 0.01,
    Ms2Tolerance = 0.025,
    RelativeAmpCutoff = 0,
    AbsoluteAmpCutoff = 0,
    SquaredWeightedDotProductCutOff = 0.6 * 0.6,
    SquaredSimpleDotProductCutOff = 0.6 * 0.6,
    SquaredReverseDotProductCutOff = 0.8 * 0.8,
    MatchedPeaksPercentageCutOff = 0.25,
    AndromedaScoreCutOff = 0.1,
    TotalScoreCutoff = 0.8,
    MinimumSpectrumMatch = 3,
    IsUseTimeForAnnotationFiltering = FALSE,
    IsUseTimeForAnnotationScoring = FALSE,
    IsUseCcsForAnnotationFiltering = FALSE,
    IsUseCcsForAnnotationScoring = FALSE
  ) {
    list(
      MassRangeBegin = MassRangeBegin,
      MassRangeEnd = MassRangeEnd,
      RtTolerance = RtTolerance,
      RiTolerance = RiTolerance,
      CcsTolerance = CcsTolerance,
      Ms1Tolerance = Ms1Tolerance,
      Ms2Tolerance = Ms2Tolerance,
      RelativeAmpCutoff = RelativeAmpCutoff,
      AbsoluteAmpCutoff = AbsoluteAmpCutoff,
      SquaredWeightedDotProductCutOff = SquaredWeightedDotProductCutOff,
      SquaredSimpleDotProductCutOff = SquaredSimpleDotProductCutOff,
      SquaredReverseDotProductCutOff = SquaredReverseDotProductCutOff,
      WeightedDotProductCutOff = sqrt(SquaredWeightedDotProductCutOff),
      SimpleDotProductCutOff = sqrt(SquaredSimpleDotProductCutOff),
      ReverseDotProductCutOff = sqrt(SquaredReverseDotProductCutOff),
      MatchedPeaksPercentageCutOff = MatchedPeaksPercentageCutOff,
      AndromedaScoreCutOff = AndromedaScoreCutOff,
      TotalScoreCutoff = TotalScoreCutoff,
      MinimumSpectrumMatch = MinimumSpectrumMatch,
      IsUseTimeForAnnotationFiltering = IsUseTimeForAnnotationFiltering,
      IsUseTimeForAnnotationScoring = IsUseTimeForAnnotationScoring,
      IsUseCcsForAnnotationFiltering = IsUseCcsForAnnotationFiltering,
      IsUseCcsForAnnotationScoring = IsUseCcsForAnnotationScoring
    )
  }



}


# entry function to interface with internal function
lipid_annotation_core <- function(){
  load("/Volumes/ExtremeSSD/Lipidomics_projects/msdial_test_data/test_obj.rda")
  mSDecResultCollections <- list(dcl)
  chromPeakFeatures <- pai$raw

  # Example parameter list (same defaults as MsRefSearchParameterBase in C#).
  example_param <- create_msref_search_parameter()

  annotate_async(mSDecResultCollections, chromPeakFeatures, params = example_param, ncores = 1)


}

{
  # # process lbm2 file from ms-dail
  # db1 <- qs::qread("/Volumes/ExtremeSSD/Lipidomics_projects/Msdial_lbm2_all.qs")
  #
  # library(DBI)
  # library(RSQLite)
  #
  # # 2. Connect to (or create) the database file
  # # If "my_data.sqlite" doesn't exist, it is created in your working directory.
  # con <- dbConnect(RSQLite::SQLite(), "/Volumes/ExtremeSSD/Lipidomics_projects/lipidomics.sqlite")
  #
  # # 3. Define some data
  # as.data.frame(db1[["isotopes"]])-> df_isotope
  #
  # # 4. Write the table to the database
  # # name: the name of the table inside the database
  # # value: the R data frame you want to save
  # dbWriteTable(con, name = "isotopes", value = df_isotope, overwrite = TRUE)
  #
  # # 5. (Optional) Verify it worked by listing tables
  # print(dbListTables(con))
  #
  # # 6. Always close the connection when finished
  # dbDisconnect(con)
  #
  # db1[["spectra"]][["comment"]]<- NULL
  # db1[["spectra"]][["fragmentation_score"]]<- NULL
  #
  # as.data.frame(db1[["spectra"]])-> df_spectra
  # con <- dbConnect(RSQLite::SQLite(), "/Volumes/ExtremeSSD/Lipidomics_projects/lipidomics.sqlite")
  # dbWriteTable(con, name = "spectra", value = df_spectra, overwrite = TRUE)
  # print(dbListTables(con))
  # dbDisconnect(con)
  # rm(a1, a2, df_isotope, df_spectra, res)
  # gc()
  #
  # df_records <- data.frame(record_index = db1[["records"]][["record_index"]],
  #                          scan_id = db1[["records"]][["scan_id"]],
  #                          precursor_mz = db1[["records"]][["precursor_mz"]],
  #                          chrom_value = db1[["records"]][["chrom_value"]],
  #                          ion_mode = db1[["records"]][["ion_mode"]],
  #                          name = db1[["records"]][["name"]],
  #                          formula_string = db1[["records"]][["formula_string"]],
  #                          formula_exact_mass = db1[["records"]][["formula_exact_mass"]],
  #                          formula_m1_isotopic_abundance = db1[["records"]][["formula_m1_isotopic_abundance"]],
  #                          formula_m2_isotopic_abundance = db1[["records"]][["formula_m2_isotopic_abundance"]],
  #                          formula_c = db1[["records"]][["formula_c"]],
  #                          formula_n = db1[["records"]][["formula_n"]],
  #                          formula_h = db1[["records"]][["formula_h"]],
  #                          formula_o = db1[["records"]][["formula_o"]],
  #                          formula_s = db1[["records"]][["formula_s"]],
  #                          formula_p = db1[["records"]][["formula_p"]],
  #                          formula_f = db1[["records"]][["formula_f"]],
  #                          formula_cl = db1[["records"]][["formula_cl"]],
  #                          formula_br = db1[["records"]][["formula_br"]],
  #                          formula_i = db1[["records"]][["formula_i"]],
  #                          formula_si = db1[["records"]][["formula_si"]],
  #                          formula_c13 = db1[["records"]][["formula_c13"]],
  #                          formula_n15 = db1[["records"]][["formula_n15"]],
  #                          formula_h2 = db1[["records"]][["formula_h2"]],
  #                          formula_o18 = db1[["records"]][["formula_o18"]],
  #                          formula_s34 = db1[["records"]][["formula_s34"]],
  #                          formula_cl37 = db1[["records"]][["formula_cl37"]],
  #                          formula_br81 = db1[["records"]][["formula_br81"]],
  #                          formula_se = db1[["records"]][["formula_se"]],
  #                          ontology = db1[["records"]][["ontology"]],
  #                          smiles = db1[["records"]][["smiles"]],
  #                          inchi_key = db1[["records"]][["inchi_key"]],
  #                          adduct_accurate_mass = db1[["records"]][["adduct_accurate_mass"]],
  #                          adduct_xmer = db1[["records"]][["adduct_xmer"]],
  #                          adduct_name = db1[["records"]][["adduct_name"]],
  #                          adduct_charge = db1[["records"]][["adduct_charge"]],
  #                          adduct_ion_mode = db1[["records"]][["adduct_ion_mode"]],
  #                          adduct_format_check = db1[["records"]][["adduct_format_check"]],
  #                          adduct_m1_intensity = db1[["records"]][["adduct_m1_intensity"]],
  #                          adduct_m2_intensity = db1[["records"]][["adduct_m2_intensity"]],
  #                          adduct_is_radical = db1[["records"]][["adduct_is_radical"]],
  #                          adduct_is_included = db1[["records"]][["adduct_is_included"]],
  #                          collision_cross_section = db1[["records"]][["collision_cross_section"]],
  #                          quant_mass = db1[["records"]][["quant_mass"]],
  #                          compound_class = db1[["records"]][["compound_class"]],
  #                          minimum_peak_height = db1[["records"]][["minimum_peak_height"]],
  #                          is_target_molecule= db1[["records"]][["is_target_molecule"]],
  #                          spectrum_peak_count = db1[["records"]][["spectrum_peak_count"]],
  #                          isotopic_peak_count = db1[["records"]][["isotopic_peak_count"]])
  #
  # spectrum_mz_raw <- sapply(db1[["records"]][["spectrum_mz"]], function(x){serialize(x, connection = NULL)})
  # spectrum_mz_raw_df <- data.frame(spectrum_mz = I(spectrum_mz_raw))
  #
  # spectrum_intensity_raw <- sapply(db1[["records"]][["spectrum_intensity"]], function(x){serialize(x, connection = NULL)})
  # spectrum_intensity_raw_df <- data.frame(spectrum_intensity = I(spectrum_intensity_raw))
  #
  # formula_element2_count_raw <- sapply(db1[["records"]][["formula_element2_count"]], function(x){serialize(x, connection = NULL)})
  # formula_element2_count_raw_df <- data.frame(formula_element2_count = I(formula_element2_count_raw))
  #
  # isotope_mass_diff_raw <- sapply(db1[["records"]][["isotope_mass_diff"]], function(x){serialize(x, connection = NULL)})
  # isotope_mass_diff_raw_df <- data.frame(isotope_mass_diff = I(isotope_mass_diff_raw))
  #
  # isotope_relative_abundance_raw <- sapply(db1[["records"]][["isotope_relative_abundance"]], function(x){serialize(x, connection = NULL)})
  # isotope_relative_abundance_raw_df <- data.frame(isotope_relative_abundance = I(isotope_relative_abundance_raw))
  #
  #
  # df_records_all <- cbind(df_records,
  #                         spectrum_mz_raw_df,
  #                         spectrum_intensity_raw_df,
  #                         formula_element2_count_raw_df,
  #                         isotope_relative_abundance_raw_df,
  #                         isotope_mass_diff_raw_df)
  #
  # con <- dbConnect(RSQLite::SQLite(), "/Volumes/ExtremeSSD/Lipidomics_projects/lipidomics.sqlite")
  # dbWriteTable(con, name = "records", value = df_records_all, overwrite = TRUE)
  # print(dbListTables(con))
  # dbDisconnect(con)
  # rm(list = ls())
  # gc()


}
