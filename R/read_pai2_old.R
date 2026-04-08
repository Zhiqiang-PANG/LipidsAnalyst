# =============================================================================
# MSDIAL5 File Readers for .pai2 and .dcl files
# .pai2 — LZ4-compressed MessagePack of List<ChromatogramPeakFeature>
# .dcl  — Custom binary format for MSDecResult records
#
# .pai2 reader requires Python + pip packages: msgpack lz4
#   Install once: reticulate::py_install(c("msgpack", "lz4"))
# .dcl  reader is pure R (no extra packages).
# =============================================================================


# =============================================================================
# Key-to-name mapping for ChromatogramPeakFeature (MessagePack integer keys)
# =============================================================================
CHROM_PEAK_KEY_MAP <- c(
  "0"  = "ChromScanIdLeft",
  "1"  = "ChromScanIdTop",
  "2"  = "ChromScanIdRight",
  "3"  = "ChromXsLeft",          # ChromXs struct
  "4"  = "ChromXsTop",           # ChromXs struct (RT, RI, Drift, Mz)
  "5"  = "ChromXsRight",         # ChromXs struct
  "6"  = "PeakHeightLeft",
  "7"  = "PeakHeightTop",
  "8"  = "PeakHeightRight",
  "9"  = "PeakAreaAboveZero",
  "10" = "PeakAreaAboveBaseline",
  "11" = "MasterPeakID",
  "12" = "PeakID",
  "13" = "ParentPeakID",
  "15" = "SeekPointToDCLFile",
  "16" = "MS1RawSpectrumIdTop",
  "17" = "MS1AccumulatedMs1RawSpectrumIdTop",
  "18" = "MS2RawSpectrumID",
  "19" = "MS2RawSpectrumID2CE",  # dict<int, double>
  "20" = "ScanID",
  "22" = "IonMode",              # 0=Positive, 1=Negative
  "24" = "Spectrum",             # list of SpectrumPeak
  "25" = "Name",
  "26" = "Formula",
  "27" = "Ontology",
  "28" = "SMILES",
  "29" = "InChIKey",
  "30" = "AdductType",
  "31" = "CollisionCrossSection",
  "33" = "MSRawID2MspIDs",
  "35" = "TextDbIDs",
  "36" = "MSRawID2MspBasedMatchResult",
  "37" = "TextDbBasedMatchResult",
  "38" = "Comment",
  "39" = "PeakCharacter",
  "40" = "PeakShape",
  "41" = "FeatureFilterStatus",
  "42" = "DriftChromFeatures",
  "43" = "Mass",                 # PrecursorMz
  "44" = "MS1RawSpectrumIdLeft",
  "45" = "MS1RawSpectrumIdRight",
  "46" = "MS1AccumulatedMs1RawSpectrumIdLeft",
  "47" = "MS1AccumulatedMs1RawSpectrumIdRight",
  "49" = "MatchResults",
  "51" = "MSDecResultIdUsed",
  "52" = "Protein",
  "53" = "ProteinGroupID"
)


# =============================================================================
# Rename integer MessagePack keys to human-readable field names.
# Works recursively on nested lists.
# =============================================================================
rename_msgpack_keys <- function(x, key_map = CHROM_PEAK_KEY_MAP) {
  if (!is.list(x)) return(x)
  nms <- names(x)
  if (!is.null(nms)) {
    mapped <- key_map[nms]
    valid  <- !is.na(mapped)
    nms[valid] <- mapped[valid]
    names(x) <- nms
  }
  lapply(x, rename_msgpack_keys, key_map = key_map)
}


# =============================================================================
# Flatten a single ChromatogramPeakFeature to a one-row data.frame
# with the most useful scalar fields only.
# =============================================================================
peak_feature_to_row <- function(feat) {
  # ChromXsTop is itself a nested MessagePack object; its key 0 = RT value,
  # key 1 = RI value, key 2 = Drift value, key 3 = Mz value (all wrapped as
  # a two-element list {value, unit} where value is at sub-key 0).
  extract_chrom_value <- function(chromxs_obj, sub_key) {
    inner <- chromxs_obj[[as.character(sub_key)]]
    if (is.list(inner)) inner[[1]] else inner
  }

  chrom_top <- feat[["ChromXsTop"]]
  rt_val  <- if (!is.null(chrom_top)) tryCatch(extract_chrom_value(chrom_top, 0), error = function(e) NA) else NA
  ri_val  <- if (!is.null(chrom_top)) tryCatch(extract_chrom_value(chrom_top, 1), error = function(e) NA) else NA
  mz_val  <- if (!is.null(chrom_top)) tryCatch(extract_chrom_value(chrom_top, 3), error = function(e) NA) else NA

  adduct_obj <- feat[["AdductType"]]
  adduct_str <- if (is.character(adduct_obj)) adduct_obj else
    if (is.list(adduct_obj)) tryCatch(adduct_obj[["AdductIonName"]] %||% adduct_obj[["1"]], error = function(e) NA) else NA

  data.frame(
    MasterPeakID          = feat[["MasterPeakID"]]          %||% NA_integer_,
    PeakID                = feat[["PeakID"]]                 %||% NA_integer_,
    RT_min                = rt_val,
    RI                    = ri_val,
    PrecursorMz           = feat[["Mass"]]                   %||% mz_val,
    PeakHeightTop         = feat[["PeakHeightTop"]]          %||% NA_real_,
    PeakAreaAboveZero     = feat[["PeakAreaAboveZero"]]      %||% NA_real_,
    IonMode               = feat[["IonMode"]]                %||% NA_integer_,
    Name                  = feat[["Name"]]                   %||% NA_character_,
    Ontology              = feat[["Ontology"]]               %||% NA_character_,
    SMILES                = feat[["SMILES"]]                 %||% NA_character_,
    InChIKey              = feat[["InChIKey"]]               %||% NA_character_,
    AdductType            = adduct_str                       %||% NA_character_,
    Comment               = feat[["Comment"]]                %||% NA_character_,
    MSDecResultIdUsed     = feat[["MSDecResultIdUsed"]]      %||% NA_integer_,
    MS2RawSpectrumID      = feat[["MS2RawSpectrumID"]]       %||% NA_integer_,
    stringsAsFactors      = FALSE
  )
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b


# =============================================================================
# read_pai2()
# Reads a .pai2 file (LZ4-compressed MessagePack) using Python via reticulate.
#
# Returns a list with two elements:
#   $raw     — the raw nested R list (all MessagePack fields with integer keys
#              renamed to human-readable names where known)
#   $peaks   — a data.frame of scalar fields for each peak feature
#
# Parameters:
#   path          — path to the .pai2 file
#   rename_keys   — if TRUE, rename integer MessagePack keys to field names
#   as_dataframe  — if TRUE, also return $peaks data.frame
# =============================================================================
read_pai2 <- function(path, rename_keys = TRUE, as_dataframe = TRUE) {
  if (!file.exists(path)) stop("File not found: ", path)

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required. Install with: install.packages('reticulate')")
  }
  library(reticulate)

  # Check Python package is available
  tryCatch(import("msgpack"), error = function(e)
    stop("Python 'msgpack' not found. Run: reticulate::py_install(c('msgpack', 'lz4'))"))

  py_run_string("
import struct, msgpack

def _unpackb(data):
    try:
        return msgpack.unpackb(data, strict_map_key=False)
    except TypeError:
        try:
            return msgpack.unpackb(data, raw=False)
        except TypeError:
            return msgpack.unpackb(data)

def read_pai2(path):
    with open(path, 'rb') as f:
        data = f.read()

    first = data[0]

    # Detect MessagePack-CSharp LZ4 ext type (type code 99)
    # ext 8  : 0xC7 | 1-byte-len  | typecode | payload
    # ext 16 : 0xC8 | 2-byte-len-BE | typecode | payload
    # ext 32 : 0xC9 | 4-byte-len-BE | typecode | payload
    if first == 0xC7:
        type_code = data[2]
        payload_start = 3
    elif first == 0xC8:
        type_code = data[3]
        payload_start = 4
    elif first == 0xC9:
        type_code = data[5]
        payload_start = 6
    else:
        return _unpackb(data)

    if type_code != 99:
        raise ValueError(f'Unknown MessagePack ext type code: {type_code}')

    payload = data[payload_start:]

    # MessagePack-CSharp ext(99) has multiple historical payload layouts.
    # Variant A (observed in MSDIAL): [0xD2][int32 BE orig_size][lz4 block]
    # Variant B:                     [int32 LE orig_size][lz4 block]
    # Variant C:                     [int32 BE orig_size][lz4 block]
    if len(payload) >= 5 and payload[0] == 0xD2:
      orig_size = struct.unpack('>I', payload[1:5])[0]
      lz4_data = payload[5:]
    elif len(payload) >= 5 and payload[0] == 0xCE:
      orig_size = struct.unpack('>I', payload[1:5])[0]
      lz4_data = payload[5:]
    else:
      orig_size = struct.unpack('<I', payload[0:4])[0]
      lz4_data = payload[4:]
      if orig_size <= 0 or orig_size > 1000000000:
        be_size = struct.unpack('>I', payload[0:4])[0]
        if 0 < be_size <= 1000000000:
          orig_size = be_size

    errors = []

    # Strategy 1: lz4.block, raw block, no stored-size prefix
    try:
        import lz4.block
        try:
            return _unpackb(lz4.block.decompress(lz4_data,
                                                  uncompressed_size=orig_size,
                                                  store_size=False))
        except TypeError:
            return _unpackb(lz4.block.decompress(lz4_data,
                                                  uncompressed_size=orig_size))
    except Exception as e:
        errors.append(f'block/no-prefix: {e}')

    # Strategy 2: lz4.frame on data after the detected prefix
    # (K4os.Compression.LZ4 may produce LZ4 frame format)
    try:
        import lz4.frame
        return _unpackb(lz4.frame.decompress(lz4_data))
    except Exception as e:
        errors.append(f'frame/skip4: {e}')

    # Strategy 3: lz4.frame on entire payload (no 4-byte skip)
    try:
        import lz4.frame
        return _unpackb(lz4.frame.decompress(payload))
    except Exception as e:
        errors.append(f'frame/full-payload: {e}')

    # Strategy 4: lz4.block with store_size=True (size prefix already in lz4_data)
    try:
        import lz4.block
        return _unpackb(lz4.block.decompress(lz4_data))
    except Exception as e:
        errors.append(f'block/store-size-true: {e}')

    diag = (
        f'All LZ4 strategies failed for {path}\\n'
        f'  first_byte=0x{first:02X}, ext_type={type_code}, '
        f'payload_start={payload_start}, orig_size={orig_size}\\n'
        f'  payload[:20] hex = {payload[:20].hex()}\\n'
        f'  lz4_data[:20] hex = {lz4_data[:20].hex()}\\n'
        f'  Errors: ' + ' | '.join(errors)
    )
    raise RuntimeError(diag)
")

  result <- py$read_pai2(path)

  # Convert Python object to R list
  result_r <- py_to_r(result)

  if (rename_keys) {
    result_r <- rename_msgpack_keys(result_r)
  }

  out <- list(raw = result_r)

  if (as_dataframe) {
    rows <- lapply(result_r, function(feat) {
      tryCatch(peak_feature_to_row(feat), error = function(e) NULL)
    })
    rows <- Filter(Negate(is.null), rows)
    out$peaks <- do.call(rbind, rows)
    rownames(out$peaks) <- NULL
  }

  out
}


# =============================================================================
# read_dcl()
# Reads a .dcl file (custom MSDIAL binary format) in pure R.
#
# Returns a list of MSDecResult records, each a list with:
#   $ScanID, $RawSpectrumID, $PrecursorMz, $IonMode,
#   $RT_min, $RI, $DriftTime_ms, $ChromXsMz,
#   $ModelPeakMz, $ModelPeakHeight, $ModelPeakArea,
#   $IntegratedHeight, $IntegratedArea,
#   $AmplitudeScore, $ModelPeakPurity, $ModelPeakQuality,
#   $SignalNoiseRatio, $EstimatedNoise,
#   $Spectrum        — data.frame(Mass, Intensity, PeakQuality)
#   $ModelPeakChromatogram — data.frame(ID, ChromX, Mass, Intensity)
#   $ModelMasses     — numeric vector
#   (if annotation info is present, also $MspID, $MspBasedMatchResult, etc.)
# =============================================================================
read_dcl <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)

  con <- file(path, "rb")
  on.exit(close(con))

  # ---- Header ---------------------------------------------------------------
  # Bytes 0-1: "DC" ASCII
  magic <- readChar(con, nchars = 2, useBytes = TRUE)
  if (magic != "DC") stop("Not a valid .dcl file (expected 'DC' header, got '", magic, "')")

  # Bytes 2-5: version (int32 LE)
  version <- readBin(con, "integer", n = 1, size = 4, endian = "little")
  if (version != 1L) stop("Unsupported DCL version: ", version, " (only version 1 is supported)")

  # Byte 6: isAnnotationInfoIncluded (bool = 1 byte)
  annot_flag <- readBin(con, "raw", n = 1)
  is_annotation <- as.integer(annot_flag) != 0L

  # Bytes 7-10: totalPeakNumber (int32 LE)
  n_peaks <- readBin(con, "integer", n = 1, size = 4, endian = "little")

  # Seek pointers: n_peaks * int64 LE (we read them but rely on sequential reading)
  seek_pts <- readBin(con, "integer", n = n_peaks * 2L, size = 4, endian = "little")
  # Pack pairs of int32 into int64 values (informational; not needed for sequential read)
  seek_pointers <- seek_pts[seq(1, length(seek_pts), 2)] +
    seek_pts[seq(2, length(seek_pts), 2)] * 2^32

  # ---- Read each MSDecResult record -----------------------------------------
  results <- vector("list", n_peaks)

  for (i in seq_len(n_peaks)) {
    r <- list()

    # -- Scan data (60 bytes) --
    # int64 SeekPoint (8), int32 ScanID (4), int32 RawSpectrumID (4),
    # float64 PrecursorMz (8), int32 IonMode (4),
    # float64 RT (8), float64 RI (8), float64 Drift (8), float64 Mz (8)
    seek_lo  <- readBin(con, "integer", n = 1, size = 4, endian = "little")
    seek_hi  <- readBin(con, "integer", n = 1, size = 4, endian = "little")
    r$SeekPoint      <- seek_lo + seek_hi * 2^32
    r$ScanID         <- readBin(con, "integer", n = 1, size = 4, endian = "little")
    r$RawSpectrumID  <- readBin(con, "integer", n = 1, size = 4, endian = "little")
    r$PrecursorMz    <- readBin(con, "double",  n = 1, size = 8, endian = "little")
    r$IonMode        <- readBin(con, "integer", n = 1, size = 4, endian = "little")
    r$RT_min         <- readBin(con, "double",  n = 1, size = 8, endian = "little")
    r$RI             <- readBin(con, "double",  n = 1, size = 8, endian = "little")
    r$DriftTime_ms   <- readBin(con, "double",  n = 1, size = 8, endian = "little")
    r$ChromXsMz      <- readBin(con, "double",  n = 1, size = 8, endian = "little")

    # -- Quant data (40 bytes) --
    # float64 x 5: ModelPeakMz, ModelPeakHeight, ModelPeakArea,
    #              IntegratedHeight, IntegratedArea
    r$ModelPeakMz      <- readBin(con, "double", n = 1, size = 8, endian = "little")
    r$ModelPeakHeight  <- readBin(con, "double", n = 1, size = 8, endian = "little")
    r$ModelPeakArea    <- readBin(con, "double", n = 1, size = 8, endian = "little")
    r$IntegratedHeight <- readBin(con, "double", n = 1, size = 8, endian = "little")
    r$IntegratedArea   <- readBin(con, "double", n = 1, size = 8, endian = "little")

    # -- Score data (20 bytes) --
    # float32 x 5: AmplitudeScore, ModelPeakPurity, ModelPeakQuality,
    #              SignalNoiseRatio, EstimatedNoise
    r$AmplitudeScore    <- readBin(con, "double", n = 1, size = 4, endian = "little")
    r$ModelPeakPurity   <- readBin(con, "double", n = 1, size = 4, endian = "little")
    r$ModelPeakQuality  <- readBin(con, "double", n = 1, size = 4, endian = "little")
    r$SignalNoiseRatio  <- readBin(con, "double", n = 1, size = 4, endian = "little")
    r$EstimatedNoise    <- readBin(con, "double", n = 1, size = 4, endian = "little")

    # -- Counters (12 bytes) --
    # int32: spectraNumber, datapointNumber, modelMassCount
    n_spectra    <- readBin(con, "integer", n = 1, size = 4, endian = "little")
    n_datapoints <- readBin(con, "integer", n = 1, size = 4, endian = "little")
    n_model_masses <- readBin(con, "integer", n = 1, size = 4, endian = "little")

    # -- Spectrum data (n_spectra * 20 bytes) --
    # Per peak: float64 Mass (8), float64 Intensity (8), int32 PeakQuality (4)
    if (n_spectra > 0L) {
      spec_raw  <- readBin(con, "raw", n = n_spectra * 20L)
      spec_mass <- readBin(spec_raw[seq(1,  n_spectra*20, 20) + 0:7 - 1],  # vectorised below
                           "double", n = n_spectra, size = 8, endian = "little")
      spec_int  <- readBin(spec_raw, "double",  n = n_spectra * 20L / 8L, endian = "little")
      spec_qual <- readBin(spec_raw, "integer", n = n_spectra * 20L / 4L, endian = "little")
      # Each record is 20 bytes = 8+8+4 -> indices in the flat vectors above:
      mass_idx  <- seq(1, by = 20/8, length.out = n_spectra)   # every 2.5 doubles -> use raw
      # Re-extract cleanly
      spec_mass <- numeric(n_spectra)
      spec_int2 <- numeric(n_spectra)
      spec_qual2 <- integer(n_spectra)
      for (j in seq_len(n_spectra)) {
        offset <- (j - 1L) * 20L
        spec_mass[j]  <- readBin(spec_raw[(offset + 1L):(offset + 8L)],  "double",  size = 8, endian = "little")
        spec_int2[j]  <- readBin(spec_raw[(offset + 9L):(offset + 16L)], "double",  size = 8, endian = "little")
        spec_qual2[j] <- readBin(spec_raw[(offset + 17L):(offset + 20L)],"integer", size = 4, endian = "little")
      }
      r$Spectrum <- data.frame(Mass = spec_mass, Intensity = spec_int2, PeakQuality = spec_qual2)
    } else {
      r$Spectrum <- data.frame(Mass = numeric(0), Intensity = numeric(0), PeakQuality = integer(0))
    }

    # -- Base peak chromatogram (n_datapoints * 28 bytes) --
    # Per point: int32 ID (4), float64 ChromX (8), float64 Mass (8), float64 Intensity (8)
    if (n_datapoints > 0L) {
      bp_raw <- readBin(con, "raw", n = n_datapoints * 28L)
      bp_id   <- integer(n_datapoints)
      bp_cx   <- numeric(n_datapoints)
      bp_mass <- numeric(n_datapoints)
      bp_int  <- numeric(n_datapoints)
      for (j in seq_len(n_datapoints)) {
        offset <- (j - 1L) * 28L
        bp_id[j]   <- readBin(bp_raw[(offset + 1L):(offset + 4L)],  "integer", size = 4, endian = "little")
        bp_cx[j]   <- readBin(bp_raw[(offset + 5L):(offset + 12L)], "double",  size = 8, endian = "little")
        bp_mass[j] <- readBin(bp_raw[(offset + 13L):(offset + 20L)],"double",  size = 8, endian = "little")
        bp_int[j]  <- readBin(bp_raw[(offset + 21L):(offset + 28L)],"double",  size = 8, endian = "little")
      }
      r$ModelPeakChromatogram <- data.frame(ID = bp_id, ChromX = bp_cx, Mass = bp_mass, Intensity = bp_int)
    } else {
      r$ModelPeakChromatogram <- data.frame(ID = integer(0), ChromX = numeric(0), Mass = numeric(0), Intensity = numeric(0))
    }

    # -- Model masses (n_model_masses * 8 bytes) --
    if (n_model_masses > 0L) {
      r$ModelMasses <- readBin(con, "double", n = n_model_masses, size = 8, endian = "little")
    } else {
      r$ModelMasses <- numeric(0)
    }

    # -- Annotation info (GC-MS only, when isAnnotationInfoIncluded = TRUE) --
    if (is_annotation) {
      ann_raw <- readBin(con, "raw", n = 73L)

      annot <- list()
      annot$MspID              <- readBin(ann_raw[1:4],   "integer", size = 4, endian = "little")
      annot$MspIDWhenOrdered   <- readBin(ann_raw[5:8],   "integer", size = 4, endian = "little")
      annot$TotalScore         <- readBin(ann_raw[9:12],  "double",  size = 4, endian = "little")
      annot$WeightedDotProduct <- readBin(ann_raw[13:16], "double",  size = 4, endian = "little")
      annot$SimpleDotProduct   <- readBin(ann_raw[17:20], "double",  size = 4, endian = "little")
      annot$ReverseDotProduct  <- readBin(ann_raw[21:24], "double",  size = 4, endian = "little")
      annot$MatchedPeaksCount  <- readBin(ann_raw[25:28], "double",  size = 4, endian = "little")
      annot$MatchedPeaksPerc   <- readBin(ann_raw[29:32], "double",  size = 4, endian = "little")
      annot$EssentialFragScore <- readBin(ann_raw[33:36], "double",  size = 4, endian = "little")
      annot$RtSimilarity       <- readBin(ann_raw[37:40], "double",  size = 4, endian = "little")
      annot$RiSimilarity       <- readBin(ann_raw[41:44], "double",  size = 4, endian = "little")
      annot$CcsSimilarity      <- readBin(ann_raw[45:48], "double",  size = 4, endian = "little")
      annot$IsotopeSimilarity  <- readBin(ann_raw[49:52], "double",  size = 4, endian = "little")
      annot$AccurateMassSim    <- readBin(ann_raw[53:56], "double",  size = 4, endian = "little")
      annot$LibraryID          <- readBin(ann_raw[57:60], "integer", size = 4, endian = "little")
      annot$LibraryIDWhenOrdered <- readBin(ann_raw[61:64], "integer", size = 4, endian = "little")
      annot$IsPrecursorMzMatch <- as.logical(ann_raw[65])
      annot$IsSpectrumMatch    <- as.logical(ann_raw[66])
      annot$IsRtMatch          <- as.logical(ann_raw[67])
      annot$IsRiMatch          <- as.logical(ann_raw[68])
      annot$IsCcsMatch         <- as.logical(ann_raw[69])
      annot$IsLipidClassMatch  <- as.logical(ann_raw[70])
      annot$IsLipidChainsMatch <- as.logical(ann_raw[71])
      annot$IsLipidPositionMatch <- as.logical(ann_raw[72])
      annot$IsOtherLipidMatch  <- as.logical(ann_raw[73])

      id_count  <- readBin(con, "integer", n = 1, size = 4, endian = "little")
      if (id_count > 0L) {
        annot$MspIDs <- readBin(con, "integer", n = id_count, size = 4, endian = "little")
      } else {
        annot$MspIDs <- integer(0)
      }
      r$MspBasedMatchResult <- annot
    }

    results[[i]] <- r
  }

  results
}


# =============================================================================
# Convenience: convert list of MSDecResults to a summary data.frame
# =============================================================================
dcl_to_dataframe <- function(dcl_list) {
  rows <- lapply(seq_along(dcl_list), function(i) {
    r <- dcl_list[[i]]
    data.frame(
      Index            = i,
      ScanID           = r$ScanID,
      RawSpectrumID    = r$RawSpectrumID,
      PrecursorMz      = r$PrecursorMz,
      RT_min           = r$RT_min,
      IonMode          = r$IonMode,
      IntegratedHeight = r$IntegratedHeight,
      IntegratedArea   = r$IntegratedArea,
      AmplitudeScore   = r$AmplitudeScore,
      SignalNoiseRatio = r$SignalNoiseRatio,
      nSpectrumPeaks   = nrow(r$Spectrum),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


# =============================================================================
# Example usage
# =============================================================================
if (FALSE) {

  # --- Reading a .pai2 file ---------------------------------------------------
  pai <- read_pai2("path/to/sample.pai2")

  # Summary data.frame (one row per peak feature)
  head(pai$peaks)

  # Full nested raw list for the first peak
  str(pai$raw[[1]], max.level = 2)

  # Get RT and precursor m/z for all peaks
  peaks_df <- pai$peaks
  hist(peaks_df$RT_min, main = "RT distribution", xlab = "RT (min)")
  plot(peaks_df$RT_min, peaks_df$PrecursorMz, pch = ".", xlab = "RT (min)", ylab = "m/z")


  # --- Reading a .dcl file ----------------------------------------------------
  dcl <- read_dcl("path/to/sample.dcl")

  # How many deconvoluted spectra?
  length(dcl)

  # Summary table
  dcl_df <- dcl_to_dataframe(dcl)
  head(dcl_df)

  # Spectrum peaks for the 5th deconvoluted result
  dcl[[5]]$Spectrum

  # Model peak chromatogram for result 1
  plot(dcl[[1]]$ModelPeakChromatogram$ChromX,
       dcl[[1]]$ModelPeakChromatogram$Intensity,
       type = "l", xlab = "RT (min)", ylab = "Intensity")

  # --- Reading multiple .dcl files for collision energies ---------------------
  dcl_files <- list.files("path/to/project/", pattern = "\\.dcl$", full.names = TRUE)
  all_dcl   <- lapply(dcl_files, read_dcl)
  names(all_dcl) <- basename(dcl_files)
}
