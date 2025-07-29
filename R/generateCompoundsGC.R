#' Generate compounds for GC-EI-HRMS workflows
#'
#' This method is adapted from generateCompoundsLibrary() to support GC-EI full scan data.
#' It uses the "MS" slot and optionally filters by Retention Index (RI).
#' 
#' @inheritParams generateCompoundsLibrary
#' @param RI Optional numeric vector of retention indices for each feature group.
#' @param RItol Retention index tolerance (default: 10 units).
#' @export

# Edited 25 Jul 2025 - Theodora - For GC-EI compound library search
#' @title generateCompoundsGC
#' @description Library search for GC-EI data with optional retention index matching
setGeneric("generateCompoundsGC", function(fGroups, MSPeakLists, MSLibrary,
                                           minSim = 0.75, minAnnSim = minSim, absMzDev = 0.002,
                                           adduct = NULL, checkIons = "adduct",
                                           specSimParams = getDefSpecSimParams(),
                                           specSimParamsLib = getDefSpecSimParams(),
                                           RIalkaneFile = NULL, RItol = 10) {
  standardGeneric("generateCompoundsGC")
})

setMethod("generateCompoundsGC", "featureGroups", function(fGroups, MSPeakLists, MSLibrary,
                                                           minSim = 0.75, minAnnSim = minSim, absMzDev = 0.002,
                                                           adduct = NULL, checkIons = "adduct",
                                                           specSimParams = getDefSpecSimParams(),
                                                           specSimParamsLib = getDefSpecSimParams(),
                                                           RIalkaneFile = NULL, RItol = 10)
{
  ac <- checkmate::makeAssertCollection()
  checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
  checkmate::assertClass(MSLibrary, "MSLibrary", add = ac)
  checkmate::assertChoice(checkIons, c("adduct", "polarity", "none"), add = ac)
  assertSpecSimParams(specSimParams, add = ac)
  assertSpecSimParams(specSimParamsLib, add = ac)
  checkmate::reportAssertions(ac)

  if (length(fGroups) == 0)
    return(compounds(algorithm = "library"))

  if (checkIons != "none")
    adduct <- checkAndToAdduct(adduct, fGroups)

  gCount <- length(fGroups)
  gInfo <- groupInfo(fGroups)
  annTbl <- annotations(fGroups)
  libRecs <- records(MSLibrary)
  libSpecs <- spectra(MSLibrary)

  # 🟡 GC-EI fallback: skip precursor filtering if largely missing
  isEIlib <- all(is.na(libRecs$PrecursorMZ)) || all(is.na(libRecs$Precursor_type))

  # 🔵 NEW: Extract RI from Comments$retention_index if present
  if (!"RI" %in% colnames(libRecs)) {
    if (!is.null(libRecs$Comments$retention_index)) {
      libRecs[, RI := suppressWarnings(
        as.numeric(sub(".*?(\\d+\\.?\\d*).*", "\\1", Comments$retention_index))
      )]
    }
  }

  # 🟡 Relaxed filters for GC-EI library
  if (!isEIlib) {
    libRecs <- libRecs[!is.na(PrecursorMZ) & !is.na(SMILES) & !is.na(InChI) & !is.na(InChIKey) & !is.na(formula)]
    if (checkIons == "adduct")
      libRecs <- libRecs[!is.na(Precursor_type)]
    else if (checkIons == "polarity")
      libRecs <- libRecs[!is.na(Ion_mode)]
  } else {
    libRecs <- libRecs[!is.na(InChIKey) & !is.na(formula)]  # 🔵 relaxed
  }

  if (!is.null(adduct) && !isEIlib)
    libRecs <- libRecs[Precursor_type == as.character(adduct)]

  if (!is.null(checkIons) && checkIons == "none")
    adduct <- NULL

  # 🔵 Optional: calculate sample RI from alkane file
  sampleRI <- NULL
  if (!is.null(RIalkaneFile)) {
    checkmate::assertDataFrame(RIalkaneFile, min.rows = 2, min.cols = 2)
    checkmate::assertNames(colnames(RIalkaneFile), must.include = c("Num", "RT"))

    calcRI <- function(rt) {
      alk <- RIalkaneFile
      if (any(is.na(rt))) return(NA_real_)
      lower <- max(which(alk$RT <= rt))
      upper <- min(which(alk$RT > rt))
      if (is.infinite(lower) || is.infinite(upper)) return(NA_real_)
      n <- alk$Num[lower]
      t1 <- alk$RT[lower]; t2 <- alk$RT[upper]
      100 * (n + (rt - t1) / (t2 - t1))
    }

    sampleRI <- vapply(gInfo$rt, calcRI, numeric(1))
    names(sampleRI) <- rownames(gInfo)
  }

  cacheDB <- openCacheDBScope()
  baseHash <- makeHash(minSim, minAnnSim, absMzDev, adduct, checkIons, specSimParams, specSimParamsLib)
  setHash <- makeHash(fGroups, MSPeakLists, MSLibrary, baseHash)
  cachedSet <- loadCacheSet("compoundsLibrary", setHash, cacheDB)
  resultHashes <- vector("character", gCount)
  resultHashCount <- 0

  printf("Processing %d feature groups with a library of %d records...\n", gCount, nrow(libRecs))

  compList <- withProg(length(fGroups), FALSE, sapply(names(fGroups), function(grp)
  {
    doProgress()
    if (is.null(MSPeakLists[[grp]]) || is.null(MSPeakLists[[grp]][["MS"]])) return(NULL)
    spec <- MSPeakLists[[grp]][["MS"]]
    if (is.null(spec) || nrow(spec) == 0) return(NULL)

    spec <- prepSpecSimilarityPL(spec, removePrecursor = specSimParams$removePrecursor,
                                 relMinIntensity = specSimParams$relMinIntensity,
                                 minPeaks = specSimParams$minPeaks)
    if (nrow(spec) == 0) return(NULL)

    cTab <- copy(libRecs)
    groupRI <- if (!is.null(sampleRI)) sampleRI[grp] else NULL

    hash <- makeHash(baseHash, spec, cTab, libSpecs[cTab$identifier])
    resultHashCount <<- resultHashCount + 1
    resultHashes[resultHashCount] <<- hash

    cached <- if (!is.null(cachedSet)) cachedSet[[hash]] else loadCacheData("compoundsLibrary", hash, cacheDB)
    if (!is.null(cached)) return(cached)

    cTab <- unifyLibNames(cTab)
    cTab[, InChIKey1 := getIKBlock1(InChIKey)]

    lspecs <- Map(libSpecs[cTab$identifier], cTab$ion_formula_mz, cTab$identifier, f = function(sp, pmz, lid) {
      ret <- copy(sp)
      ret[, ID := seq_len(.N)]
      ret <- assignPrecursorToMSPeakList(ret, pmz)
      prepSpecSimilarityPL(ret, removePrecursor = specSimParamsLib$removePrecursor,
                           relMinIntensity = specSimParamsLib$relMinIntensity,
                           minPeaks = specSimParamsLib$minPeaks)
    })

    lspecs <- pruneList(lspecs, checkZeroRows = TRUE)
    cTab <- cTab[identifier %in% names(lspecs)]
    if (nrow(cTab) == 0) return(NULL)

    sims <- specDistRect(list(spec), lspecs, specSimParams$method, specSimParams$shift, 0, 0,
                         specSimParams$mzWeight, specSimParams$intWeight, specSimParams$absMzDev)

    cTab[, c("score", "libMatch") := sims[1, ]]

    # 🔵 NEW: filter by RI if both are available
    if (!is.null(groupRI) && "RI" %in% colnames(cTab)) {
      cTab <- cTab[abs(RI - groupRI) <= RItol]
    }

    cTab <- cTab[numGTE(score, minSim)]
    setorderv(cTab, "score", -1)
    cTab <- unique(cTab, by = "InChIKey1")
    cTab[, explainedPeaks := sapply(lspecs[identifier], nrow)]
    cTab[, SpectrumType := "MS"]  # 🔵 force GC-EI type
    cTab[, database := "library"]

    saveCacheData("compoundsLibrary", cTab, hash, cacheDB)
    return(cTab)
  }, simplify = FALSE))

  if (is.null(cachedSet))
    saveCacheSet("compoundsLibrary", resultHashes[seq_len(resultHashCount)], setHash, cacheDB)

  compList <- pruneList(compList, checkZeroRows = TRUE)
  printf("Loaded %d compounds from %d features (%.2f%%).\n", sum(unlist(lapply(compList, nrow))),
         length(compList), if (gCount == 0) 0 else length(compList) * 100 / gCount)

  return(compounds(groupAnnotations = compList, scoreTypes = c("score", "libMatch"),
                   scoreRanges = sapply(compList, function(ct) list(score = range(ct$score),
                                                                   libMatch = range(ct$libMatch)), simplify = FALSE),
                   algorithm = "library"))
})
