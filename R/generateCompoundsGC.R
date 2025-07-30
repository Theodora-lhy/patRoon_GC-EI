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
  if (!is.null(RIalkaneFile))
    checkmate::assertDataFrame(RIalkaneFile, min.rows = 2, col.names = "named",
                               must.include = c("Num", "RT.min"), add = ac)
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

  # Try to extract RI values from Comments column if present
  if (!"Retention_index" %in% names(libRecs)) {
    if ("Comments" %in% names(libRecs)) {
      comments_vec <- as.character(libRecs$Comments)
      extracted_RI <- suppressWarnings(as.numeric(sub(".*[Rr]etention[_ ]index:? ?", "", comments_vec)))
      libRecs$Retention_index <- ifelse(is.na(extracted_RI), NA_real_, extracted_RI)
    } else {
      libRecs$Retention_index <- NA_real_
    }
  }

  # GC-EI fallback: skip precursor filtering if largely missing
  isEIlib <- all(is.na(libRecs$PrecursorMZ)) || all(is.na(libRecs$Precursor_type))

  # Relax filters for GC-EI libraries
  if (!isEIlib) {
    libRecs <- libRecs[!is.na(PrecursorMZ) & !is.na(SMILES) & !is.na(InChI) & !is.na(InChIKey) & !is.na(formula)]
    if (checkIons == "adduct")
      libRecs <- libRecs[!is.na(Precursor_type)]
    else if (checkIons == "polarity")
      libRecs <- libRecs[!is.na(Ion_mode)]
  } else {
    libRecs <- libRecs[!is.na(InChIKey) & !is.na(formula)]  # Relaxed filter
  }

  if (!is.null(adduct) && !isEIlib)
    libRecs <- libRecs[Precursor_type == as.character(adduct)]

  if (!is.null(checkIons) && checkIons == "none")
    adduct <- NULL

  # 🟢 Compute RI for each feature group (if RT and alkane table provided)
  if (!is.null(RIalkaneFile)) {
    rtVec <- gInfo$ret / 60  # Convert seconds to minutes
    fRI <- sapply(rtVec, function(rt) {
      before <- RIalkaneFile[RIalkaneFile$RT.min <= rt, ]
      after <- RIalkaneFile[RIalkaneFile$RT.min > rt, ]
      if (nrow(before) == 0 || nrow(after) == 0) return(NA_real_)
      i <- which.max(before$RT.min)
      j <- which.min(after$RT.min)
      c1 <- before$Num[i]
      c2 <- after$Num[j]
      t1 <- before$RT.min[i]
      t2 <- after$RT.min[j]
      100 * c1 + 100 * (rt - t1)/(t2 - t1) * (c2 - c1)
    })
    gInfo[, calcRI := fRI]
  }

  cacheDB <- openCacheDBScope()
  baseHash <- makeHash(minSim, minAnnSim, absMzDev, adduct, checkIons, specSimParams, specSimParamsLib)
  setHash <- makeHash(fGroups, MSPeakLists, MSLibrary, baseHash)
  cachedSet <- loadCacheSet("compoundsLibrary", setHash, cacheDB)
  resultHashes <- vector("character", gCount)
  resultHashCount <- 0

  printf("Processing %d feature groups with a library of %d records...\n", gCount, nrow(libRecs))

  compList <- withProg(gCount, FALSE, sapply(seq_len(gCount), function(i)
  {
    doProgress()
    grp <- names(fGroups)[i]
    if (is.null(MSPeakLists[[grp]]) || is.null(MSPeakLists[[grp]][["MS"]])) return(NULL)
    spec <- MSPeakLists[[grp]][["MS"]]
    if (nrow(spec) == 0) return(NULL)

    spec <- prepSpecSimilarityPL(spec, removePrecursor = specSimParams$removePrecursor,
                                 relMinIntensity = specSimParams$relMinIntensity,
                                 minPeaks = specSimParams$minPeaks)
    if (nrow(spec) == 0) return(NULL)

    cTab <- copy(libRecs)

    # 🔵 Apply RI filtering
    if (!is.null(RIalkaneFile)) {
      fRI <- gInfo$calcRI[i]
      if (!is.na(fRI)) {
        cTab <- cTab[is.na(Retention_index) | abs(Retention_index - fRI) <= RItol]
      }
    }

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
    cTab <- cTab[numGTE(score, minSim)]
    setorderv(cTab, "score", -1)
    cTab <- unique(cTab, by = "InChIKey1")
    cTab[, explainedPeaks := sapply(lspecs[identifier], nrow)]
    cTab[, database := "library"]
    cTab[, SpectrumType := "MS"]  # force spectrum type for GC-EI

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

