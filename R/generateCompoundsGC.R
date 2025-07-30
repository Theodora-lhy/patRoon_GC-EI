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
setGeneric("generateCompoundsGC", function(fGroups, MSPeakLists, MSLibrary, minSim = 0.75,
                                           minAnnSim = minSim, absMzDev = 0.002, adduct = NULL,
                                           checkIons = "adduct", specSimParams = getDefSpecSimParams(),
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
  aapply(checkmate::assertNumber, . ~ minSim + minAnnSim + absMzDev, lower = 0, finite = TRUE, fixed = list(add = ac))
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

  # 🟠 Relax filtering for GC-EI libraries
  libRecs <- libRecs[!is.na(InChIKey) & !is.na(formula)]

  if (!isEIlib) {
    if (checkIons == "adduct")
      libRecs <- libRecs[!is.na(Precursor_type)]
    else if (checkIons == "polarity")
      libRecs <- libRecs[!is.na(Ion_mode)]
  }

  if (!is.null(adduct) && !isEIlib)
    libRecs <- libRecs[Precursor_type == as.character(adduct)]

  if (!is.null(checkIons) && checkIons == "none")
    adduct <- NULL

  # 🔵 Retention Index calculation from alkane dictionary (if provided)
  if (!is.null(RIalkaneFile)) {
    colnames(RIalkaneFile) <- tolower(colnames(RIalkaneFile))
    if (!all(c("num", "rt.min") %in% colnames(RIalkaneFile)))
      stop("Alkane RI file must contain columns 'Num' and 'RT.min' (case insensitive).")

    alkaneDict <- data.table::data.table(
      Num = as.numeric(RIalkaneFile[["num"]]),
      RT = as.numeric(RIalkaneFile[["rt.min"]])
    )
    alkaneDict <- alkaneDict[order(RT)]

    if (!"retentionTime" %in% colnames(gInfo)) {
      warning("No 'retentionTime' in groupInfo; skipping RI calculation.")
      fRI <- rep(NA_real_, nrow(gInfo))
    } else {
      fRI <- sapply(gInfo$retentionTime, function(rt) {
        lower <- which.max(alkaneDict$RT[alkaneDict$RT <= rt])
        upper <- which.min(alkaneDict$RT[alkaneDict$RT >= rt])
        if (length(lower) == 0 || length(upper) == 0) return(NA_real_)
        if (lower == upper) return(alkaneDict$Num[lower] * 100)
        rt_low <- alkaneDict$RT[lower]
        rt_high <- alkaneDict$RT[upper]
        n_low <- alkaneDict$Num[lower]
        n_high <- alkaneDict$Num[upper]
        RI <- 100 * (n_low + (rt - rt_low) / (rt_high - rt_low) * (n_high - n_low))
        return(RI)
      })
    }
    gInfo[, calcRI := fRI]
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

    if (is.null(MSPeakLists[[grp]]) || is.null(MSPeakLists[[grp]][["MS"]]))
      return(NULL)
    spec <- MSPeakLists[[grp]][["MS"]]
    if (is.null(spec) || nrow(spec) == 0)
      return(NULL)

    featMZ <- if (!is.null(gInfo[["mz"]])) gInfo[grp, mz] else NA_real_

    spec <- prepSpecSimilarityPL(spec, removePrecursor = specSimParams$removePrecursor,
                                 relMinIntensity = specSimParams$relMinIntensity,
                                 minPeaks = specSimParams$minPeaks)
    if (nrow(spec) == 0)
      return(NULL)

    cTab <- copy(libRecs)

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
    if (nrow(cTab) == 0)
      return(NULL)

    sims <- specDistRect(list(spec), lspecs, specSimParams$method, specSimParams$shift, 0, 0,
                         specSimParams$mzWeight, specSimParams$intWeight, specSimParams$absMzDev)

    cTab[, c("score", "libMatch") := sims[1, ]]
    cTab <- cTab[numGTE(score, minSim)]

    # 🔵 Optional: Retention Index matching filter
    if (!is.null(RIalkaneFile) && !is.null(gInfo$calcRI[[grp]]) && "Retention_index" %in% names(cTab)) {
      fRI <- gInfo$calcRI[[grp]]
      if (!is.na(fRI)) {
        libRI <- suppressWarnings(as.numeric(cTab$Retention_index))
        cTab <- cTab[!is.na(libRI) & abs(libRI - fRI) <= RItol]
      }
    }

    if (nrow(cTab) == 0)
      return(NULL)

    cTab[, SpectrumType := "MS"]  # Override type for GC-EI
    cTab[, explainedPeaks := sapply(lspecs[identifier], nrow)]
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
