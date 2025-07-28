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

# GENERIC
setGeneric("generateCompoundsGC", function(fGroups, MSPeakLists, MSLibrary, minSim = 0.75,
                                             minAnnSim = minSim, absMzDev = 0.002, adduct = NULL,
                                             checkIons = "adduct", specSimParams = getDefSpecSimParams(),
                                             specSimParamsLib = getDefSpecSimParams(), RI = NULL, RItol = 10) {
  standardGeneric("generateCompoundsGC")
})

# METHOD
setMethod("generateCompoundsGC", "featureGroups", function(fGroups, MSPeakLists, MSLibrary, 
                                                           minSim = 0.75, minAnnSim = minSim, absMzDev = 0.002,
                                                           adduct = NULL, checkIons = "adduct", 
                                                           specSimParams = getDefSpecSimParams(), 
                                                           specSimParamsLib = getDefSpecSimParams(),
                                                           RI = NULL, RItol = 10) 
{
  # Set GC-EI-safe default spectrum cleaning
  specSimParams$removePrecursor <- FALSE
  specSimParams$relMinIntensity <- 0.005
  specSimParams$minPeaks <- 5

  specSimParamsLib$removePrecursor <- FALSE
  specSimParamsLib$relMinIntensity <- 0.005
  specSimParamsLib$minPeaks <- 5

  # Validations
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
  annTbl <- annotations(fGroups)
  libRecs <- records(MSLibrary)
  libSpecs <- spectra(MSLibrary)

  # 🔁 GC-EI: Relax field requirements
  libRecs <- libRecs[!is.na(InChIKey) & !is.na(formula)]  # Removed SMILES, InChI, PrecursorMZ, etc.

  # 🔁 GC-EI: Optional RI filtering
  if (!is.null(RI) && "RI" %in% names(libRecs)) {
    libRecs[, RI := as.numeric(RI)]
  }

  if (!is.null(checkIons) && checkIons == "none")
    adduct <- NULL

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

    spec <- prepSpecSimilarityPL(spec, removePrecursor = specSimParams$removePrecursor,
                                 relMinIntensity = specSimParams$relMinIntensity,
                                 minPeaks = specSimParams$minPeaks)
    if (nrow(spec) == 0)
      return(NULL)

    cTab <- copy(libRecs)

    # 🔁 Apply RI filtering
    if (!is.null(RI) && "RI" %in% names(cTab)) {
      thisRI <- RI[grp]
      if (!is.na(thisRI))
        cTab <- cTab[abs(RI - thisRI) <= RItol]
    }

    if (nrow(cTab) == 0)
      return(NULL)

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
    setorderv(cTab, "score", -1)
    cTab <- unique(cTab, by = "InChIKey1")

    cTab[, explainedPeaks := sapply(lspecs[identifier], nrow)]
    cTab[, database := "library"]
    cTab[, SpectrumType := "MS"]  # 🔵 Force MS type for GC-EI full scan

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
