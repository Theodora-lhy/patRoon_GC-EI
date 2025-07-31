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
  
  if (!is.null(RIalkaneFile)) {
    checkmate::assertDataFrame(RIalkaneFile, min.rows = 2, col.names = "named", add = ac)
    checkmate::assertSubset(c("Num", "RT.min"), colnames(RIalkaneFile), add = ac)
    RIalkaneFile <- as.data.table(RIalkaneFile)
    RIalkaneFile[, RT.sec := RT.min * 60]
  }
  
  checkmate::reportAssertions(ac)

  if (length(fGroups) == 0)
    return(compounds(algorithm = "library"))

  if (checkIons != "none")
    adduct <- checkAndToAdduct(adduct, fGroups)

  gCount <- length(fGroups)
  gInfo <- as.data.table(groupInfo(fGroups))

# Check if retention time exists under any known name
if (!"ret" %in% names(gInfo)) {
  rtCol <- grep("^ret(ention)?|^rt$", names(gInfo), value = TRUE, ignore.case = TRUE)
  if (length(rtCol) == 1) setnames(gInfo, rtCol, "ret")
  else stop("Retention time column not found in feature group info (needed to compute RI)")
}

 
  annTbl <- annotations(fGroups)
  libRecs <- records(MSLibrary)
  libSpecs <- spectra(MSLibrary)

  # 🔵 Try extracting RI from MSP metadata if present
  if ("Comments" %in% colnames(libRecs)) {
    if (!"Retention_index" %in% colnames(libRecs)) {
      libRecs[, Retention_index := as.numeric(stringr::str_extract(Comments, "(?<=Retention_index: )[0-9.]+"))]
    }
  }

  isEIlib <- all(is.na(libRecs$PrecursorMZ)) || all(is.na(libRecs$Precursor_type))

  if (!isEIlib) {
    libRecs <- libRecs[!is.na(PrecursorMZ) & !is.na(SMILES) & !is.na(InChI) & !is.na(InChIKey) & !is.na(formula)]
    if (checkIons == "adduct")
      libRecs <- libRecs[!is.na(Precursor_type)]
    else if (checkIons == "polarity")
      libRecs <- libRecs[!is.na(Ion_mode)]
  } else {
    libRecs <- libRecs[!is.na(InChIKey) & !is.na(formula)]
  }

  if (!is.null(adduct) && !isEIlib)
    libRecs <- libRecs[Precursor_type == as.character(adduct)]

  if (!is.null(checkIons) && checkIons == "none")
    adduct <- NULL

  # 🟢 Calculate RI for features if alkane file is provided
  if (!is.null(RIalkaneFile)) {
    featRTs <- fGroups[, ret]  # in seconds
    calcRI <- approx(x = RIalkaneFile$RT.sec, y = RIalkaneFile$Num * 100, xout = featRTs, rule = 2)$y
    gInfo[, RI := calcRI]
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

    spec <- prepSpecSimilarityPL(spec, removePrecursor = specSimParams$removePrecursor,
                                 relMinIntensity = specSimParams$relMinIntensity,
                                 minPeaks = specSimParams$minPeaks)
    if (nrow(spec) == 0)
      return(NULL)

    cTab <- copy(libRecs)

    # 🔵 Filter by RI if applicable
    if (!is.null(RIalkaneFile) && "RI" %in% colnames(gInfo) && "Retention_index" %in% colnames(cTab)) {
      fRI <- gInfo[group == grp, RI]
      if (!is.na(fRI))
        cTab <- cTab[!is.na(Retention_index) & abs(Retention_index - fRI) <= RItol]
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
    cTab[, SpectrumType := "MS"]  # 🟢 Explicitly label GC-EI spectra

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
