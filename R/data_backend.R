# Where the data lives (2.6.0).
#
# Bioconductor caps a software package at 5 MB per file and ~10 MB per
# tarball; eight of this package's layers are larger than that. They therefore
# live in TWO possible homes, resolved per resource in this order:
#
#   1. options(cpgdirection.data_dir = "<dir>") -- an explicit local directory
#      holding the files (csv.gz/csv/rds). For air-gapped HPC nodes, pinned
#      analyses, and tests.
#   2. the package's own inst/extdata -- the "fat" build installed from
#      GitHub (BiocManager::install("teindor/cpgdirection")), which bundles
#      everything and needs no network at run time.
#   3. the cpgdirectionData ExperimentHub data package -- the "thin"
#      Bioconductor build, which ships no large file and fetches each layer
#      from the Hub on first use (cached locally by ExperimentHub thereafter).
#
# Every loader in the package goes through .cpgd_data(), so the two builds
# run the same code and return identical tables; only the transport differs.

# resources that the thin build serves from cpgdirectionData
.CPGD_HUB_RESOURCES <- c(
  "lookup_blood_hg19",
  "lookup_nasal_epithelium_hg19",
  "lookup_solid_tissue_hg19",
  "smr_directions",
  "cpg_positions_hg19",
  "brain_directions",
  "saliva_bridge_scores",
  "epicv2_probe_gene_annotation")

# fread() on a .gz path delegates to the R.utils package, which this package
# does not (and should not) depend on: the Bioconductor builders do not have
# it, and cpgdirection 2.99.4 failed R CMD check on every platform for exactly
# that reason. Decompress with base R into a temporary file first, then let
# fread parse the plain text. Same result on every machine, no extra
# dependency.
.cpgd_fread <- function(path, ...) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)
    inn <- gzfile(path, "rb")
    out <- file(tmp, "wb")
    tryCatch({
      while (length(chunk <- readBin(inn, "raw", 8L * 1024L * 1024L)))
        writeBin(chunk, out)
    }, finally = { close(inn); close(out) })
    path <- tmp
  }
  data.table::fread(path, showProgress = FALSE, ...)
}

.cpgd_read_data_file <- function(path) {
  if (grepl("\\.rds$", path)) {
    data.table::as.data.table(readRDS(path))
  } else {
    .cpgd_fread(path)
  }
}

# `base` is the resource name without extension, e.g. "lookup_blood_hg19".
# Returns a data.table, or NULL when required = FALSE and nothing is found.
.cpgd_data <- function(base, required = TRUE) {
  # 1. explicit local directory
  dir <- getOption("cpgdirection.data_dir")
  if (!is.null(dir) && length(dir) == 1L && nzchar(dir)) {
    for (f in file.path(dir, paste0(base, c(".csv.gz", ".csv", ".rds")))) {
      if (file.exists(f)) return(.cpgd_read_data_file(f))
    }
  }
  # 2. bundled extdata (fat build; always the case for the small resources)
  for (ext in c(".csv.gz", ".csv")) {
    p <- system.file("extdata", paste0(base, ext), package = "cpgdirection")
    if (nzchar(p) && file.exists(p)) return(.cpgd_read_data_file(p))
  }
  # 3. the Hub-backed data package (thin build)
  if (base %in% .CPGD_HUB_RESOURCES &&
      requireNamespace("cpgdirectionData", quietly = TRUE)) {
    out <- tryCatch(cpgdirectionData::cpgdData(base),
                    error = function(e) {
                      if (isTRUE(required)) {
                        stop("could not retrieve '", base, "' from ",
                             "cpgdirectionData: ", conditionMessage(e),
                             call. = FALSE)
                      }
                      NULL
                    })
    if (!is.null(out)) return(data.table::as.data.table(out))
    return(NULL)
  }
  if (!isTRUE(required)) return(NULL)
  stop("The data layer '", base, "' is not available. This installation ",
       "bundles no copy, and the cpgdirectionData package is not installed. ",
       "Either reinstall the full package ",
       "(BiocManager::install(\"teindor/cpgdirection\")), install ",
       "cpgdirectionData, or point options(cpgdirection.data_dir=) at a ",
       "directory holding the data files.", call. = FALSE)
}

# TRUE when the resource can be obtained from any tier; used by tests
.cpgd_has_data <- function(base) {
  !is.null(tryCatch(.cpgd_data(base, required = FALSE), error = function(e) NULL))
}
