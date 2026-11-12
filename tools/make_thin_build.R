#!/usr/bin/env Rscript
# Produce the Bioconductor submission tarball ("thin" build).
#
# The GitHub repository stays FAT: it bundles every data layer so
# BiocManager::install("teindor/cpgdirection") works offline at run time.
# Bioconductor caps files at 5 MB and tarballs at ~10 MB, so the submission
# tarball instead retrieves the eight large layers from the cpgdirectionData
# ExperimentHub package. Both builds run the same code -- R/data_backend.R
# resolves each layer from bundled extdata when present and from
# cpgdirectionData otherwise -- so this script only:
#   1. copies the package to a temporary directory,
#   2. deletes the extdata files that cpgdirectionData serves,
#   3. moves cpgdirectionData from Suggests to Imports,
#   4. runs R CMD build.
#
# Usage: Rscript tools/make_thin_build.R  (from the package root)

hub_served <- c(
  "lookup_blood_hg19.csv.gz",
  "lookup_nasal_epithelium_hg19.csv.gz",
  "lookup_solid_tissue_hg19.csv.gz",
  "smr_directions.csv.gz",
  "cpg_positions_hg19.csv.gz",
  "brain_directions.csv.gz",
  "saliva_bridge_scores.csv.gz",
  "epicv2_probe_gene_annotation.csv.gz",
  "cpg_coordinates.csv.gz")

stopifnot(file.exists("DESCRIPTION"))
tmp <- file.path(tempdir(), "thin_build")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)
src <- normalizePath(".")
stopifnot(file.copy(src, tmp, recursive = TRUE))
pkgdir <- file.path(tmp, basename(src))

removed <- file.path(pkgdir, "inst", "extdata", hub_served)
unlink(removed)
cat("dropped", sum(!file.exists(removed)), "Hub-served layers\n")

d <- read.dcf(file.path(pkgdir, "DESCRIPTION"))
sug <- d[1, "Suggests"]
sug <- gsub("\\s*cpgdirectionData\\s*,?", "", sug)
sug <- gsub(",\\s*$", "", trimws(sug))
d[1, "Suggests"] <- sug
d[1, "Imports"] <- paste0(trimws(d[1, "Imports"]), ",\n    cpgdirectionData")
write.dcf(d, file.path(pkgdir, "DESCRIPTION"), keep.white = names(d[1, ]))

res <- system2("R", c("CMD", "build", "--no-build-vignettes", pkgdir))
stopifnot(res == 0)
tb <- list.files(pattern = "^cpgdirection_.*\\.tar\\.gz$")
tb <- tb[which.max(file.mtime(tb))]
cat(sprintf("thin tarball: %s (%.1f MB)\n", tb, file.size(tb) / 1e6))
if (file.size(tb) > 10e6) warning("tarball exceeds 10 MB")
