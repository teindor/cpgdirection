# The data backend (2.6.0): every loader resolves its layer through
# .cpgd_data(), which tries (1) options(cpgdirection.data_dir=), then (2) the
# bundled extdata, then (3) the cpgdirectionData ExperimentHub package. The
# fat and thin builds therefore run identical code.

test_that("the local-directory override wins over the bundled copy", {
  d <- tempfile("cpgd_data_dir")
  dir.create(d)
  fake <- data.table::data.table(cpg_id = "cg_fake", chr = "chrZ", pos = 42L)
  data.table::fwrite(fake, file.path(d, "cpg_positions_hg19.csv.gz"))
  old <- options(cpgdirection.data_dir = d)
  on.exit(options(old), add = TRUE)
  got <- cpgdirection:::.cpgd_data("cpg_positions_hg19")
  expect_equal(got$pos, 42L)   # the override, not the bundled table
})

test_that("rds files are accepted in the override directory", {
  d <- tempfile("cpgd_data_dir")
  dir.create(d)
  saveRDS(data.frame(cpg_id = "cg_fake2", x = 7L),
          file.path(d, "smr_directions.rds"))
  old <- options(cpgdirection.data_dir = d)
  on.exit(options(old), add = TRUE)
  got <- cpgdirection:::.cpgd_data("smr_directions")
  expect_s3_class(got, "data.table")
  expect_equal(got$x, 7L)
})

test_that("a missing resource reports every retrieval route", {
  expect_error(cpgdirection:::.cpgd_data("no_such_layer"),
               "cpgdirection.data_dir")
  expect_false(cpgd_has_data("no_such_layer"))
})

test_that("required = FALSE degrades to NULL instead of stopping", {
  expect_null(cpgdirection:::.cpgd_data("no_such_layer", required = FALSE))
})

test_that("the 1.1.0 SMR layer is read under its Hub name, never the older table", {
  # cpgdirectionData serves two SMR tables: smr_directions (read by 0.99.x and
  # 1.0.x) and smr_directions_v2 (the 1.1.0 layer). A directory holding both,
  # like the Hub upload folder, must give the second.
  expect_equal(cpgdirection:::.cpgd_hub_name("smr_directions"), "smr_directions_v2")
  expect_equal(cpgdirection:::.cpgd_hub_name("cpg_coordinates"), "cpg_coordinates")
  d <- tempfile("cpgd_data_dir")
  dir.create(d)
  saveRDS(data.frame(cpg_id = "cg_old", direction = 1L),
          file.path(d, "smr_directions.rds"))
  saveRDS(data.frame(cpg_id = "cg_new", direction = 1L, smr_accuracy = 0.9),
          file.path(d, "smr_directions_v2.rds"))
  old <- options(cpgdirection.data_dir = d)
  on.exit(options(old), add = TRUE)
  got <- cpgdirection:::.cpgd_data("smr_directions")
  expect_equal(got$cpg_id, "cg_new")
  expect_true("smr_accuracy" %in% names(got))
})

test_that("the Hub route asks cpgdirectionData for smr_directions_v2", {
  skip_if_not_installed("cpgdirectionData")
  skip_if_not("smr_directions_v2" %in% cpgdirectionData::cpgdData_resources()$name,
              "cpgdirectionData older than 0.99.8")
  skip_if(nzchar(system.file("extdata", "smr_directions.csv.gz",
                             package = "cpgdirection")),
          "a bundled SMR table takes precedence over the Hub route")
  d <- tempfile("cpgdData_local")
  dir.create(d)
  saveRDS(data.frame(cpg_id = "cg_old", direction = 1L),
          file.path(d, "smr_directions.rds"))
  saveRDS(data.frame(cpg_id = "cg_new", direction = 1L, smr_accuracy = 0.9),
          file.path(d, "smr_directions_v2.rds"))
  old <- options(cpgdirection.data_dir = NULL, cpgdirectionData.local = d)
  on.exit(options(old), add = TRUE)
  got <- cpgdirection:::.cpgd_data("smr_directions")
  expect_s3_class(got, "data.table")
  expect_equal(got$cpg_id, "cg_new")
})
