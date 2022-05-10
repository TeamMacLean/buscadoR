test_that("buscar classifies properly", {
  expect_equal(2 * 2, 5 )
})


test_that("condense reduces data to single line properly", {
  expect_equal(3, 4)
})

test_that("as.DrawProteins returns right dataframe",{
  expect_equal(2,3)
})

test_that("do_searches returns proper result when given new protein file and restart file does not exist",{

  busc <- new_buscador(
    protein_file = system.file("tests", "testthat", "phobius_ex2.fa", package="buscadoR"),
    restart_file = file.path("test_new_protein.rds"),
    progress = TRUE,
    email = "dan.maclean@tsl.ac.uk",
    pfam_eval_cutoff = 1,
    blast_eval_cutoff = 1,
    wait = 5,
    maxchecktime = 120
  )

  busc <- do_searches(busc)
  expect_output(str(busc$phobius), "'data.frame':	1 obs. of  4 variables:")
  expect_output(str(busc$pfam), "'data.frame':	4 obs. of  12 variables:")
  expect_output(str(busc$ecto), "'data.frame':	234 obs. of  13 variables:")
  unlink("test_new_protein.rds")

})

test_that("do_searches returns proper result when given existing restart_file", {

  #Need to create a restart file first...
  busc <- new_buscador(
    protein_file = system.file("tests", "testthat", "phobius_ex2.fa", package="buscadoR"),
    restart_file = file.path("test_restart.rds"),
    progress = TRUE,
    email = "dan.maclean@tsl.ac.uk",
    pfam_eval_cutoff = 1,
    blast_eval_cutoff = 1,
    wait = 1,
    maxchecktime = 1 #short check time so job shouldn't complete on server
  )
  busc <- do_searches(busc)
  #now test that the object got populated and failed to complete
  # so is in status 'submission'

  expect_true(file.exists("test_restart.rds"))
  expect_s3_class(busc$pfam, "data.frame" )
  expect_true(attr(busc$pfam, "status") == "submission")
  expect_s3_class(busc$ecto, "data.frame")
  expect_type(busc$pfam_progress, "list")

  expect_null(busc$lrr_rp)
  expect_null(busc$lrr_rk)
  expect_null(busc$lrr_rp_rk_with_ecto)
  expect_null(busc$non_lrr_rp)
  expect_null(busc$non_lrr_rk)

  ## restart and protein file exist, test restart works
  ## send result to PFAM and check processing goes ok
  busc <- new_buscador(
    protein_file = system.file("tests", "testthat", "phobius_ex2.fa", package="buscadoR"),
    restart_file = file.path("test_restart.rds"),
    progress = TRUE,
    email = "dan.maclean@tsl.ac.uk",
    pfam_eval_cutoff = 1,
    blast_eval_cutoff = 1,
    wait = 5,
    maxchecktime = 360 ## should be long enough, but warn if it isn't
  )

  busc <- tryCatch( {do_searches(busc)},
                    error =  function(e) {
                      unlink("test_restart.rds")
                      sprintf("Failed in do_searches with restarted submission: %s", e )
                      },
                    warning = function(w) {sprint("Warning in do_searches with rstarted submission: %s", w)}
  )

  expect_output(str(busc$phobius), "'data.frame':	1 obs. of  4 variables:")
  expect_output(str(busc$pfam), "'data.frame':	4 obs. of  12 variables:")
  expect_output(str(busc$ecto), "'data.frame':	234 obs. of  13 variables:")
  expect_true(attr(busc$pfam, "status") == "complete")
  expect_null(busc$pfam_progress)
  expect_null(busc$lrr_rp)
  expect_null(busc$lrr_rk)
  expect_null(busc$lrr_rp_rk_with_ecto)
  expect_null(busc$non_lrr_rp)
  expect_null(busc$non_lrr_rk)
  unlink("test_restart.rds")
})
