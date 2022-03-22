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
    restart_file = file.path("test_restart.rds"),
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
  unlink("test_restart.rds")

})

test_that("do_searches returns proper result when given existing restart_file", {
  #TODO
})
