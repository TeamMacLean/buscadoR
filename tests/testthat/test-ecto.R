test_that("ecto domain blast changes results for different inputs", {
  a <- get_ecto(system.file("tests", "testthat", "s1.fa", package="buscadoR"))
  expect_output(str(a), "data.frame':	177 obs. of  12 variables:")
  b <- get_ecto(system.file("tests", "testthat", "s2.fa", package="buscadoR"))
  expect_output(str(b), "data.frame':	210 obs. of  12 variables:")
  expect_false(all(dim(a) == dim(b)))
})

