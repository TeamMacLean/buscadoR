check_internet <- function(){
  if (!RCurl::url.exists("google.com")){
    skip("internet not available.")
  }
}
check_phobius_up <- function() {
  if (!RCurl::url.exists("https://phobius.sbc.su.se/")){
    skip("Phobius website not available")
  }

}


test_that("Phobius web search function returns correct result", {
  skip_on_cran()
  skip_on_ci()
  check_internet()
  check_phobius_up()

  ## 0 TM
  p <- get_phobius(system.file("tests", "testthat", "phobius_ex1.fa", package="buscadoR"), progress=TRUE)
  expect_output(str(p), "'data.frame':	1 obs. of  6 variables:")
  pr <- process_phobius(p)
  expect_output(str(pr), "'data.frame':	0 obs. of  4 variables:")
  expect_equal(names(pr), c("Name", "cut_site", "tm_start", "tm_end"))

  ##exactly 1 TM
  p <- get_phobius(system.file("tests", "testthat", "phobius_ex2.fa", package="buscadoR"), progress=TRUE)
  expect_output(str(p), "'data.frame':	1 obs. of  6 variables:")
  pr <- process_phobius(p)
  expect_output(str(pr), "'data.frame':	1 obs. of  4 variables:")
  expect_equal(names(pr), c("Name", "cut_site", "tm_start", "tm_end"))
})

