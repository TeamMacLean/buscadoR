check_internet <- function(){
  if (!RCurl::url.exists("google.com")){
    skip("internet not available.")
  }
}
check_pfam_up <- function() {
  if (!RCurl::url.exists("https://www.ebi.ac.uk/Tools/services/rest/pfamscan")){
    skip("PFAMscan API not available")
  }

}
test_that("Completed PFAM API search returns correct result", {
  skip_on_cran()
  skip_on_ci()
  check_internet()
  check_pfam_up()
  aastrset <- Biostrings::readAAStringSet(system.file("tests", "testthat", "pfam_ex1.fa", package="buscadoR"))
  pfam <- submit_pfam(tostrvec(aastrset), email="dan.maclean@tsl.ac.uk", progress=TRUE,wait = 5, maxchecktime=360)
  expect_true(attr(pfam, "status") == "complete" )
  expect_output(str(pfam), "'data.frame':	6 obs. of  9 variables:")

  pr <- process_pfam(pfam,1e-6)
  expect_output(str(pr), "'data.frame':	6 obs. of  12 variables:")
  expect_equal(names(pr), c("seq_name","hit","acc","eval","type","seq_from", "seq_to","hit_from", "hit_to","base_acc", "b_type","pfam_length"))

})

test_that("Incomplete PFAM API search returns correct result", {
  skip_on_cran()
  skip_on_ci()
  check_internet()
  check_pfam_up()
  aastrset <- Biostrings::readAAStringSet(system.file("tests", "testthat", "pfam_ex1.fa", package="buscadoR"))
  ## The very short wait times mean that the search at PFAM shouldn't have time to complete, but it might...
  pfam <- submit_pfam(tostrvec(aastrset), email="dan.maclean@tsl.ac.uk", progress=TRUE,wait = 1, maxchecktime=1)
  expect_true(attr(pfam, "status") == "submission" )
  expect_output(str(pfam), "'data.frame':	1 obs. of  4 variables:")
  expect_equal(names(pfam), c("id","done","retrieved","time_posted"))

})

test_that("Restarted PFAM API search returns correct result", {
  skip_on_cran()
  skip_on_ci()
  check_internet()
  check_pfam_up()
  aastrset <- Biostrings::readAAStringSet(system.file("tests", "testthat", "pfam_ex1.fa", package="buscadoR"))


  ## The very short wait times mean that the search at PFAM shouldn't have time to complete, but it might...
  pfam <- submit_pfam(tostrvec(aastrset), email="dan.maclean@tsl.ac.uk", progress=TRUE,wait = 1, maxchecktime=1)
  busc <- list(
    pfam = pfam,
    pfam_progress = vector(mode="list", length = length(pfam$id))
  )
  names(busc$pfam_progress) <- pfam$id

  message("waiting 10 seconds for search to complete on PFAM server")
  time_start <- lubridate::now()
  while (lubridate::now() - time_start < 10){} #wait for 10 seconds

  busc <- retrieve_pfam(busc)
  expect_true(attr(busc$pfam, "status") == "complete" )
})
