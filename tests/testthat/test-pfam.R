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
test_that("PFAM API search function returns correct result", {
  skip_on_cran()
  skip_on_ci()
  check_internet()
  check_pfam_up()
  p <- get_pfam(system.file("tests", "testthat", "pfam_ex1.fa", package="buscadoR"), email="dan.maclean@tsl.ac.uk", progress=TRUE,wait = 5, maxchecktime=360)
  expect_output(str(p), "'data.frame':	6 obs. of  9 variables:")
  pr <- process_pfam(p,1e-6)
  expect_output(str(pr), "'data.frame':	6 obs. of  12 variables:")
  expect_equal(names(pr), c("seq_name","hit","acc","eval","type","seq_from", "seq_to","hit_from", "hit_to","base_acc", "b_type","pfam_length"))
})
