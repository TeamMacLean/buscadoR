test_that("fasta splits properly", {
  fnames <-  split_fasta("multi.fa", ".", 2)
  expect_equal(length(fnames), 3)
  unlink(fnames)
})
