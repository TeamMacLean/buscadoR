#' Run BLASTP using protein file as database and ectodomains as query
#'
#' @param file protein fasta file
#' @param progress show a progress bar
#' @export
get_ecto <- function(file, progress = FALSE) {


  if (missing(progress)) {
    progress <- FALSE
  }
  if (length(progress) > 1) {
    progress <- FALSE
    warning("progress should be of length 1, setting to default:
        progress = FALSE",
            call. = FALSE
    )
  }
  if (!is.logical(progress)) {
    progress <- as.logical(progress)
    warning("progress is not logical, converting using 'as.logical'",
            call. = FALSE
    )
  }
  if (is.na(progress)) {
    progress <- FALSE
    warning("progress was set to NA, setting to default: progress = FALSE",
            call. = FALSE
    )
  }
  if (length(file) > 1) {
    stop("one fasta file per function call can be supplied",
         call. = FALSE
    )
  }
  if (file.exists(file)) {
    file_name <- file
  } else {
    stop("cannot find file in the specified path",
         call. = FALSE
    )
  }

  if (progress) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = 5,
      style = 3
    )
  }

  db_dir <- tempdir()
  db_file <- paste0(db_dir, "/", basename(file))
  file.copy(file, db_dir)
  rBLAST::makeblastdb(db_file, dbtype="prot")
  db <- rBLAST::blast(db=db_file,type = "blastp")
  if (progress) {
    utils::setTxtProgressBar(pb, 2)
  }

  #load("data/ecto_dom_seqs.rda")
  res <- stats::predict(db, ecto_dom_seqs)
  if (progress) {
    utils::setTxtProgressBar(pb, 4)
  }



  if (progress) {
    close(pb)
  }
  return(res)
}



