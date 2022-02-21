#' Query PFAM web server.
#'
#' @param file A character file path of FASTA protein amino acid sequences a
#' @param progress Boolean, whether to show the progress bar, at default set to
#' FALSE.
#' @param email Character, a valid user email address, needed for PFAMscan use
#' @param eval PFAM hit eval cutoff
#' @return  A data frame with columns:
#' \describe{
#' \item{seq_name}{Character, name of the submitted sequence.}
#' \item{hit}{Character, the name of a PFAM hit.}
#' \item{acc}{Character, the accession of a PFAM hit.}
#' \item{eval}{Double, E-value of the hit.}
#' \item{type}{Character, PFAM type of the hit}
#' \item{seq_from}{Integer, start of match on sequence}
#' \item{seq_to}{Integer, end of match on sequence}
#' \item{hit_from}{Integer, start of match on hit}
#' \item{hit_to}{Integer, end of match on hit}
#' }
#'
#'
#'
#' @source \url{https://www.ebi.ac.uk/Tools/common/tools/help/index.html}
#'
#' @references Madeira F, Park YM, Lee J, et al. The EMBL-EBI search and sequence analysis
#' tools APIs in 2019. Nucleic Acids Research. 2019 Jul;47(W1):W636-W641.
#' \doi{10.1093/nar/gkz268}.
#' @import seqinr
#' @import httr
#' @import stringr
#' @import xml2
#' @importFrom methods is



get_pfam <- function(file, email, progress = FALSE, eval=20) {
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
  tmppfam <- tempdir()
  file_list <- split_fasta(
    path_in = file_name,
    path_out = tmppfam,
    num_seq = 100 #max pfam allows
  )
  len <- length(file_list)
  if (grepl("temp_", file_name)) {
    unlink(file_name)
  }
  if (progress) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = len,
      style = 3
    )
  }
  collected_res <- vector("list", len)
  for (i in seq_along(len)) {
    fasta_str <- readr::read_file(file_list[i])
    r <- NULL
    if (progress){
      r <- pfamscanr::pfamscan(
        fasta_str,
        email,
        evalue=eval
      ) } else {
        r <- suppressMessages(
          pfamscanr::pfamscan(
            fasta_str,
            email,
            evalue = eval
          )
        )
      }
    if (is.null(r) ){
      stop("PFAMScan failed")
    }
    collected_res[[i]] <- data.frame(
      seq_name = r$seq$name,
      hit = r$name,
      acc = r$acc,
      eval = as.numeric(r$evalue),
      type = r$type,
      seq_from = as.numeric(r$seq$from),
      seq_to = as.numeric(r$seq$to),
      hit_from = as.numeric(r$hmm$from),
      hit_to = as.numeric(r$hmm$to)
    )

    unlink(file_list[i])
    if (progress) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (progress) {
    close(pb)
  }
  collected_res <- do.call(
    rbind,
    collected_res
  )

  return(collected_res)
}

#' process raw pfam result into dataframe
#'
#' @param pr PFAM result from `get_pfam()`
#' @return dataframe
#' @importFrom rlang .data
process_pfam <- function(pr,pfam_eval) {
  pr %>%
  dplyr::filter(.data$eval < pfam_eval) %>%
    dplyr::mutate(base_acc = substr(.data$acc, 1,7 ),
                  seq_to = as.numeric(.data$seq_to),
                  seq_from = as.numeric(.data$seq_from),
                  b_type = dplyr::if_else(.data$base_acc %in% lrr_pfams, "LRR_PFAM",
                                          dplyr::if_else(.data$base_acc %in% non_lrr_pfams, "NON_LRR_PFAM",
                                                         dplyr::if_else(.data$base_acc %in% kinase_pfams, "KINASE_PFAM", "Other"))),
                  pfam_length = .data$seq_to - .data$seq_from
    ) %>%
    dplyr::distinct()

}
