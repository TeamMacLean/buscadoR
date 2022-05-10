#' Query PFAM web server.
#'
#' @param file A character vector of FASTA protein amino acid sequences
#' @param progress Boolean, whether to report progress bar, at default set to
#' FALSE.
#' @param email Character, a valid user email address, needed for PFAMscan use
#' @param eval PFAM hit eval cutoff
#' @return  The buscador$pfam object populated with submission query ids and marked with class 'submission'
#' or if search completed in time a data frame with class 'complete':
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
#' @importFrom rlang .data
submit_pfam <- function(fasta_str, email, progress = FALSE, eval=20, maxchecktime=120, wait=10) {
 ### progress doesn't work


    r <- NULL
    pfam_submission <- NULL
    time_now <- lubridate::now()
    if (progress){
      r <- pfamscanr::pfamscan(
        fasta_str,
        email,
        evalue=eval,
        maxchecktime = maxchecktime,
        wait = wait
      ) } else {
        r <- suppressMessages(
          pfamscanr::pfamscan(
            fasta_str,
            email,
            evalue = eval,
            maxchecktime = maxchecktime,
            wait = wait
          )
        )
      }
      pfam_submission <- pfamscanr::pfamscanr_queries$dat # pull out of environment object belonging to session
      attr(pfam_submission, "status") <- "submission"
    if (is.null(r) ){ # if the search didnt return in time.
      return(
        pfam_submission %>%
               dplyr::filter(.data$time_posted > time_now )
        # This is done to stop the list of submissions growing if the user has a persistent R session
        # and has run the code more than once.
        # The environment from pfamscanr holds the submissions for the entire session, not just the
        # last run of the code
        # A single submission id will represent 100 proteins (or the remainder of the set)
        )
      #stop("PFAMScan failed, null obtained from PFAM. This is most likely due to a timeout, please try increasing 'wait' and 'maxchecktime' ")
    } else {
      collected_res <- rename_pfam_cols(r)
      attr(collected_res, "status") <- "complete"
      return(collected_res)
    }
}

#' reformat pfam result columns
#'
#' @param r pfam result from pfamscanr
#'
rename_pfam_cols <- function(r) {
  data.frame(
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

}

#' check the PFAM server for an update on job progress
#'
#' looks at the buscador$pfam_progress object and checks progress of ID on PFAM.
#' If all IDs are complete and retrieved, then the individual dataframes are combined
#' and put into buscador$pfam and given class "complete". If not the
#' buscador$pfam_progress is updated as far as possible but the class remains
#' "submission"
#' @param busc buscador object
#' @return buscador object with updated pfam_progress slot or buscador object with
#' complete pfam slot
retrieve_pfam <- function(busc) {

  n <- names(busc$pfam_progress)
  busc$pfam_progress = lapply(n, function(pfam){ pfamscanr::pfamscan_retrieve(pfam)})
  names(busc$pfam_progress) <- n


  done <- sum(sapply(busc$pfam_progress, function(x) !is.null(x), simplify = TRUE))

  message(paste0("Done ", done, " of ", length(busc$pfam_progress), " PFAMScan jobs at ebi.ac.uk"))
  if (done == length(busc$pfam_progress)){
    ## make the big dataframe
    res <- lapply(busc$pfam_progress, rename_pfam_cols)
    collected_res = do.call(rbind, res)
    busc$pfam = collected_res
    attr(busc$pfam, "status") <- "complete"
  } else {
    message(paste0("PFAMScan jobs not complete at ebi.ac.uk. Please restart later."))
  }


  return(busc)
}

#' process raw pfam result into dataframe
#'
#' @param pr complete PFAM result
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
