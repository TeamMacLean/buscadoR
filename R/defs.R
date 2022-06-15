#' report the definitions used for each receptor type
#' @param which the type to show a definition for; one of
#' 'all', 'lrr_rp', 'lrr_rk', 'non_lrr_rp', 'non_lrr_rk', 'lrr_rp_ecto'
#' @return character vector
#' @export
definiciones <- function(which="all") {

  types <- c('all', 'lrr_rp', 'lrr_rk', 'non_lrr_rp', 'non_lrr_rk', 'lrr_rp_ecto')

  if (! which %in% types){
    stop('Unknown definition requested\nChoose one of:\n', paste(types, collapse=",", sep=","))
  }

  defs <- list(

    lrr_rp = "lrr_rp =
     * A signal peptide according to Phobius.
     * Exactly one Transmembrane Domain according to Phobius.
     * At least one of the `lrr_pfams` according to PFAMscan.
     * The Transmembrane domain should be closer to the C terminal than the end of the pfam hit.",

    lrr_rk = "\nlrr_rk =
     * An lrr_rp.
     * At least one `kinase_pfam` according to PFAMscan.
     * The `kinase_pfam` should be at least 250 AAs long.
     * The `kinase_pfam` closer to the C terminal than the TM domain.",

    non_lrr_rp = "\nnon_lrr_rp =
     * Not an lrr_rp or an lrr_rk.
     * A signal peptide according to Phobius.
     * Exactly one Transmembrane Domain according to Phobius.
     * At least one `non_lrr_pfam` according to PFAMscan.
     * `non_lrr_pfam` closer to the N terminal than the start of the Transmembrane domain.",

    non_lrr_rk = "\nnon_lrr_rk =
     * A non_lrr_rp.
     * At least one `kinase_pfam` according to PFAMscan.
     * The `kinase_pfam` should be at least 250 AAs long.
     * The `kinase_pfam` closer to the C terminal than the TM domain.",

    lrr_rp_ecto = "\nlrr_rp_ecto =
    * An lrr_rp or an lrr_rk
    * A BLAST hit of > 50% identity to `arabidopsis_ectodomains`
    * The BLAST hit closer to the C terminus than the signal peptide
    * The BLAST hit closer to the N terminus than the TM domain"

  )

  if (which == "all"){
    cat( paste0(defs) )
  }
  else {
    cat(defs[[which]])
  }

}
