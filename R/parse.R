#' process raw hmmer result into dataframe
#'
#' @param hmmer_raw complete hmmer result
#' @return dataframe
#' @importFrom rlang .data
parse_raw_hmmer <- function(hmmer_raw, pfam_eval) {
  rhmmer::read_domtblout(hmmer_raw) %>%
    dplyr::select(.data$query_name, .data$domain_name, .data$domain_accession, .data$domain_ievalue, .data$ali_from, .data$ali_to, .data$hmm_from, .data$hmm_to) %>%
    dplyr::rename(seq_name = .data$query_name,
                  hit = .data$domain_name,
                  acc = .data$domain_accession,
                  eval = .data$domain_ievalue,
                  seq_from = .data$ali_from,
                  seq_to = .data$ali_to,
                  hit_from = .data$hmm_from,
                  hit_to = .data$hmm_to
                  ) %>%
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


#' process raw deeptmhmm result into basic dataframe
#'
#' @param dtm complete deeptmhmm result
#' @return dataframe
#' @importFrom rlang .data
parse_raw_deeptmhmm <- function(dtm_raw) {
  readLines(dtm_raw) %>%
    stringr::str_replace("//", "#") %>%
    readr::read_table( comment="#", col_names = c("seq_name", "type", "start", "end", "nowt")  ) %>%
    dplyr::select(-.data$nowt)
}



#' get vector of seq_names with 1 tm domain and 1 signal peptide
#'
#' @param dtm parsed deeptmhmm result from `process_raw_deeptmhmm`
#' @return character vector
#' @importFrom rlang .data
find_tm_proteins <- function(dtm) {
  r <- dtm %>%
  dplyr::group_by(.data$seq_name) %>%
    tidyr::nest() %>%
    dplyr::transmute( .data$seq_name,
      tm_domains = purrr::map_dbl(data, function(x){sum(x$type == "TMhelix")}),
      signal_peptides = purrr::map_dbl(data, function(x){sum(x$type == "signal") } )
      ) %>%
    dplyr::filter(.data$tm_domains == 1, .data$signal_peptides == 1)

  if (nrow(r) == 0) warning("No proteins with one signal peptide and one TM helix found")

  return(r$seq_name)

}

#' get filtered deeptmhmm hits
#' @param dtm parsed deeptmhmm result from `process_raw_deeptmhmm`
#' @return character vector
#' @importFrom rlang .data
process_deeptmhmm <- function(dtm) {

  tm_proteins <- find_tm_proteins(dtm)

  if (length(tm_proteins) == 0){
    return(
      tibble::tibble(
        seq_name = character(),
        cut_site = numeric(),
        tm_start = numeric(),
        tm_end = numeric()
      )
    )
  }

  dtm %>% dplyr::filter(seq_name %in% tm_proteins) %>%
    tidyr::unite(coord, .data$start, .data$end ) %>%
    tidyr::pivot_wider(names_from = .data$type, values_from = .data$coord) %>%
    dplyr::select(.data$seq_name, .data$signal, .data$TMhelix) %>%
    tidyr::separate(.data$signal, into = c("sig_start", "cut_site"), convert = TRUE) %>%
    tidyr::separate(.data$TMhelix, into = c("tm_start","tm_end"), convert = TRUE) %>%
    dplyr::select(.data$seq_name, .data$cut_site, .data$tm_start, .data$tm_end)
}

#' parse blast result file
#' @param ecto_raw blast result file (`-outfmt 7`)
#' @return dataframe
parse_raw_ecto <- function(ecto_raw) {

  readr::read_tsv(ecto_raw,
                  comment="#",
                  col_names = c("ecto", "seq_name", "percent_id",
                                "alignment_length", "mismatches", "gap_opens",
                                "ecto_start", "ecto_end", "seq_start",
                                "seq_end", "evalue", "bit_score"),
                  show_col_types = FALSE)
}

is_empty <- function(fname) {
  p <- readLines(fname) %>%
    stringr::str_replace("//", "#")
  p <- suppressWarnings( readr::read_table(p, comment = "#", col_types = readr::cols() ) )
  if (nrow(p) == 0){
    return(TRUE)
  }
  return(FALSE)

}
