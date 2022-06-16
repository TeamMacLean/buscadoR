
#' combine the search files and perform the RLK annotation
#' @param deeptmhmm results file from deeptmhmm (as .gz or plain text)
#' @param hmmer results file from hmmer --domtblout (as .gz or plain text)
#' @param blast results file from blastp -outfmt 6 (as. gz or plain text)
#' @param fasta original fasta file used in searches
#' @param hmmer_eval_cutoff evalue cutoff for hmmer search
#' @param progress show progress bar (not implemented)
#'
#' @return busc object
#' @export
buscar <- function(deeptmhmm = NULL, hmmer=NULL, blast=NULL, fasta=NULL, hmmer_eval_cutoff = 20,  progress=TRUE ) {



  if (is_empty(deeptmhmm) ) stop("deeptmhmm results file contains no records.")
  if (!any(grepl("signal", readLines(deeptmhmm)))) stop("deeptmhmm file contains no records for signal peptides")
  if (!any(grepl("TMhelix", readLines(deeptmhmm)))) stop("deeptmhmm file contains no records of TM Helices")
  if (is_empty(hmmer)) stop("hmmer results file contains no records")
  if (is_empty(blast)) stop("blast file contains no records")

  busc <- new_buscador(
      deeptmhmm_file = deeptmhmm,
      hmmer_file = hmmer,
      blast_file = blast,
      fasta_file=fasta,
      hmmer_eval_cutoff = hmmer_eval_cutoff
    )


  busc <- get_classes(busc)


 return(busc)
}


get_classes <- function(busc) {


  deep_hmmer <- dplyr::left_join(busc$deeptmhmm, busc$hmmer, by = c("seq_name")) %>%
    dplyr::distinct()


  ## get the lrr_rp class
  busc$lrr_rp <- deep_hmmer %>%
    dplyr::filter(.data$b_type == "LRR_PFAM", .data$tm_start > .data$seq_to)


  # # get the lrr_rk class
  lrr_rk <- deep_hmmer %>%
    dplyr::filter(.data$seq_name %in% busc$lrr_rp$seq_name,
                  .data$b_type == "KINASE_PFAM",
                  .data$pfam_length >= 250, .data$tm_start < .data$seq_from) %>%
    dplyr::distinct()
  rows_to_add <- dplyr::filter(busc$lrr_rp, .data$seq_name %in% lrr_rk$seq_name )
  busc$lrr_rk <- dplyr::bind_rows(rows_to_add, lrr_rk)

  # # remove the lrr_rps that became the lrr_rks from the lrr_rp object
  busc$lrr_rp <- busc$lrr_rp %>%
    dplyr::filter(!.data$seq_name %in% busc$lrr_rk$seq_name)


  # get the non_lrr_rp class
  busc$non_lrr_rp <- deep_hmmer %>%
    dplyr::filter(! .data$seq_name %in% c(busc$lrr_rp$seq_name, busc$lrr_rk),
                  .data$b_type == "NON_LRR_PFAM", .data$seq_from < .data$tm_start ) %>%
    dplyr::distinct()

  # get the non_lrr_rk class
  busc$non_lrr_rk <- deep_hmmer %>%
    dplyr::filter(.data$seq_name %in% busc$non_lrr_rp$seq_name, .data$b_type == "KINASE_PFAM",
                  .data$pfam_length >= 250, .data$tm_start < .data$seq_from) %>%
    dplyr::distinct()

  ## remove the non_lrr_rp that became non_lrr_rk
  busc$non_lrr_rp <- busc$non_lrr_rp %>%
    dplyr::filter(! .data$seq_name %in% busc$non_lrr_rk$seq_name)

  # # get the ecto domain class
  deep_hmmer_ecto <- dplyr::left_join(deep_hmmer, busc$ecto, by = c("seq_name") ) %>%
    dplyr::distinct()

   busc$lrr_rp_rk_with_ecto <- deep_hmmer_ecto %>%
     dplyr::filter( .data$seq_name %in% c(busc$lrr_rp$seq_name, busc$lrr_rk),
                    .data$percent_id > 50, .data$seq_start > .data$cut_site, .data$seq_start < .data$tm_start ) %>%
     tidyr::unite(pfam_coord, .data$seq_from:.data$seq_to, sep="-", remove=FALSE)

  return(busc)
}



#' create the buscador object
new_buscador <- function(...){


  x <- list(...)

  x['lrr_rp'] = list(NULL)
  x['lrr_rk'] = list(NULL)
  x['lrr_rp_rk_with_ecto'] = list(NULL)
  x['non_lrr_rp'] = list(NULL)
  x['non_lrr_rk'] = list(NULL)

  if(! is.null(x[['fasta_file']])) {
    x[['aastringset']] = Biostrings::readAAStringSet(x$fasta_file)
  }

  x[['hmmer']] = parse_raw_hmmer(x$hmmer_file, x$hmmer_eval_cutoff)
  x[['deeptmhmm']] = parse_raw_deeptmhmm(x$deeptmhmm_file ) %>% process_deeptmhmm()
  x[['ecto']] = parse_raw_ecto(x$blast_file)





  b <- structure(x,
                 class = "buscador"
  )
}

#' turn the tidy long format dataframe into a wider sequence per line dataframe
#'
#' @param df dataframe to condense
#' @return dataframe
#' @importFrom rlang .data
condense <- function(df) {
  df %>%
    tidyr::unite(pfam_coord, .data$seq_from:.data$seq_to, sep="-") %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$seq_name) %>%
    dplyr::summarise(
      sp_cut_site = .data$cut_site,
      .data$tm_start,
      .data$tm_end,
      pfams_hit = paste0(.data$hit, collapse=";"),
      pfams_acc = paste0(.data$acc, collapse=";"),
      pfams_loc = paste0(.data$pfam_coord, collapse=";")
    ) %>%
    dplyr::distinct( )%>%
    dplyr::ungroup()

}

#' make a table of the number of each type of receptor from an internet search
#'
#' @param b an object of class `buscador`, usually the search results from `buscar()`
#' @return knitr::kable object
#' @export
#' @importFrom rlang .data
mesa <- function(b){

  if ("buscador" %in% class(b)){
    df <- as.data.frame(b)
  }

  if (! "b_type" %in% colnames(df) ){
    stop("column b_type not found in data, cannot make table")
  }


  blank <- data.frame(
    b_type = c("lrr_rp", "lrr_rk", "non_lrr_rp", "non_lrr_rk", "lrr_rp_rk_with_ecto")
  )
  dplyr::group_by(df, .data$b_type ) %>%
    dplyr::summarise(count = dplyr::n() ) %>%
    dplyr::ungroup() %>%
    dplyr::full_join( blank, by="b_type") %>%
    knitr::kable(caption = "BuscadoR RLK Finding results")

}


#' generic for converting `buscador` search object to a dataframe
#' @param x `buscador` object to coerce
#' @param ... parameters for other functions
#' @export
#' @importFrom rlang .data
as.data.frame.buscador <- function(x,...){

  con_lrr_rp <- x$lrr_rp %>%
    condense() %>%
    dplyr::mutate(b_type = "lrr_rp")

  con_lrr_rk <- x$lrr_rk %>%
    condense() %>%
    dplyr::mutate(b_type = "lrr_rk")

  con_non_lrr_rp <- x$non_lrr_rp %>%
    condense() %>%
    #dplyr::rename(non_lrr_rp_pfams_hit = pfams_hit, non_lrr_rp_pfams_acc = pfams_acc, non_lrr_rp_pfams_loc = pfams_loc) %>%
    dplyr::mutate(b_type = "non_lrr_rp")

  con_non_lrr_rk <- x$non_lrr_rk %>%
    condense() %>%
    #dplyr::rename(non_lrr_rp_with_rk_pfams_hit = pfams_hit, non_lrr_rp__with_rk_pfams_acc = pfams_acc, non_lrr_rp__with_pfams_loc = pfams_loc) %>%
    dplyr::mutate(b_type = "non_lrr_rk")

  con_lrr_rp_rk_ecto <- x$lrr_rp_rk_with_ecto %>%
    tidyr::unite(hit_coord, .data$hit_from, .data$hit_to, sep="-" ) %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$seq_name) %>%
    dplyr::summarise(
      sp_cut_site = .data$cut_site,
      .data$tm_start,
      .data$tm_end,
      pfams_hit = paste0(.data$hit, collapse=";"),
      pfams_acc = paste0(.data$acc, collapse=";"),
      pfams_loc = paste0(.data$pfam_coord, collapse=";"),
      ectos_hit = paste0(.data$ecto, collapse=";"),
      ectos_coord = paste0(.data$hit_coord, collapse=";")

    ) %>%
    dplyr::distinct( ) %>%
    dplyr::ungroup()  %>%
    dplyr::mutate(b_type = "lrr_rp_rk_with_ecto")

  dplyr::bind_rows(
    list(
      con_lrr_rp, con_lrr_rk, con_non_lrr_rp, con_non_lrr_rk, con_lrr_rp_rk_ecto
    )
  ) %>%
    dplyr::distinct()
}

#' write the annotated receptor sequences to a fasta file
#'
#' @param b `buscador` search object returned from `buscar()`
#' @param out_file_path fasta file path and name to write
#' @export
#' @importFrom rlang .data
write_seqs <- function(b, out_file_path) {

  ninfo <- data.frame(
    seq_id = as.data.frame(b)$seq_name,
    b_type = as.data.frame(b)$b_type
  )
  seq_df <- seq_to_df(b)
  seq_df <- dplyr::left_join(ninfo, seq_df, by="seq_id") %>%
    dplyr::mutate(seq_id = paste0(.data$seq_id, "|", .data$b_type))
  seqv <- seq_df$sequence
  names(seqv) <- seq_df$seq_id
  sset <- Biostrings::AAStringSet(seqv)
  Biostrings::writeXStringSet(sset, out_file_path)


}

#' convert the aastrset to a dataframe
#' @param b `buscador` search object returned from `buscar()`
seq_to_df <- function(b){

  ids <- stringr::str_split(names(b$aastringset), " ",simplify = TRUE)[,1]
  data.frame(
    seq_id = ids,
    sequence =  as.character(b$aastringset, use.names=FALSE),
    seq_length = Biostrings::width(b$aastringset)
  )
}

#' get a dataframe of receptor type data in `drawProteins` format
#' for pretty drawing. From a `buscador` search object from `buscar()`
#'
#' @param b buscador search object returned from `buscar`
#' @param which the receptor type to return, one of 'lrr_rp', 'lrr_rk', 'non_lrr_rp', 'non_lrr_rk', 'lrr_rp_ecto'
#' @return dataframe
#' @export
#' @importFrom rlang .data
as.drawProteins <- function(b, which="lrr_rp") {

  seq_df <- seq_to_df(b) %>%
    dplyr::filter(.data$seq_id %in% b[[which]]$seq_name) %>%
    dplyr::rename("entryName" = "seq_id", "length"="seq_length" ) %>%
    dplyr::select(.data$entryName, .data$length) %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="CHAIN", begin=1, end=.data$length, description=which) %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)

  orders <- 1:length(seq_df$entryName)
  names(orders) <- seq_df$entryName

  tm_df <- b[[which]] %>%
    dplyr::rename("entryName"="seq_name", "begin"="tm_start", "end"="tm_end") %>%
    dplyr::distinct() %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="DOMAIN", length=.data$end-.data$begin, description="transmembrane domain") %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)

  domain_df <- b[[which]] %>%
    dplyr::rename("entryName"="seq_name", "begin"="seq_from", "end"="seq_to", "description"="hit") %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="DOMAIN", length=.data$end-.data$begin)  %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)

  return(dplyr::bind_rows(seq_df, tm_df, domain_df) %>%
           dplyr::mutate(order = orders[.data$entryName])
  )

}



#' draw each found proteins of a given receptor type.
#'
#' @param b `buscador` search object returned from `buscar()`
#' @param which the receptor type to return, one of 'lrr_rp', 'lrr_rk', 'non_lrr_rp', 'non_lrr_rk', 'lrr_rp_ecto'
#' @param label_domains write a label on the domain
#' @return ggplot2
#' @export
dibujar <- function(b, which="lrr_rp", label_domains=FALSE) {
  pd <- as.drawProteins(b, which=which)
  drawProteins::draw_canvas(pd) %>%
    drawProteins::draw_chains(pd) %>%
    drawProteins::draw_domains(pd, label_domains = label_domains) +
    ggplot2::theme_bw(base_size = 20) +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
                   panel.grid.major=ggplot2::element_blank()) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank()) +
    ggplot2::theme(panel.border = ggplot2::element_blank())
}

