
#' check and combine the search files and perform the protein annotation
#' @param deeptmhmm results file from deeptmhmm (as .gz or plain text)
#' @param hmmer results file from hmmer --domtblout (as .gz or plain text)
#' @param blast results file from blastp -outfmt 6 (as. gz or plain text)
#' @param fasta original fasta file used in searches
#' @param hmmer_eval_cutoff evalue cutoff for hmmer search
#'
#' @return busc object
#' @export
buscar <- function(deeptmhmm = NULL, hmmer=NULL, blast=NULL, fasta=NULL, hmmer_eval_cutoff = 20) {



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


 return(busc)
}



#' work out whether a sequence has kinase pfams according to criteria
#' @param sn vector of sequence names with TM and SP
#' @param b buscador object
#' @importFrom rlang .data
has_kinase_pfam <- function(sn, b) {
  seqs_w_kinase_pfam_after_tm <- dplyr::left_join(b$pfams, b$tm_signal_pep, by = c("seq_name") ) %>%
    dplyr::filter(.data$b_type == "KINASE_PFAM",
                  .data$pfam_length > 250,
                  .data$seq_from > .data$tm_start)

  sn %in% seqs_w_kinase_pfam_after_tm$seq_name

}

#' work out whether a sequence has lrr pfams according to criteria
#' @param sn vector of sequence names with TM and SP
#' @param b buscador object
#' @importFrom rlang .data
has_lrr_pfam <- function(sn, b){

  seqs_w_lrr_pfam_before_tm <- dplyr::left_join(b$pfams, b$tm_signal_pep, by = c("seq_name") ) %>%
    dplyr::filter(.data$b_type == "LRR_PFAM",
                  .data$seq_to < .data$tm_start)

  sn %in% seqs_w_lrr_pfam_before_tm$seq_name

}
#' work out whether a sequence has other pfams according to criteria
#' @param sn vector of sequence names with TM and SP
#' @param b buscador object
#' @importFrom rlang .data
has_other_pfam <- function(sn,b) {
  seqs_w_other_pfam_before_tm <- dplyr::left_join(b$pfams, b$tm_signal_pep, by = c("seq_name") ) %>%
    dplyr::filter(.data$b_type == "OTHER_PFAM",
                  .data$seq_to < .data$tm_start)
  sn %in% seqs_w_other_pfam_before_tm$seq_name
}

#' work out whether a sequence has lrr blast according to criteria
#' @param sn vector of sequence names with TM and SP
#' @param b buscador object
#' @importFrom rlang .data
has_at_lrr_blast <- function(sn, b) {
  seqs_w_blast_hit_before_tm <- dplyr::left_join(b$blasts, b$tm_signal_pep, by = c("seq_name")) %>%
    dplyr::filter(.data$b_type == "LRR_BLAST",
                  .data$seq_end < .data$tm_start)

  sn %in% seqs_w_blast_hit_before_tm$seq_name
}

#' work out whether a sequence has other blast according to criteria
#' @param sn vector of sequence names with TM and SP
#' @param b buscador object
#' @importFrom rlang .data
has_at_other_blast <- function(sn, b){
  seqs_w_blast_hit_before_tm <- dplyr::left_join(b$blasts, b$tm_signal_pep, by = c("seq_name")) %>%
    dplyr::filter(.data$b_type == "OTHER_BLAST",
                  .data$seq_end < .data$tm_start)

  sn %in% seqs_w_blast_hit_before_tm$seq_name
}
#' work out whether a sequence has unspec blast according to criteria
#' @param sn vector of sequence names with TM and SP
#' @param b buscador object
#' @importFrom rlang .data
has_at_unspec_blast <- function(sn, b){

  seqs_w_blast_hit_before_tm <- dplyr::left_join(b$blasts, b$tm_signal_pep, by = c("seq_name")) %>%
    dplyr::filter(.data$b_type == "UNSPEC_BLAST",
                  .data$seq_end < .data$tm_start)

  sn %in% seqs_w_blast_hit_before_tm$seq_name
}

#' join the classification from the def matrix to the attributes determined
#' by the blast hit
#'
#' @param b buscador pbkect
#' @importFrom rlang .data
classify_protein <- function(b) {
  dplyr::left_join(b$matrix, class_def_matrix, by = c("SP", "TM", "kinase_pfam",
                                                      "lrr_pfam", "other_pfam",
                                                      "at_lrr_blast",
                                                      "at_other_blast",
                                                      "at_unspec_blast")) %>%
    dplyr::select(.data$seq_name, .data$group, .data$evidence)

}
#' create the buscador object, do the classifications
#' @param deeptmhmm_file results file from deeptmhmm (as .gz or plain text)
#' @param hmmer_file results file from hmmer --domtblout (as .gz or plain text)
#' @param blast_file results file from blastp -outfmt 6 (as. gz or plain text)
#' @param fasta_file original fasta file used in searches
#' @param hmmer_eval_cutoff evalue cutoff for hmmer search
#'
#' @return busc object
#' @importFrom rlang .data
new_buscador <- function(deeptmhmm_file = NULL, hmmer_file=NULL, blast_file=NULL, fasta_file=NULL, hmmer_eval_cutoff = 20){


  x <- list(deeptmhmm_file = deeptmhmm_file,
            hmmer_file = hmmer_file,
            blast_file = blast_file,
            fasta_file=fasta_file,
            hmmer_eval_cutoff = hmmer_eval_cutoff)


  x[['tm_signal_pep']] <- parse_raw_deeptmhmm(x$deeptmhmm_file ) %>% process_deeptmhmm()
  x[['pfams']] <- parse_raw_hmmer(x$hmmer_file, x$hmmer_eval_cutoff)
  x[['blasts']] <- parse_raw_blast(x$blast_file)



  x[['matrix']] <- tibble::tibble(seq_name = x[['tm_signal_pep']][['seq_name']],
                                SP = TRUE,
                                TM = TRUE
                                ) %>%
    dplyr::mutate(kinase_pfam = has_kinase_pfam(.data$seq_name, x),
                  lrr_pfam = has_lrr_pfam(.data$seq_name, x),
                  other_pfam = has_other_pfam(.data$seq_name, x),
                  at_lrr_blast = has_at_lrr_blast(.data$seq_name, x),
                  at_other_blast = has_at_other_blast(.data$seq_name, x),
                  at_unspec_blast = has_at_unspec_blast(.data$seq_name, x)
                  )

  x[['classes']] <- classify_protein(x)

  if(! is.null(x[['fasta_file']])) {
    x[['aastringset']] = Biostrings::readAAStringSet(x$fasta_file)
  }



  structure(x, class = "buscador")
}

#' turn the tidy long format blast dataframe into a wider sequence per line dataframe
#'
#' @param b buscador object
#' @return dataframe
#' @importFrom rlang .data
condense_blasts <- function(b) {
  b$blasts %>%
    dplyr::filter(.data$seq_name %in% b$classes$seq_name ) %>%
    tidyr::unite("blast_coord", .data$seq_start:.data$seq_end, sep="-") %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$seq_name) %>%
    dplyr::summarise(
#      sp_cut_site = .data$cut_site,
#      .data$tm_start,
#      .data$tm_end,
      blast_hit = paste0(.data$ecto, collapse=";"),
#      pfams_acc = paste0(.data$acc, collapse=";"),
      blast_loc = paste0(.data$blast_coord, collapse=";")
    ) %>%
    dplyr::distinct( )%>%
    dplyr::ungroup()

}


#' turn the tidy long format pfams dataframe into a wider sequence per line dataframe
#'
#' @param b buscador objects
#' @return dataframe
#' @importFrom rlang .data
condense_pfams <- function(b) {
  b$pfams %>%
    dplyr::filter(.data$seq_name %in% b$classes$seq_name) %>%
    tidyr::unite("pfam_coord", .data$seq_from:.data$seq_to, sep="-") %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$seq_name) %>%
    dplyr::summarise(
      #      sp_cut_site = .data$cut_site,
      #      .data$tm_start,
      #      .data$tm_end,
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

  if (! "buscador" %in% class(b)){
    stop("needs a buscador search object")
  }

  dplyr::group_by(b$classes, .data$group ) %>%
    dplyr::tally() %>%
    #dplyr::ungroup() %>%
    #dplyr::full_join( blank, by="b_type") %>%
    knitr::kable(caption = "BuscadoR RLK Finding results")

}


#' converts `buscador` search object to a tibble
#' @param b `buscador` object to coerce
#' @export
as_tibble <- function(b){

  pfam_df <- condense_pfams(b)

  blast_df <- condense_blasts(b)

  dplyr::left_join(b$classes, b$tm_signal_pep, by=c("seq_name")) %>%
   dplyr::left_join( pfam_df, by=c("seq_name") ) %>%
    dplyr::left_join(blast_df, by=c("seq_name"))

}

#' write the annotated receptor sequences to a fasta file
#'
#' @param b `buscador` search object returned from `buscar()`
#' @param out_file_path fasta file path and name to write
#' @export
#' @importFrom rlang .data
write_seqs <- function(b, out_file_path) {


  seq_df <- seq_to_df(b)
  seq_df <- dplyr::left_join(b$classes, seq_df, by=c("seq_name" = "seq_id") ) %>%
    dplyr::mutate(seq_id = paste0(.data$seq_name, "|", .data$group, "|", .data$evidence))
  seqv <- seq_df$sequence
  names(seqv) <- seq_df$seq_id
  sset <- Biostrings::AAStringSet(seqv)
  Biostrings::writeXStringSet(sset, out_file_path)


}

#' convert the aastrset to a dataframe
#' @param b `buscador` search object returned from `buscar()`
seq_to_df <- function(b){

  ids <- stringr::str_split(names(b$aastringset), " ",simplify = TRUE)[,1]
  tibble::tibble(
    seq_id = ids,
    sequence =  as.character(b$aastringset, use.names=FALSE),
    seq_length = Biostrings::width(b$aastringset)
  )
}

#' get a dataframe of receptor type data in `drawProteins` format
#' for pretty drawing. From a `buscador` search object from `buscar()`
#'
#' @param b buscador search object returned from `buscar`
#' @param which the receptor type to return, one of
#' @return dataframe
#' @export
#' @importFrom rlang .data
as_drawProteins <- function(b, which="LRR-RP") {

  keeps <- dplyr::filter(b$classes, .data$group == which)$seq_name
  seq_df <- seq_to_df(b) %>%
    dplyr::filter(.data$seq_id %in% keeps) %>%
    dplyr::rename("entryName" = "seq_id", "length"="seq_length" ) %>%
    dplyr::select(.data$entryName, .data$length) %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="CHAIN", begin=1, end=.data$length, description=which) %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)

  orders <- 1:length(seq_df$entryName)
  names(orders) <- seq_df$entryName

  tm_df <- dplyr::filter(b$tm_signal_pep, b$tm_signal_pep$seq_name %in% keeps) %>%
    dplyr::rename("entryName"="seq_name", "begin"="tm_start", "end"="tm_end") %>%
    dplyr::distinct() %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="DOMAIN", length=.data$end-.data$begin, description="transmembrane domain") %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)



  domain_df <-dplyr::filter(b$pfams, b$pfams$seq_name %in% keeps) %>%
    dplyr::rename("entryName"="seq_name", "begin"="seq_from", "end"="seq_to", "description"="hit") %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="DOMAIN", length=.data$end-.data$begin)  %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)

  blast_df <- dplyr::filter(b$blasts, b$blasts$seq_name %in% keeps) %>%
    dplyr::rename("entryName"="seq_name", "begin"="seq_start", "end"="seq_end", "description"="ecto") %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="DOMAIN", length=.data$end-.data$begin)  %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)

  return(dplyr::bind_rows(seq_df, tm_df, domain_df, blast_df) %>%
           dplyr::mutate(order = orders[.data$entryName])
  )

}



#' draw each found proteins of a given receptor type.
#'
#' @param b `buscador` search object returned from `buscar()`
#' @param which the receptor type to return, one of
#' @param label_domains write a label on the domain
#' @return ggplot2
#' @export
dibujar <- function(b, which="LRR-RP", label_domains=FALSE) {
  pd <- as_drawProteins(b, which=which)
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

