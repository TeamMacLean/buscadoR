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


#' perform internet searches and classify given proteins into receptor types
#'
#' runs the internet search using the Phobius server, the PFAM API and local BLAST. Then collates the results
#' and classifies proteins into receptor types as defined in `definiciones()`
#'
#' @param file_path path to protein fasta file
#' @param progress show a progress bar as job runs (not recommended for non interactive sessions)
#' @param email a valid email address
#' @param pfam_eval_cutoff exclude pfam hits with eval over this value (20)
#' @param blast_eval_cutoff exclude blast hits with eval over this value (1e-6)
#' @param ... other arguments for internet search functions - not yet implemented
#' @return list of tidy dataframes of each type of receptor and its Phobius, PFAM and BLAST hits, class 'busco'
#' @export
#' @importFrom rlang .data
buscar <- function(file_path, progress=FALSE, email = NULL, pfam_eval_cutoff=20, blast_eval_cutoff=1e-6, ...){

  #load("data/pfam_data.rda")
  searches <- do_searches(file_path, progress, email,
                          pfam_eval=pfam_eval_cutoff, blast_eval=blast_eval_cutoff, ...)


  searches$file_path <- file_path


  ## add pfam info to all passing signal and tm domain carrying proteins
  phobius_pfam <- dplyr::left_join(searches$phobius, searches$pfam, by = c("Name" = "seq_name")) %>%
    dplyr::distinct()


  ## get the lrr_rp class
  searches$lrr_rp <- phobius_pfam %>%
    dplyr::filter(.data$b_type == "LRR_PFAM", .data$tm_start > .data$seq_to)


  # get the lrr_rk class
   lrr_rk <- phobius_pfam %>%
     dplyr::filter(.data$Name %in% searches$lrr_rp$Name,
                   .data$b_type == "KINASE_PFAM",
                   .data$pfam_length >= 250, .data$tm_start < .data$seq_from) %>%
     dplyr::distinct()
   rows_to_add <- dplyr::filter(searches$lrr_rp, .data$Name %in% lrr_rk$Name )
    searches$lrr_rk <- dplyr::bind_rows(rows_to_add, lrr_rk)

  # remove the lrr_rps that became the lrr_rks from the lrr_rp object
   searches$lrr_rp <- searches$lrr_rp %>%
    dplyr::filter(!.data$Name %in% searches$lrr_rk$Name)


  # get the non_lrr_rp class
   searches$non_lrr_rp <- phobius_pfam %>%
     dplyr::filter(! .data$Name %in% c(searches$lrr_rp$Name, searches$lrr_rk),
                   .data$b_type == "NON_LRR_PFAM", .data$seq_from < .data$tm_start ) %>%
     dplyr::distinct() #%>%

  # get the non_lrr_rk class
   searches$non_lrr_rk <- phobius_pfam %>%
     dplyr::filter(.data$Name %in% searches$non_lrr_rp$Name, .data$b_type == "KINASE_PFAM",
                   .data$pfam_length >= 250, .data$tm_start < .data$seq_from) %>%
     dplyr::distinct()

    ## remove the non_lrr_rp that became non_lrr_rk
   searches$non_lrr_rp <- searches$non_lrr_rp %>%
     dplyr::filter(! .data$Name %in% searches$non_lrr_rk)

  # get the ecto domain class
   phobius_pfam_ecto <- dplyr::left_join(phobius_pfam, searches$ecto, by = c("Name" = "SubjectID") ) %>%
     dplyr::distinct()

   searches$lrr_rp_rk_with_ecto <- phobius_pfam_ecto %>%
     dplyr::filter( .data$Name %in% c(searches$lrr_rp$Name, searches$lrr_rk),
                    .data$Perc.Ident > 50, .data$S.start > .data$cut_site, .data$S.start < .data$tm_start
                    ) %>%
     tidyr::unite(.data$pfam_coord, .data$seq_from:.data$seq_to, sep="-", remove=FALSE)

  class(searches) <- "busco"
  return(searches)
}

#' turn the tidy long format dataframe into a wider sequence per line dataframe
#'
#' @param df dataframe to condense
#' @return dataframe
#' @importFrom rlang .data
condense <- function(df) {
  df %>%
    tidyr::unite(.data$pfam_coord, .data$seq_from:.data$seq_to, sep="-") %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$Name) %>%
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
#' @param b an object of class `busco`, usually the search results from `buscar()`
#' @return knitr::kable object
#' @export
#' @importFrom rlang .data
mesa <- function(b){

  if ("busco" %in% class(b)){
    df <- as.data.frame(b)
  }

  if (! "b_type" %in% colnames(df) ){
    stop("column b_type not found in data, cannot make table")
  }

  dplyr::group_by(df, .data$b_type ) %>%
    dplyr::summarise(count = dplyr::n() ) %>%
    dplyr::ungroup() %>%
    knitr::kable(caption = "BuscadoR RLK Finding results")

}

#' run searches over the internet and local BLAST.
#'
#' Runs searches for proteins at Phobius, PFAM and local BLAST.
#' Reformat returned results from searches and apply basic quality filters
#'
#' @param file_path path to protein fasta file
#' @param progress show a progress bar as job runs (not recommended for non interactive sessions)
#' @param email a valid email address
#' @param pfam_eval exclude pfam hits with eval over this value (20)
#' @param blast_eval exclude blast hits with eval over this value (1e-6)
#' @return list of dataframes of search results
#' @importFrom rlang .data
do_searches <- function(file_path, progress=FALSE, email=NULL, pfam_eval=20, blast_eval=1e-6){

  result <- list(
    phobius = NULL,
    pfam = NULL,
    ecto = NULL
  )

  if (progress) message("Starting phobius search")

  result$phobius <- get_phobius(file_path, progress=progress) %>%
    dplyr::filter(.data$is.phobius == TRUE, .data$tm == 1) %>%
    dplyr::mutate(
      tm_start = as.numeric( stringr::str_split(
        stringr::str_split(.data$prediction, "[oi]", simplify = TRUE)[,2],
        "-", simplify=TRUE)[,1] ) ,
      tm_end = as.numeric( stringr::str_split(
        stringr::str_split(.data$prediction, "[oi]", simplify = TRUE)[,2],
        "-", simplify=TRUE)[,2] ),
    ) %>%
    dplyr::select(.data$Name, .data$cut_site, .data$tm_start, .data$tm_end) %>%
    dplyr::distinct()

  if (progress) message(paste("Phobius complete, found", length(result$phobius$Name), "proteins with Signal Peptide and single TM domain\n"))

  if (progress) message("Starting PFAM with found proteins")

  filtered_protein_tmpfile <- keep_phobius_hits(result$phobius, file_path)
  result$pfam <- get_pfam(filtered_protein_tmpfile, email, progress=progress, eval=pfam_eval) %>%
    dplyr::filter(.data$eval < .data$pfam_eval) %>%
    dplyr::mutate(base_acc = substr(.data$acc, 1,7 ),
                  seq_to = as.numeric(.data$seq_to),
                  seq_from = as.numeric(.data$seq_from),
                  b_type = dplyr::if_else(.data$base_acc %in% lrr_pfams, "LRR_PFAM",
                                          dplyr::if_else(.data$base_acc %in% non_lrr_pfams, "NON_LRR_PFAM",
                                                         dplyr::if_else(.data$base_acc %in% kinase_pfams, "KINASE_PFAM", "Other"))),
                  pfam_length = .data$seq_to - .data$seq_from
    ) %>%
    dplyr::distinct()

  if (progress) message("PFAM complete\nStarting BLAST for ectodomains with found proteins\n")
  result$ecto <- get_ecto(file_path, progress=progress) %>%
    dplyr::filter(.data$E < blast_eval) %>%
    tidyr::unite(.data$hit_coord, .data$S.start:.data$S.end, sep="-", remove=FALSE)

  return(result)
}

#' Extract proteins with Phobius hits from input file and return a tmpfile fasta
#' of those proteins for subsequent searches
#'
#' @param df dataframe of Phobius hits
#' @param file_path path to protein fasta file
#' @return tmpfile path
keep_phobius_hits <- function(df, file_path) {
  p <- Biostrings::readAAStringSet(file_path)
  p <- p[df$Name]
  tmpfile <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(p, tmpfile)
  return(tmpfile)
}


#' get the dataframe of PFAM results from a `busco` search object from `buscar()`
#'
#' @param b busco search object returned from `busco`
#' @return dataframe
#' @export
pfam_results <- function(b) {
  b$pfam
}

#' get the dataframe of Phobius results from a `busco` search object from `buscar()`
#'
#' @param b busco search object returned from `busco`
#' @return dataframe
#' @export
phobius_results <- function(b) {
  b$phobius
}
#' get the dataframe of Ecto domain results from a `busco` search object from `buscar()`
#'
#' @param b busco search object returned from `busco`
#' @return dataframe
#' @export
ecto_results <- function(b) {
  b$ecto
}

#' get the dataframe of lrr_rp results from a `busco` search object from `buscar()`
#'
#' @param b busco search object returned from `busco`
#' @return dataframe
#' @export
lrr_rp <- function(b) {
  b$lrr_rp
}

#' get the dataframe of lrr_rk results from a `busco` search object from `buscar()`
#'
#' @param b busco search object returned from `busco`
#' @return dataframe
#' @export
lrr_rk <- function(b) {
  b$lrr_rk
}

#' get the dataframe of non_lrr_rp results from a `busco` search object from `buscar()`
#'
#' @param b busco search object returned from `busco`
#' @return dataframe
#' @export
non_lrr_rp <- function(b) {
  b$non_lrr_rp
}

#' get the dataframe of non_lrr_rk results from a `busco` search object from `buscar()`
#'
#' @param b busco search object returned from `busco`
#' @return dataframe
#' @export
non_lrr_rk <- function(b) {
  b$non_lrr_rk
}

#' get the dataframe of lrr_rp_rk_with_ecto results from a `busco` search object from `buscar()`
#'
#' @param b busco search object returned from `busco`
#' @return dataframe
#' @export
lrr_rp_rk_with_ecto <- function(b) {
  b$lrr_rp_rk_with_ecto
}

#' get a dataframe of receptor type data in `drawProteins` format
#' for pretty drawing. From a `busco` search object from `buscar()`
#'
#' @param b busco search object returned from `busco`
#' @param which the receptor type to return, one of 'lrr_rp', 'lrr_rk', 'non_lrr_rp', 'non_lrr_rk', 'lrr_rp_ecto'
#' @return dataframe
#' @export
#' @importFrom rlang .data
as.drawProteins <- function(b, which="lrr_rp") {

  seq_df <- seq_to_df(b) %>%
    dplyr::filter(.data$seq_id %in% b[[which]]$Name) %>%
    dplyr::rename("entryName" = "seq_id", "length"="seq_length" ) %>%
    dplyr::select(.data$entryName, .data$length) %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="CHAIN", begin=1, end=.data$length, description=which) %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)

  orders <- 1:length(seq_df$entryName)
  names(orders) <- seq_df$entryName

  tm_df <- b[[which]] %>%
    dplyr::rename("entryName"="Name", "begin"="tm_start", "end"="tm_end") %>%
    dplyr::distinct() %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="DOMAIN", length=.data$end-.data$begin, description="transmembrane domain") %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)

  domain_df <- b[[which]] %>%
    dplyr::rename("entryName"="Name", "begin"="seq_from", "end"="seq_to", "description"="hit") %>%
    dplyr::mutate(accession=.data$entryName, taxid=2712, order=1, type="DOMAIN", length=.data$end-.data$begin)  %>%
    dplyr::select(.data$type, .data$description, .data$begin, .data$end, .data$length, .data$accession, .data$entryName, .data$taxid, .data$order)

  return(dplyr::bind_rows(seq_df, tm_df, domain_df) %>%
         dplyr::mutate(order = orders[.data$entryName])
         )

}



#' draw each found proteins of a given receptor type.
#'
#' @param b `busco` search object returned from `buscar()`
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

#' load the protein fasta file in to a dataframe
#' @param b `busco` search object returned from `buscar()`
seq_to_df <- function(b){

  seqs <- Biostrings::readAAStringSet(filepath = b$file_path)
  data.frame(
    seq_id = names(seqs),
    sequence =  as.character(seqs, use.names=FALSE),
    seq_length = Biostrings::width(seqs)
  )
}

#' write the receptor sequences to a fasta file
#'
#' @param b `busco` search object returned from `buscar()`
#' @param out_file_path fasta file path and name to write
#' @export
#' @importFrom rlang .data
write_seqs <- function(b, out_file_path) {

  ninfo <- data.frame(
    seq_id = as.data.frame(b)$Name,
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

#' generic for converting `busco` searh object to a dataframe
#' @param x `busco` object to coerce
#' @param ... parameters for other functions
#' @export
#' @importFrom rlang .data
as.data.frame.busco <- function(x,...){

  con_lrr_rp <- x$lrr_rp %>%
    condense() %>%
    #dplyr::rename(lrr_pfams_hit = pfams_hit, lrr_pfams_acc = pfams_acc, lrr_pfams_loc = pfams_loc) %>%
    dplyr::mutate(b_type = "lrr_rp")

  con_lrr_rk <- x$lrr_rk %>%
    condense() %>%
    #dplyr::rename(kinase_pfams_hit = pfams_hit, kinase_pfams_acc = pfams_acc, kinase_pfams_loc = pfams_loc) %>%
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
         dplyr::distinct() %>%
         dplyr::group_by(.data$Name) %>%
         dplyr::summarise(
           sp_cut_site = .data$cut_site,
           .data$tm_start,
           .data$tm_end,
           pfams_hit = paste0(.data$hit, collapse=";"),
           pfams_acc = paste0(.data$acc, collapse=";"),
           pfams_loc = paste0(.data$pfam_coord, collapse=";"),
           ectos_hit = paste0(.data$QueryID, collapse=";"),
           ectos_coord = paste0(.data$hit_coord, collapse=";")

       ) %>%
       dplyr::distinct( ) %>%
       dplyr::ungroup()  %>%
       dplyr::mutate(b_type = "lrr_rp_rk_with_ecto")

    dplyr::bind_rows(
      list(
        con_lrr_rp, con_lrr_rk, con_non_lrr_rp, con_lrr_rk, con_lrr_rp_rk_ecto
      )
    ) %>%
      dplyr::distinct()
}
