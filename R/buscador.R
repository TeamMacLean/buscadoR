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


#' @export
buscar <- function(file_path, progress=FALSE, email = NULL, pfam_eval_cutoff=20, blast_eval_cutoff=1e-6, ...){

  load("data/pfam_data.rda")
  searches <- do_searches(file_path, progress, email,
                          pfam_eval=pfam_eval_cutoff, blast_eval=blast_eval_cutoff, ...)


  searches$file_path <- file_path


  ## add pfam info to all passing signal and tm domain carrying proteins
  phobius_pfam <- dplyr::left_join(searches$phobius, searches$pfam, by = c("Name" = "seq_name")) %>%
    dplyr::distinct()


  ## get the lrr_rp class
  searches$lrr_rp <- phobius_pfam %>%
    dplyr::filter(b_type == "LRR_PFAM", tm_start > seq_to)


  # get the lrr_rk class
   lrr_rk <- phobius_pfam %>%
     dplyr::filter(Name %in% searches$lrr_rp$Name,
                   b_type == "KINASE_PFAM",
                   pfam_length >= 250, tm_start < seq_from) %>%
     dplyr::distinct()
   rows_to_add <- dplyr::filter(searches$lrr_rp, Name %in% lrr_rk$Name )
    searches$lrr_rk <- dplyr::bind_rows(rows_to_add, lrr_rk)

  # remove the lrr_rps that became the lrr_rks from the lrr_rp object
   searches$lrr_rp <- searches$lrr_rp %>%
    dplyr::filter(!Name %in% searches$lrr_rk$Name)


  # get the non_lrr_rp class
   searches$non_lrr_rp <- phobius_pfam %>%
     dplyr::filter(! Name %in% c(searches$lrr_rp$Name, searches$lrr_rk),
                   b_type == "NON_LRR_PFAM", seq_from < tm_start ) %>%
     dplyr::distinct() #%>%

  # get the non_lrr_rk class
   searches$non_lrr_rk <- phobius_pfam %>%
     dplyr::filter(Name %in% searches$non_lrr_rp$Name, b_type == "KINASE_PFAM",
                   pfam_length >= 250, tm_start < seq_from) %>%
     dplyr::distinct()

    ## remove the non_lrr_rp that became non_lrr_rk
   searches$non_lrr_rp <- searches$non_lrr_rp %>%
     dplyr::filter(! Name %in% searches$non_lrr_rk)

  # get the ecto domain class
   phobius_pfam_ecto <- dplyr::left_join(phobius_pfam, searches$ecto, by = c("Name" = "SubjectID") ) %>%
     dplyr::distinct()

   searches$lrr_rp_rk_with_ecto <- phobius_pfam_ecto %>%
     dplyr::filter( Name %in% c(searches$lrr_rp$Name, searches$lrr_rk),
                    Perc.Ident > 50, S.start > cut_site, S.start < tm_start
                    ) %>%
     tidyr::unite(pfam_coord, seq_from:seq_to, sep="-", remove=FALSE)

  class(searches) <- "busco"
  return(searches)
}


condense <- function(df) {
  df %>%
    tidyr::unite(pfam_coord, seq_from:seq_to, sep="-") %>%
    dplyr::distinct() %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(
      sp_cut_site = cut_site,
      tm_start,
      tm_end,
      pfams_hit = paste0(hit, collapse=";"),
      pfams_acc = paste0(acc, collapse=";"),
      pfams_loc = paste0(pfam_coord, collapse=";")
    ) %>%
    dplyr::distinct( )%>%
    dplyr::ungroup()

}

#' @export
mesa <- function(df){

  if ("busco" %in% class(df)){
    df <- as.data.frame(df)
  }

  if (! "b_type" %in% colnames(df) ){
    stop("column b_type not found in data, cannot make table")
  }

  dplyr::group_by(df, b_type ) %>%
    dplyr::summarise(count = dplyr::n() ) %>%
    dplyr::ungroup() %>%
    knitr::kable(caption = "BuscadoR RLK Finding results")

}

do_searches <- function(file_path, progress=FALSE, email=NULL, pfam_eval=20, blast_eval=1e-6){

  result <- list(
    phobius = NULL,
    pfam = NULL,
    ecto = NULL
  )

  if (progress) message("Starting phobius search")

  result$phobius <- get_phobius(file_path, progress=progress) %>%
    dplyr::filter(is.phobius == TRUE, tm == 1) %>%
    dplyr::mutate(
      tm_start = as.numeric( stringr::str_split(
        stringr::str_split(prediction, "[oi]", simplify = TRUE)[,2],
        "-", simplify=TRUE)[,1] ) ,
      tm_end = as.numeric( stringr::str_split(
        stringr::str_split(prediction, "[oi]", simplify = TRUE)[,2],
        "-", simplify=TRUE)[,2] ),
    ) %>%
    dplyr::select(Name, cut_site, tm_start, tm_end) %>%
    dplyr::distinct()

  if (progress) message(paste("Phobius complete, found", length(result$phobius$Name), "proteins with Signal Peptide and single TM domain\n"))

  if (progress) message("Starting PFAM with found proteins")

  filtered_protein_tmpfile <- keep_phobius_hits(result$phobius, file_path)
  result$pfam <- get_pfam(filtered_protein_tmpfile, email, progress=progress, eval=pfam_eval) %>%
    dplyr::filter(eval < pfam_eval) %>%
    dplyr::mutate(base_acc = substr(acc, 1,7 ),
                  seq_to = as.numeric(seq_to),
                  seq_from = as.numeric(seq_from),
                  b_type = dplyr::if_else(base_acc %in% lrr_pfams, "LRR_PFAM",
                                          dplyr::if_else(base_acc %in% non_lrr_pfams, "NON_LRR_PFAM",
                                                         dplyr::if_else(base_acc %in% kinase_pfams, "KINASE_PFAM", "Other"))),
                  pfam_length = seq_to - seq_from
    ) %>%
    dplyr::distinct()

  if (progress) message("PFAM complete\nStarting BLAST for ectodomains with found proteins\n")
  result$ecto <- get_ecto(file_path, progress=progress) %>%
    dplyr::filter(E < blast_eval) %>%
    tidyr::unite(hit_coord, S.start:S.end, sep="-", remove=FALSE)

  return(result)
}


keep_phobius_hits <- function(df, file_path) {
  p <- Biostrings::readAAStringSet(file_path)
  p <- p[df$Name]
  tmpfile <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(p, tmpfile)
  return(tmpfile)
}



#' @export
pfam_results <- function(b) {
  b$pfam
}

#' @export
phobius_results <- function(b) {
  b$phobius
}
#' @export
ecto_results <- function(b) {
  b$ecto
}

#' @export
as.drawProteins <- function(b, which="lrr_rp") {

  seq_df <- seq_to_df(b) %>%
    dplyr::filter(seq_id %in% b[[which]]$Name) %>%
    dplyr::rename("entryName" = "seq_id", "length"="seq_length" ) %>%
    dplyr::select(entryName, length) %>%
    dplyr::mutate(accession=entryName, taxid=2712, order=1, type="CHAIN", begin=1, end=length, description=which) %>%
    dplyr::select(type, description, begin, end, length, accession, entryName, taxid, order)

  orders <- 1:length(seq_df$entryName)
  names(orders) <- seq_df$entryName

  tm_df <- b[[which]] %>%
    dplyr::rename("entryName"="Name", "begin"="tm_start", "end"="tm_end") %>%
    dplyr::distinct() %>%
    dplyr::mutate(accession=entryName, taxid=2712, order=1, type="DOMAIN", length=end-begin, description="transmembrane domain") %>%
    dplyr::select(type, description, begin, end, length, accession, entryName, taxid, order)

  domain_df <- b[[which]] %>%
    dplyr::rename("entryName"="Name", "begin"="seq_from", "end"="seq_to", "description"="hit") %>%
    dplyr::mutate(accession=entryName, taxid=2712, order=1, type="DOMAIN", length=end-begin)  %>%
    dplyr::select(type, description, begin, end, length, accession, entryName, taxid, order)

  return(dplyr::bind_rows(seq_df, tm_df, domain_df) %>%
         dplyr::mutate(order = orders[entryName])
         )

}

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

seq_to_df <- function(b){

  seqs <- Biostrings::readAAStringSet(filepath = b$file_path)
  data.frame(
    seq_id = names(seqs),
    sequence =  as.character(seqs, use.names=FALSE),
    seq_length = Biostrings::width(seqs)
  )
}
#' @export
write_seqs <- function(b, out_file_path) {

  ninfo <- data.frame(
    seq_id = as.data.frame(b)$Name,
    b_type = as.data.frame(b)$b_type
  )
  seq_df <- seq_to_df(b)
  seq_df <- dplyr::left_join(ninfo, seq_df, by="seq_id") %>%
    dplyr::mutate(seq_id = paste0(seq_id, "|", b_type))

  seqv <- seq_df$sequence
  names(seqv) <- seq_df$seq_id
  sset <- Biostrings::AAStringSet(seqv)
  Biostrings::writeXStringSet(sset, out_file_path)


}
#' @export
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
         dplyr::group_by(Name) %>%
         dplyr::summarise(
           sp_cut_site = cut_site,
           tm_start,
           tm_end,
           pfams_hit = paste0(hit, collapse=";"),
           pfams_acc = paste0(acc, collapse=";"),
           pfams_loc = paste0(pfam_coord, collapse=";"),
           ectos_hit = paste0(QueryID, collapse=";"),
           ectos_coord = paste0(hit_coord, collapse=";")

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
