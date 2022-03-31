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

#' create the buscador object
new_buscador <- function(...){


  x <- list(...)
  x['phobius'] = list(NULL)
  x['pfam'] = list(NULL)
  x['pfam_progress'] = list(NULL)
  x['ecto'] = list(NULL)
  x['lrr_rp'] = list(NULL)
  x['lrr_rk'] = list(NULL)
  x['lrr_rp_rk_with_ecto'] = list(NULL)
  x['non_lrr_rp'] = list(NULL)
  x['non_lrr_rk'] = list(NULL)
  x[['aastringset']] = Biostrings::readAAStringSet(x$protein_file)


  b <- structure(x,
  class = "buscador"
  )
}
#' perform internet searches and classify given proteins into receptor types
#'
#' runs the internet search using the Phobius server, the PFAM API and local BLAST. Then collates the results
#' and classifies proteins into receptor types as defined in `definiciones()`
#'
#' @param protein_file path to protein fasta file (ignored if restart file has some searches). Do not use if you want to restart
#' @param restart_file file to use as restart file, will perform the search actions and automatically save each step as you go. if file exists will be appended with new search results. search actions already performed and saved in the file will not be replaced
#' @param progress show a progress bar as job runs (not recommended for non interactive sessions)
#' @param email a valid email address
#' @param pfam_eval_cutoff exclude pfam hits with eval over this value (20)
#' @param blast_eval_cutoff exclude blast hits with eval over this value (1e-6)
#' @param wait seconds to wait between PFAM checks (5)
#' @param maxchecktime seconds to stop checking PFAM and abort (120)
#' @return buscador object list of tidy dataframes of each type of receptor and its Phobius, PFAM and BLAST hits
#' @export
#' @importFrom rlang .data
buscar <- function(protein_file=NULL, restart_file=NULL, progress=FALSE, email = NULL, pfam_eval_cutoff=20, blast_eval_cutoff=1e-6, wait=1, maxchecktime=1){

  busc <- NULL

  if ( (!is.null(protein_file) && file.exists(protein_file))  & !file.exists(restart_file)){ #we are starting from scratch
   message(sprintf("restart_file %s doesn't exist, interpreting this as new job with %s as input", restart_file, protein_file ))

    busc <- new_buscador(
      protein_file = protein_file,
      restart_file = restart_file,
      progress = progress,
      email = email,
      pfam_eval_cutoff = pfam_eval_cutoff,
      blast_eval_cutoff = blast_eval_cutoff,
      wait = wait,
      maxchecktime = maxchecktime
    )


  } else if ( file.exists(restart_file) ){ # we have a restart file so need to pick up
    message(sprintf("restart_file %s found, interpreting this as a restarted job with existing protein set.", restart_file))
    busc <- readRDS(restart_file)
  } else {
    stop("Couldn't find restart_file or protein_file.")
  }

  busc <- do_searches(busc)

  if (is.null(busc$phobius)) warning(sprintf("No phobius result in buscador object.")); return(busc)
  if (attr(busc$pfam, "status") == "submission" ) warning(sprintf("PFAMScan job not yet complete on ebi.ac.uk  Please try a restart later.")); return(busc)
  if (is.null(busc$ecto)) warning(sprintf("No ectodomain result in buscador object.")); return(busc)


  ## add pfam info to all passing signal and tm domain carrying proteins
  phobius_pfam <- dplyr::left_join(busc$phobius, busc$pfam, by = c("Name" = "seq_name")) %>%
    dplyr::distinct()


  ## get the lrr_rp class
  busc$lrr_rp <- phobius_pfam %>%
    dplyr::filter(.data$b_type == "LRR_PFAM", .data$tm_start > .data$seq_to)


  # get the lrr_rk class
   lrr_rk <- phobius_pfam %>%
     dplyr::filter(.data$Name %in% busc$lrr_rp$Name,
                   .data$b_type == "KINASE_PFAM",
                   .data$pfam_length >= 250, .data$tm_start < .data$seq_from) %>%
     dplyr::distinct()
   rows_to_add <- dplyr::filter(busc$lrr_rp, .data$Name %in% lrr_rk$Name )
    busc$lrr_rk <- dplyr::bind_rows(rows_to_add, lrr_rk)

  # remove the lrr_rps that became the lrr_rks from the lrr_rp object
   busc$lrr_rp <- busc$lrr_rp %>%
    dplyr::filter(!.data$Name %in% busc$lrr_rk$Name)


  # get the non_lrr_rp class
   busc$non_lrr_rp <- phobius_pfam %>%
     dplyr::filter(! .data$Name %in% c(busc$lrr_rp$Name, busc$lrr_rk),
                   .data$b_type == "NON_LRR_PFAM", .data$seq_from < .data$tm_start ) %>%
     dplyr::distinct() #%>%

  # get the non_lrr_rk class
   busc$non_lrr_rk <- phobius_pfam %>%
     dplyr::filter(.data$Name %in% busc$non_lrr_rp$Name, .data$b_type == "KINASE_PFAM",
                   .data$pfam_length >= 250, .data$tm_start < .data$seq_from) %>%
     dplyr::distinct()

    ## remove the non_lrr_rp that became non_lrr_rk
   busc$non_lrr_rp <- busc$non_lrr_rp %>%
     dplyr::filter(! .data$Name %in% busc$non_lrr_rk)

  # get the ecto domain class
   phobius_pfam_ecto <- dplyr::left_join(phobius_pfam, busc$ecto, by = c("Name" = "SubjectID") ) %>%
     dplyr::distinct()

   busc$lrr_rp_rk_with_ecto <- phobius_pfam_ecto %>%
     dplyr::filter( .data$Name %in% c(busc$lrr_rp$Name, busc$lrr_rk),
                    .data$Perc.Ident > 50, .data$S.start > .data$cut_site, .data$S.start < .data$tm_start
                    ) %>%
     tidyr::unite(pfam_coord, .data$seq_from:.data$seq_to, sep="-", remove=FALSE)

  return(busc)
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

  if ("buscador" %in% class(b)){
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

#' convert aastrset to a tempfile on disk
#' @param aastrset Biostrings::aastrset object
#' @return filename
totempfile <- function(aastrset) {
  tmpfile <- tempfile(pattern = "buscador_", tmpdir = tempdir(), fileext = ".fa")
  Biostrings::writeXStringSet(aastrset, tmpfile)
  return(tmpfile)
}

#' convert aastrset to a fasta character vector
#' @param aastrset Biostrings::aastrset object
#' @return character vector with fasta strings
tostrvec <- function(aastrset) {
  paste0(">", names(as.character(aastrset)), "\n", as.character(aastrset))
}

#' run searches over the internet and local BLAST.
#'
#' Runs searches for proteins at Phobius, PFAM and local BLAST.
#' Reformat returned results from searches and apply basic quality filters
#'
#' if the search at PFAM does not complete the buscador$pfam_progress object is
#' populated with submission IDs for a later restart run and the restart_file written.
#' In later runs the PFAM server is queried to collect all completed runs. This process should
#' be repeated until all the queries are complete. Once all queries are complete the buscador$pfam object gets the
#' dataframe of all PFAM results and is marked with class "complete"
#'
#' @param busc buscador object from `buscar`
#' @return search populated buscador object
#' @importFrom rlang .data
do_searches <- function(busc){


  if( is.null(busc$phobius)){
    if (busc$progress) message("Starting phobius search")

      busc$phobius = get_phobius(totempfile(busc$aastringset), progress=busc$progress) %>% ##set file
        process_phobius()
      busc$aastringset = busc$aastringset[busc$phobius$Name] #reduce to only proteins passing Phobius
      saveRDS(busc, busc$restart_file)
      if (length(busc$phobius$Name) == 0) stop("No proteins with Signal Peptide found, quitting.")

    if (busc$progress) message(paste("Phobius complete, found", length(busc$phobius$Name), "proteins with Signal Peptide and single TM domain\n"))
  } else {
    if (busc$progress) message(paste("Phobius result found in restart_file, not redoing"))
  }

  if ( is.null(busc$pfam)){
    if (busc$progress) message("Starting PFAMScan jobs at ebi.ac.uk.")

    busc$pfam = submit_pfam(tostrvec(busc$aastringset), email=busc$email, progress=busc$progress, eval=busc$pfam_eval, wait=busc$wait, maxchecktime=busc$maxchecktime)

    if (attr(busc$pfam, "status") == "submission"){ ## if initial submission didnt complete in time and a restart will be needed
      busc$pfam_progress = lapply(busc$pfam$id, function(x) { NULL })
      names(busc$pfam_progress) <- busc$pfam$id

    } else if (attr(busc$pfam, "status") == "complete"){

      busc$pfam = busc$pfam %>%
        process_pfam(busc$pfam_eval)
      busc$pfam_progess = NULL
      attr(busc$pfam, "status") <- "complete"

    }

    saveRDS(busc, busc$restart_file)

  } else if (attr(busc$pfam, "status") == "submission") {
    if (busc$progress) message("Checking PFAMScan at ebi.ac.uk for status of previously submitted jobs")

    busc <- retrieve_pfam(busc)

    if (attr(busc$pfam, "status") == "complete"){
      busc$pfam = busc$pfam %>%
        process_pfam(busc$pfam_eval)
      busc$pfam_progess = NULL
      attr(busc$pfam, "status") <- "complete"
    }

    saveRDS(busc, busc$restart_file)

  } else if (attr(busc$pfam, "status") == "complete") {
    if (busc$progress) message("PFAMScan result found in restart_file, not redoing")
  }

  if (is.null(busc$ecto)){
    if (busc$progress) message("Starting ectodomain BLAST on local machine.")
    busc$ecto = get_ecto(totempfile(busc$aastringset), progress=busc$progress) %>%
    dplyr::filter(.data$E < busc$blast_eval) %>%
    tidyr::unite(hit_coord, .data$S.start:.data$S.end, sep="-", remove=FALSE)

    saveRDS(busc, busc$restart_file)

  } else {
    if (busc$progress) message("Ectodomain result found in restart_file, not redoing")
  }
  return(busc)
}


#' check the status of the buscador object
#'
#' are all the searches done and the results compiled
#'
#' @param busc a buscador object
#' @export
completed <- function(busc) {
  data.frame(
   result = c(
     "phobius_search",
     "pfam_retrieval",
     "ectodomain_search",
     "lrr_rp_annotation",
     "lrr_rk_annotation",
     "lrr_rp_rk_with_ecto_annotation",
     "non_lrr_rp_annotation",
     "non_lrr_rk_annotation"
   ),
   completed = c(
    !is.null(busc$phobius),
    attr(busc$pfam, "status") == "complete",
    !is.null(busc$ecto),
    !is.null(busc$lrr_rp),
    !is.null(busc$lrr_rk),
    !is.null(busc$lrr_rp_rk_with_ecto),
    !is.null(busc$non_lrr_rp),
    !is.null(busc$non_lrr_rk)
  )
  )
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
#' @param b buscador search object returned from `buscar`
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

#' convert the aastrset to a dataframe
#' @param b `buscador` search object returned from `buscar()`
seq_to_df <- function(b){

  data.frame(
    seq_id = names(b$aastringset),
    sequence =  as.character(b$aastringset, use.names=FALSE),
    seq_length = Biostrings::width(b$aastringset)
  )
}

#' write the annotated receptor sequences to a fasta file
#'
#' @param b `buscador` search object returned from `buscar()`
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

#' generic for converting `buscador` search object to a dataframe
#' @param x `buscador` object to coerce
#' @param ... parameters for other functions
#' @export
#' @importFrom rlang .data
as.data.frame.buscador <- function(x,...){

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
