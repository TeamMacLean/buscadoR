#' @export
definiciones <- function(which="all") {

  types <- c('all', 'lrr_rp', 'lrr_rk', 'non_lrr_rp', 'non_lrr_rk', 'lrr_rp_ecto', 'lrr_rk_ecto')

  if (! which %in% types){
    stop('Unknown definition requested\nChoose one of:\n', paste(types, collapse=",", sep=","))
  }

  defs <- list(

    lrr_rp = "lrr_rp =
     * Exactly one Transmembrane Domain according to Phobius.
     * A signal peptide according to Phobius.
     * At least one of the `lrr_pfams` according to PFAMscan.
     * The Transmembrane domain should be closer to the C terminal than the end of the pfam hit.",

    lrr_rk = "\nlrr_rk =
     * An llr_rp.
     * At least one `kinase_pfam` according to PFAMscan.
     * The `kinase_pfam` should be at least 250 AAs long.
     * The `kinase_pfam` closer to the C terminal than the TM domain.",

    non_lrr_rp = "\nnon_lrr_rp =
     * Not an lrr_rp or an lrr_rk.
     * At least one `non_lrr_pfam` according to PFAMscan.
     * Closer to the N terminal than the start of the Transmembrane domain.",

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
buscar <- function(file_path, quiet=TRUE, email = NULL, pfam_eval_cutoff=1e-6, blast_cut_off=1e-6, ...){

  searches <- do_searches(file_path, quiet, email, ...)
  load("data/pfam_data.rda")


  phobius <- searches$phobius %>% #find proteins with 1 TM domain and signal peptide
    dplyr::filter(is.phobius == TRUE,
                  tm == 1) %>%
    dplyr::mutate(
      tm_start = stringr::str_split(
                  stringr::str_split(prediction, "[oi]", simplify = TRUE)[,2],
                  "-", simplify=TRUE)[,1],
      tm_end = stringr::str_split(
        stringr::str_split(prediction, "[oi]", simplify = TRUE)[,2],
        "-", simplify=TRUE)[,2],

      ) %>%

    dplyr::select(Name, cut_site, tm_start, tm_end) %>%
    dplyr::distinct()

  searches$prc_phob <- phobius


  ## find the pfams
  pfam <- searches$pfam %>%
    dplyr::filter(eval < pfam_eval_cutoff) %>%
    dplyr::mutate(base_acc = substr(acc, 1,7 ),
                  seq_to = as.numeric(seq_to),
                  seq_from = as.numeric(seq_from),
                  b_type = dplyr::if_else(base_acc %in% lrr_pfams, "LRR_PFAM",
                           dplyr::if_else(base_acc %in% non_lrr_pfams, "NON_LRR_PFAM",
                           dplyr::if_else(base_acc %in% kinase_pfams, "KINASE_PFAM", "Other"))),
                  pfam_length = seq_to - seq_from
                  ) %>%
    dplyr::distinct()



  ## add pfam info to all passing signal and tm domain carrying proteins
  phobius_pfam <- dplyr::left_join(phobius, pfam, by = c("Name" = "seq_name")) %>%
    dplyr::distinct()

  ## find the ectos
  ecto <- searches$ecto %>%
    dplyr::filter(E < blast_cut_off) %>%
    tidyr::unite(hit_coord, S.start:S.end, sep="-", remove=FALSE)


  ## get the lrr_rp class
  searches$lrr_rp <- phobius_pfam %>%
    dplyr::filter(b_type == "LRR_PFAM", tm_start > seq_to) %>%
    condense() %>%
    dplyr::rename(lrr_pfams_hit = pfams_hit, lrr_pfams_acc = pfams_acc, lrr_pfams_loc = pfams_loc)


  # get the lrr_rk class
   searches$lrr_rp_with_rk <- phobius_pfam %>%
     dplyr::filter(Name %in% searches$lrr_rp$Name, b_type == "KINASE_PFAM", pfam_length >= 250, tm_start < seq_from) %>%
     dplyr::distinct() %>%
     condense() %>%
     dplyr::rename(kinase_pfams_hit = pfams_hit, kinase_pfams_acc = pfams_acc, kinase_pfams_loc = pfams_loc)

  # remove the lrr_rps that became the lrr_rks from the lrr_rp object
  searches$lrr_rp <- searches$lrr_rp %>%
    dplyr::filter(! Name %in% searches$lrr_rp_with_rk)


  # get the non_lrr_rp class
   searches$non_lrr_rp <- phobius_pfam %>%
     dplyr::filter(! Name %in% c(searches$lrr_rp$Name, searches$lrr_rp_with_rk),
                   b_type == "NON_LRR_PFAM", seq_from < tm_start ) %>%
     dplyr::distinct() %>%
     condense() %>%
     dplyr::rename(non_lrr_rp_pfams_hit = pfams_hit, non_lrr_rp_pfams_acc = pfams_acc, non_lrr_rp_pfams_loc = pfams_loc)

  # get the non_lrr_rk class
   searches$non_lrr_rp_with_rk <- phobius_pfam %>%
     dplyr::filter(Name %in% searches$non_lrr_rp$Name, b_type == "KINASE_PFAM",
                   pfam_length >= 250, tm_start < seq_from) %>%
     dplyr::distinct() %>%
     condense() %>%
     dplyr::rename(non_lrr_rp_with_rk_pfams_hit = pfams_hit, non_lrr_rp__with_rk_pfams_acc = pfams_acc, non_lrr_rp__with_pfams_loc = pfams_loc)

    ## remove the non_lrr_rp that became non_lrr_rk
   searches$non_lrr_rp <- searches$non_lrr_rp %>%
     dplyr::filter(! Name %in% searches$non_lrr_rp_with_rk)

  # get the ecto domain class
   phobius_pfam_ecto <- dplyr::left_join(phobius_pfam, ecto, by = c("Name" = "SubjectID") ) %>%
     dplyr::distinct()

   searches$lrr_rp_rk_with_ecto <- phobius_pfam_ecto %>%
     dplyr::filter( Name %in% c(searches$lrr_rp$Name, searches$lrr_rp_with_rk),
                    Perc.Ident > 50, S.start > cut_site, S.start < tm_start
                    ) %>%
     tidyr::unite(pfam_coord, seq_from:seq_to, sep="-") %>%
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
     dplyr::distinct( )%>%
     dplyr::ungroup()


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
  if (! "b_type" %in% colnames(df) ){
    stop("column b_type not found in data, cannot make table")
  }

  table(df$b_type)

}

do_searches <- function(file_path, quiet=TRUE, email=NULL){

  result <- list(
    phobius = NULL,
    pfam = NULL,
    ecto = NULL
  )
  if (quiet){
    result$phobius = get_phobius(file_path, progress=FALSE)
    result$pfam = get_pfam(file_path, email, progress=FALSE)
    result$ecto = get_ecto(file_path, progress=FALSE)
  } else {
    message("Starting phobius search\n")
    result$phobius = get_phobius(file_path, progress=TRUE)
    message("phobius search completed\n")
    message("Starting PFAMscan search\n")
    result$pfam = get_pfam(file_path, email, progress=TRUE)
    message("PFAMscan search completed\n")
    message("Starting BLAST for ectodomains\n")
    result$ecto = get_ecto(file_path, progress=TRUE)
    message("BLAST for ectodomains completed\n")
  }
  return(result)
}
