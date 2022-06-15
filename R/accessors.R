#' get the dataframe of PFAM results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
hmmer_results <- function(b) {
  b$hmmer
}

#' get the dataframe of deeptmhmm results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
deeptmhmm_results <- function(b) {
  b$deeptmhmm
}
#' get the dataframe of ecto domain results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
ecto_results <- function(b) {
  b$ecto
}

#' get the dataframe of lrr_rp results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
lrr_rp <- function(b) {
  b$lrr_rp
}

#' get the dataframe of lrr_rk results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
lrr_rk <- function(b) {
  b$lrr_rk
}

#' get the dataframe of non_lrr_rp results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
non_lrr_rp <- function(b) {
  b$non_lrr_rp
}

#' get the dataframe of non_lrr_rk results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
non_lrr_rk <- function(b) {
  b$non_lrr_rk
}

#' get the dataframe of lrr_rp_rk_with_ecto results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
lrr_rp_rk_with_ecto <- function(b) {
  b$lrr_rp_rk_with_ecto
}



